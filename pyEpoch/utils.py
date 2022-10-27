import numpy as np
import pandas as pd
import scanpy as sc
import skfda
from pygam import GAM, s,l
from skfda import FDataGrid
import skfda.preprocessing.smoothing.kernel_smoothers as ks
import igraph as igraph
import scipy




def makeExpMat(adata):
    expMat = pd.DataFrame(adata.X.T, index = adata.var_names, columns = adata.obs_names).copy()
    return expMat




def makeSampTab(adata):
    sampTab = adata.obs.copy()
    return sampTab



def gamFit(expMat,genes,celltime):

    genes2=(set(genes) & set(expMat.index))
    def abcd(input_data):
        z=pd.DataFrame()
        z["z"]=input_data.values
        z["t"]=celltime.values
        z.index=expMat.columns
        X=celltime.values.reshape((celltime.shape[0],1))
        y=z["z"].values

        gam=GAM(l(0)).fit(X,y)
        p=gam.statistics_['p_values'][0]
        return p
    ans=expMat.loc[genes2][celltime.index].apply(abcd,axis=1)
    return ans




def grnKsmooth(adata,BW=.25,pseudotime_column=None):
    #cells=ccells
    #BW=.1
    
    newadata = adata.copy()

    if 'cells' in newadata.uns.keys():
        cellids = newadata.uns['cells'].index.tolist()
        newadata = newadata[cellids,:]
    else:
        cells=pd.DataFrame()
        cells["cell_name"]=newadata.obs.index
        cells["pseudotime"]=newadata.obs[pseudotime_column].values
        cells.index=cells["cell_name"]
        cells=cells.sort_values(by="pseudotime")
        newadata.uns['cells'] = cells

    if 'dgenes' in newadata.uns.keys():
        newadata = newadata[:,newadata.uns['dgenes']]


    expDat = makeExpMat(newadata)
    cells = newadata.uns['cells']

    # genes=adata.var.index
    # expDat=pd.DataFrame(adata.X).T
    # expDat.columns=adata.obs.index
    # expDat.index=genes
    # expDat=expDat.loc[expDat.sum(axis=1)!=0]
    
    BW=min(BW, max(cells["pseudotime"])-min(cells["pseudotime"])/10)
    t1=pd.DataFrame(cells["pseudotime"])
    t1.index=cells["cell_name"]
    t1=t1.sort_values(by='pseudotime', ascending=True)
    #expDat.iloc[t1.index]

    # order expDat
    expDat=expDat[list(t1.index)]

    ans=pd.DataFrame(columns=np.arange(expDat.shape[1]))
    for i in np.arange(expDat.shape[0]):
        y=expDat.iloc[i].values
        x=t1["pseudotime"].values
        fd = FDataGrid(sample_points=[x],data_matrix=[y])
        smoother = ks.NadarayaWatsonSmoother(smoothing_parameter=BW)
        smoothed = smoother.fit_transform(fd)
        a=smoothed.data_matrix.round(10)
        each_row=[]
        for j in a:
            for k in j:
                for l in k:
                    each_row.append(l)
        ans=pd.concat([ans,pd.DataFrame(each_row).T])

    ans.index=expDat.index
    ans.columns=expDat.columns

    adata.uns['smoothed_expression'] = ans

    print("Done. Smoothed expression stored in .uns['smoothed_expression'].")

    return adata


# def add_interactions_types(adata, grn="dynamic_GRN"):
#     if type(adata.uns[grn])==pd.core.frame.DataFrame:
#         GRN=adata.uns[grn]
#         if "corr" not in GRN.columns:
#             sys.exit("Missing 'corr' information. Run reconstruction first.")
#         else:
#             true_false=GRN["corr"]>0
#             activation=[]
#             for i in true_false:
#                 if i==True:
#                     activation.append("activation")
#                 else:
#                     activation.append("repression")
#             GRN["interaction"]=activation
#         return GRN

#     elif type(adata.uns[grn])==dict:
#         GRN=adata.uns[grn]
#         for i in GRN:
#             if "corr" not in GRN[i].columns:
#                 sys.exit("Missing 'corr' information. Run reconstruction first.")
#             else:
#                 true_false=GRN[i]["corr"]>0
#                 activation=[]
#                 for j in true_false:
#                     if j==True:
#                         activation.append("activation")
#                     else:
#                         activation.append("repression")
#             GRN[i]["interaction"]=activation
#         return GRN




# New and updated... 
# If correlation is missing, compute and add it in
def add_interaction_type(adata, grn="dynamic_GRN"):
    if type(adata.uns[grn])==pd.core.frame.DataFrame:
        GRN=adata.uns[grn]
        if "corr" not in GRN.columns:
            print("No correlation stored, computing correlation.")
            GRN = add_corr(GRN,adata)

        true_false=GRN["corr"]>0
        activation=[]
        for i in true_false:
            if i==True:
                activation.append("activation")
            else:
                activation.append("repression")
        GRN["interaction"]=activation
        adata.uns[grn] = GRN

    elif type(adata.uns[grn])==dict:
        GRN=adata.uns[grn]
        for i in GRN:
            if "corr" not in GRN[i].columns:
                print("No correlation stored, computing correlation.")
                GRN[i] = add_corr(GRN[i],adata)

            true_false=GRN[i]["corr"]>0
            activation=[]
            for j in true_false:
                if j==True:
                    activation.append("activation")
                else:
                    activation.append("repression")
            GRN[i]["interaction"]=activation
        adata.uns[grn] = GRN

    return adata



# a little side function...
def add_corr(grnDF,adata):

    if scipy.sparse.issparse(adata.X):
        expX = pd.DataFrame(adata.X.todense())
    else:
        expX=pd.DataFrame(adata.X)
    
    expX.index=adata.obs_names
    expX.columns=adata.var_names

    corrmat=expX.corr()
    corrmat['TG'] = corrmat.index
    corrdf=pd.melt(corrmat, id_vars=["TG"])
    corrdf.columns=["TG","TF","corr"]

    newdf = grnDF.merge(corrdf,on=['TG','TF'],how="left")

    return newdf





def find_communities(adata,grn="dynamic_GRN",use_weights=False,weight_column=None,communities_slot=None):
    #use_weights=True
    #weight_column="weighted_score"

    if type(adata.uns[grn])==pd.core.frame.DataFrame:

        if communities_slot is None:
            communities_slot = "static_communities"

        GRN=adata.uns[grn]

        if use_weights==True:
            weight_value=GRN[weight_column]
        else:
            weight_value=None

        tuples = [tuple(x) for x in GRN.values]
        Gm = igraph.Graph.TupleList(tuples, directed = False)
        
        genes_ordered=[]
        for i in Gm.vs:
            genes_ordered.append(i["name"])

        communities=Gm.community_multilevel(weights=weight_value)

        communities_names=0
        communities_genes=[]
        communities_names_list=[]
        for i in communities:
            communities_names=communities_names+1
            for j in i:
                communities_genes.append(genes_ordered[j])
                communities_names_list.append(j)
        communities_dataframe=pd.DataFrame()
        communities_dataframe["genes"]=communities_genes
        communities_dataframe["communities"]=communities_list

        adata.uns[communities_slot]=communities_dataframe
        print("Communities saved in .uns[" + communities_slot + "].")

        return adata

    elif type(adata.uns[grn])==dict:

        if communities_slot is None:
            communities_slot = "dynamic_communities"

        GRN=adata.uns[grn]

        communities_epoch_dataframes=[]
        for i in GRN:
            
            if use_weights==True:
                weight_value=GRN[i][weight_column]
            else:
                weight_value=None
            

            tuples = [tuple(x) for x in GRN[i].values]
            Gm = igraph.Graph.TupleList(tuples, directed = False)
            genes_ordered=[]
            for j in Gm.vs:
                genes_ordered.append(j["name"])

            communities=Gm.community_multilevel(weights=weight_value)

            communities_names=0
            communities_genes=[]
            communities_names_list=[]
            for k in communities:
                communities_names=communities_names+1
                for j in k:
                    communities_genes.append(genes_ordered[j])
                    communities_names_list.append(communities_names)
            communities_dataframe=pd.DataFrame()
            communities_dataframe["genes"]=communities_genes
            communities_dataframe["communities"]=communities_names_list
            
            communities_epoch_dataframes.append(communities_dataframe)

        l1=GRN.keys()
        res=dict(zip(l1,communities_epoch_dataframes))
        
        adata.uns[communities_slot]=res
        print("Communities saved in .uns[" + communities_slot + "].")
    
        return adata





















