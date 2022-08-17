import numpy as np
import pandas as pd
import sys
import igraph as igraph
from .utils import *



def edge_uniqueness(grn_DFs,tfs,weight_column):
    #grn_DFs=[net1,net2]  
    #weight_column="weighted_score"
    #tfs=mmTFs

    # grn_DFs = 
    genes=[]
    for i in grn_DFs:
        values=list(i["TG"])+list(i["TF"])
        genes=genes+values
    genes=np.unique(genes)
    tfs=list(set(genes) & set(tfs))

    adj_list=[]
    for i in grn_DFs:
        tuples = [tuple(x) for x in i.values]
        Gm = igraph.Graph.TupleList(tuples, directed = True, edge_attrs = [weight_column])
        vs=igraph.VertexSeq(Gm)

        vertex=[]
        for j in Gm.vs:
            vertex.append(j["name"])

        addvtcs=list(set(genes) - set(vertex))
        Gm.add_vertices(len(addvtcs))
        Gm.vs["name"]=vertex+addvtcs

        adj=pd.DataFrame(Gm.get_adjacency(attribute=weight_column),index=Gm.vs["name"],columns=Gm.vs["name"])
        adj=adj.loc[tfs].T
        adj_list.append(adj)
    res=edge_rank(adj_list)
    return res



def edge_rank(grnMats):
    #grnMats=adj_list

    i=grnMats[0].sort_index(key=lambda x: x.str.lower())
    TFs=[]
    for j in i.columns:
        TFs=TFs+[j]*i.shape[0]
    TGs=[]
    for j in i.index:
        TGs=TGs+[j]*i.shape[1]
    full_df=pd.DataFrame()
    full_df["TGs"]=TGs
    full_df["TFs"]=TFs

    count=0
    for i in grnMats:
        i=i.sort_index(key=lambda x: x.str.lower())
        count=count+1
        full_df["network"+str(count)]=i.values.ravel()


    mins=[]
    maxes=[]
    for i in np.arange(full_df.shape[0]):
        mins.append(min(full_df.loc[i].values[2:]))
        maxes.append(max(full_df.loc[i].values[2:]))
    diff=list(np.array(maxes)-np.array(mins))
    full_df["min"]=mins
    full_df["max"]=maxes
    full_df["diff"]=diff

    full_df=full_df.sort_values(by='diff',ascending=False)
    return full_df


# def dynamic_difference_network(adata, adata1, condition, types, diff_thresh=3, condition_thresh=6):
#     #edgeDF=res
#     #epochs=[adata1,adata2]
#     #condition="network1"
#     #types="on"
#     #diff_thresh=7.5
#     #condition_thresh=10
    
#     edgeDF=[adata.uns["grnDF"], adata1.uns["grnDF"]]
#     edgeDF=edgeDF.loc[edgeDF["diff"]!=0]
    
#     epochs1 = adata.uns['epochs'].pop('mean_expression')
#     epochs2 = adata1.uns['epochs'].pop('mean_expression')
    
#     epochs_combined=epochs1

#     for i in epochs2:
#             epochs_combined[i]=set(epochs_combined[i]).union(set(i))
            
#     print(epochs_combined)
        
#     dynamic_edges=[]
#     for i in epochs_combined:
#         dynamic_edges.append(edgeDF.loc[edgeDF["TFs"].isin(epochs_combined[i])])

#     edgeDF=dynamic_edges
#     conditions=[]
#     for i in list(set(list(edgeDF[0].columns))-set(['TGs', 'TFs','min', 'max', 'diff'])):
#         conditions.append(i)  

#     if condition not in conditions:
#         sys.exit("Condition not in network.")

#     diffnet=[]
#     for i in edgeDF:
#         diffnet.append(i[["TGs","TFs",condition,"diff"]])

#     for i in np.arange(len(diffnet)):
#         diffnet[i]=diffnet[i].loc[diffnet[i]["diff"]>diff_thresh]
#     #print(diffnet)
    
#     if types=="on":
#         for i in np.arange(len(diffnet)):
#             diffnet[i]=diffnet[i].loc[diffnet[i][condition]>=condition_thresh]
#     elif types=="off":
#         for i in np.arange(len(diffnet)):
#             diffnet[i]=diffnet[i].loc[diffnet[i][condition]<condition_thresh]
#     else:
#         sys.exit("Types should be 'on' or 'off'")
#     #print(diffnet)
#     return diffnet

def dynamic_difference_network(edgeDF, adata, adata1, condition, types, diff_thresh=3, condition_thresh=6):
    #edgeDF=res
    #epochs=[adata1,adata2]
    #condition="network1"
    #types="on"
    #diff_thresh=7.5
    #condition_thresh=10
    
    edgeDF=edgeDF.loc[edgeDF["diff"]!=0]
    
    epochs1 = adata.uns['epochs']
    epochs2 = adata1.uns['epochs']
    
    epochs_combined=epochs1

    for i in epochs2:
            epochs_combined[i]=set(epochs_combined[i]).union(set(epochs2[i]))
                    
    dynamic_edges=[]
    for i in epochs_combined:
        dynamic_edges.append(edgeDF.loc[edgeDF["TFs"].isin(epochs_combined[i])])

    edgeDF=dynamic_edges
    conditions=[]
    for i in list(set(list(edgeDF[0].columns))-set(['TGs', 'TFs','min', 'max', 'diff'])):
        conditions.append(i)  

    if condition not in conditions:
        sys.exit("Condition not in network.")

    diffnet=[]
    for i in edgeDF:
        diffnet.append(i[["TGs","TFs",condition,"diff"]])

    for i in np.arange(len(diffnet)):
        diffnet[i]=diffnet[i].loc[diffnet[i]["diff"]>diff_thresh]
    #print(diffnet)
    
    if types=="on":
        for i in np.arange(len(diffnet)):
            diffnet[i]=diffnet[i].loc[diffnet[i][condition]>=condition_thresh]
    elif types=="off":
        for i in np.arange(len(diffnet)):
            diffnet[i]=diffnet[i].loc[diffnet[i][condition]<condition_thresh]
    else:
        sys.exit("Types should be 'on' or 'off'")
    #print(diffnet)
    return diffnet



# def dynamic_difference_network(edgeDF, epochs, condition, types, diff_thresh=3, condition_thresh=6):
#     #edgeDF=res
#     #epochs=[adata1,adata2]
#     #condition="network1"
#     #types="on"
#     #diff_thresh=7.5
#     #condition_thresh=10 

#     edgeDF=edgeDF.loc[edgeDF["diff"]!=0]


#     epochs_combined={}
#     for i in epochs:
#         for j in i.uns["epochs"]:
#             if j not in epochs_combined.keys():
#                 epochs_combined[j]=i.uns["epochs"][j]
#             else:
#                 epochs_combined[j]=epochs_combined[j]+i.uns["epochs"][j]
#     print(epochs_combined)
#     for i in epochs_combined:
#         epochs_combined[i]=np.unique(epochs_combined[i])
        
#     dynamic_edges=[]
#     for i in epochs_combined:
#         dynamic_edges.append(edgeDF.loc[edgeDF["TFs"].isin(epochs_combined[i])])

#     edgeDF=dynamic_edges
#     conditions=[]
#     for i in list(set(list(edgeDF[0].columns))-set(['TGs', 'TFs','min', 'max', 'diff'])):
#         conditions.append(i)  

#     if condition not in conditions:
#         sys.exit("Condition not in network.")

#     diffnet=[]
#     for i in edgeDF:
#         diffnet.append(i[["TGs","TFs",condition,"diff"]])

#     for i in np.arange(len(diffnet)):
#         diffnet[i]=diffnet[i].loc[diffnet[i]["diff"]>diff_thresh]
#     #print(diffnet)
    
#     if types=="on":
#         for i in np.arange(len(diffnet)):
#             diffnet[i]=diffnet[i].loc[diffnet[i][condition]>=condition_thresh]
#     elif types=="off":
#         for i in np.arange(len(diffnet)):
#             diffnet[i]=diffnet[i].loc[diffnet[i][condition]<condition_thresh]
#     else:
#         sys.exit("Types should be 'on' or 'off'")
#     #print(diffnet)
#     return diffnet




# def add_interactions_types(diffnet, types, grnDF_on, grnDF_offlist):
#     diffnet = pd.DataFrame(diffnet)
    



def add_interactions_types(adata):
    if type(adata.uns["dynamic_GRN"])==pd.core.frame.DataFrame:
        dynamic_GRN=adata.uns["dynamic_GRN"]
        if "corr" not in dynamic_GRN.columns:
            sys.exit("Missing 'corr' information. Run reconstruction first.")
        else:
            true_false=dynamic_GRN["corr"]>0
            activation=[]
            for i in true_false:
                if i==True:
                    activation.append("activation")
                else:
                    activation.append("repression")
            dynamic_GRN["interaction"]=activation
        return dynamic_GRN

    elif type(adata.uns["dynamic_GRN"])==dict:
        dynamic_GRN=adata.uns["dynamic_GRN"]
        for i in dynamic_GRN:
            if "corr" not in dynamic_GRN[i].columns:
                sys.exit("Missing 'corr' information. Run reconstruction first.")
            else:
                true_false=dynamic_GRN[i]["corr"]>0
                activation=[]
                for j in true_false:
                    if j==True:
                        activation.append("activation")
                    else:
                        activation.append("repression")
            dynamic_GRN[i]["interaction"]=activation
        return dynamic_GRN

    
def find_communities(grn,use_weights=False,weight_column=None):
    #use_weights=True
    #weight_column="weighted_score"
    if type(adata.uns["dynamic_GRN"])==pd.core.frame.DataFrame:
        dynamic_GRN=adata.uns["dynamic_GRN"]

        if use_weights==True:
            weight_value=dynamic_GRN[weight_column]
        else:
            weight_value=None

        tuples = [tuple(x) for x in dynamic_GRN.values]
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
        return communities_dataframe

    elif type(adata.uns["dynamic_GRN"])==dict:
        dynamic_GRN=adata.uns["dynamic_GRN"]

        communities_epoch_dataframes=[]
        for i in dynamic_GRN:
            
            if use_weights==True:
                weight_value=dynamic_GRN[i][weight_column]
            else:
                weight_value=None
            

            tuples = [tuple(x) for x in dynamic_GRN[i].values]
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
        return communities_epoch_dataframes