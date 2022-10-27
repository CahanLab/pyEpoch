#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import scanpy as sc
import math
import sys
from scipy import stats
import statsmodels.stats.multitest as multi
import igraph as igraph
from .utils import *



def define_epochs(adata,method='pseudotime',num_epochs=2,pseudotime_cuts=None,group_assignments=None):
    
    #genes=adata.var.index
    #expDat=pd.DataFrame(adata.X).T
    #expDat.columns=adata.obs.index
    #expDat.index=genes
    #expDat=expDat.loc[expDat.sum(axis=1)!=0]
    
    #dyn_expDat=expDat.loc[dgenes,:]

    cells = adata.uns['cells'].copy()
    
    if method=="cell_order":
        t1=cells["pseudotime"]
        t1.index=list(cells["cell_name"].values)
        t1=t1.sort_values(ascending=True)
        chunk_size=math.floor(len(t1)/num_epochs)

        epochs=[]
        count=1
        for i in np.arange(num_epochs):
            epochs[(((i)*chunk_size)):((i+1)*chunk_size)]=[count]*chunk_size
            count=count+1

        if len(epochs)!= len(t1):
            for i in np.arange(len(t1)-len(epochs)):
                epochs.append(epochs[len(epochs)-1])

        numbers=np.arange(num_epochs)+1
        epoch_names=[]
        for i in epochs:
            for j in numbers:
                if i==j:
                    epoch_names.append("epoch"+str(j))

        cells["epoch"]=epoch_names
        adata.uns['cells'] = cells
        print("Done. Updated ordered cells with epoch assignments stored in .uns['cells'].")
        return adata

    if method=="pseudotime":
        if (pseudotime_cuts is not None):
            a = split_epochs_by_pseudotime(adata,pseudotime_cuts)
        else:
            pseudotime_cuts=np.linspace(0, np.max(cells["pseudotime"]),num_epochs+1)
            pseudotime_cuts=pseudotime_cuts[1:]
            pseudotime_cuts=pseudotime_cuts[:(len(pseudotime_cuts)-1)]
            pseudotime_cuts=pseudotime_cuts.tolist()
            a=split_epochs_by_pseudotime(adata,pseudotime_cuts)
        print("Done. Updated ordered cells with epoch assignments stored in .uns['cells'].")
        return a

    if method=="group":
        if group_assignments is None:
            sys.exit("Must provide group_assignments for group method.")
        a=split_epochs_by_group(adata,group_assignments)
        print("Done. Updated ordered cells with epoch assignments stored in .uns['cells'].")
        return a



# In[2]:

def split_epochs_by_pseudotime(adata,pseudotime_cuts,epoch_names2=None):
    sampTab=adata.uns['cells'].copy()
    if max(pseudotime_cuts)>max(sampTab["pseudotime"]):
        sys.exit("Cuts must be within pseudotime.")
    if (epoch_names2!=None):
        if (len(epoch_names2)!=len(pseudotime_cuts)+1):
            sys.exit("Length of epoch_names must be equal to 1+length(cuts).")
    if epoch_names2==None:
        epoch_names2=[]
        for i in np.arange(len(pseudotime_cuts)+1):
            epoch_names2.append("epoch"+str(i))
    pseudotime_cuts.insert(0,-.1)
    pseudotime_cuts.append(np.max(sampTab["pseudotime"])+.1)
    epoch=[]
    for i in np.arange(len(pseudotime_cuts)-1):
        b=sampTab.loc[(sampTab['pseudotime']>=pseudotime_cuts[i]) & (sampTab['pseudotime']<pseudotime_cuts[i+1])]
        for j in np.arange(b.shape[0]):
            epoch.append("epoch"+str(i+1))
    sampTab["epoch"]=epoch
    adata.uns['cells'] = sampTab
    return adata

# In[3]:


#split_epochs_by_group<-function(dynRes,assignment)
def split_epochs_by_group(adata,group_assignments):
    adata.uns['cells']["epoch"]=group_assignments
    return adata


# In[4]:


#assign_epochs<-function(expSmoothed,dynRes,method='active_expression',pThresh_dyn=0.05,pThresh_DE=0.05,toScale=FALSE)
#epoch_assignments<-assign_epochs(expDat[dgenes,],xdyn)
def assign_epochs(adata, method="active_expression",pThresh_dyn=.05,pThresh_DE=.05,toScale=False):
    #method="active_expresion"
    #pThresh_dyn=.05
    #pThresh_DE=.05
    #toScale=False

    exp = makeExpMat(adata)
    dgenes = adata.uns['dgenes'].copy()
    exp = exp.loc[dgenes]

    cells = adata.uns['cells'].copy()

    # genes=adata.var.index
    # expDat=pd.DataFrame(adata.X).T
    # expDat.columns=adata.obs.index
    # expDat.index=genes
    # expDat=expDat.loc[expDat.sum(axis=1)!=0]
    # exp=expDat.loc[adata.uns['genes'].loc[adata.uns['genes']["dynamic_pval"]<pThresh_dyn].index,]

    if toScale==True:
        save_genes=exp.index
        save_cells=exp.columns
        exp=pd.DataFrame(preprocessing.scale(exp))
        exp.index=save_genes
        exp.columns=save_cells

    epoch_names=list(cells["epoch"].unique())
    epochs={}
    for i in epoch_names:
        epochs[i]=None

    navg=math.ceil(len(exp.columns)*.05)
    # compute thresholds for each gene
    # for each gene, order expression --- compute top as average of top 5%, bottom as average of bottom 5% 
    # set threshold (active/inactive) as midpoint (or maybe 1/3 in case gene is expressed at different levels) between top and bottom
    thresholds=pd.DataFrame()
    thresholds["genes"]=exp.index
    #thresholds["thresh"]=[0]*len(exp.index)
    thresholds.index=thresholds["genes"]
    #thresholds

    list_of_thresholds=[]
    for i in list(thresholds.index):
        #i=list(thresholds.index)[0]
        profile=exp.loc[i,].sort_values(ascending=True)
        bottom=np.mean(profile[0:navg])
        top=np.mean(profile[-(navg+1):])#mistake in the R code
        thresh=((top-bottom)*.33)+bottom
        list_of_thresholds.append(thresh)
    thresholds["thresh"]=list_of_thresholds
    mean_expression=pd.DataFrame(columns=["gene","epoch","mean_expression"])

    if method=="active_expression":
        for i in list(epochs.keys()):
            chunk_cells=list(cells.loc[cells["epoch"]==i]["cell_name"])
            chunk=exp[chunk_cells]
            chunk_df=pd.DataFrame(np.mean(chunk,axis=1),columns=["means"],index=np.mean(chunk,axis=1).index)
            chunk_df=pd.concat([chunk_df,thresholds],axis=1)
            chunk_df["active"]=chunk_df["means"]>chunk_df["thresh"]
            epochs[i]=list(chunk_df.loc[chunk_df["active"]==True]["genes"])
            a=pd.DataFrame(chunk.index,columns=["gene"],index=chunk.index)
            b=pd.DataFrame([i]*len(chunk.index),columns=["epoch"],index=chunk.index)
            c=pd.DataFrame(np.mean(chunk,axis=1),columns=["mean_expression"],index=chunk.index)
            concatted=pd.concat([a,b,c],axis=1)
            mean_expression=pd.concat([mean_expression,concatted],axis=0)

    else:
        for i in list(epochs.keys()):
            #i=list(epochs.keys())[0]
            chunk_cells=list(cells.loc[cells["epoch"]==i]["cell_name"])
            chunk=exp[chunk_cells]
            column_names=list(np.setdiff1d(list(exp.columns),chunk_cells))
            background=exp[column_names]
            diffres=pd.DataFrame(columns=["gene","mean_diff","pval"])

            for j in list(exp.index):
                #j=list(exp.index)[0]    
                t=stats.ttest_ind(chunk.loc[j,:],background.loc[j,:],equal_var=False)
                single_data=(pd.DataFrame([j,np.mean(chunk.loc[j,:])-np.mean(background.loc[j,:]),t.pvalue],index=["gene","mean_diff","pval"])).T
                diffres=pd.concat([diffres,single_data],axis=0)
                #failed methods
                #res=stat()
                #chunk1=pd.DataFrame(chunk.loc[j,:])
                #chunk1["label"]=["chunk"]*len(chunk1.index)
                #background1=pd.DataFrame(background.loc[j,:])
                #background1["label"]=["background"]*len(background1.index)
                #t_test_data=pd.concat([chunk1,background1],axis=0)
                #res.ttest(df=t_test_data,xfac="label",res=j,test_type=2,evar=False)
                #summary, results = rp.ttest(group1=chunk.loc[j,:],group1_name="chunk",group2=background.loc[j,:],group2_name="background")

            p_adj=multi.multipletests(diffres["pval"],method="fdr_bh")[1]
            diffres["padj"]=p_adj
            diffres.index=diffres["gene"]
            diffres=diffres.loc[(diffres["mean_diff"]>0)] #filter out genes that are on
            epochs[i]=list(diffres.loc[diffres["padj"]<pThresh_DE]["gene"])

            chunk_df=pd.DataFrame(np.mean(chunk,axis=1),columns=["means"],index=np.mean(chunk,axis=1).index)
            chunk_df=pd.concat([chunk_df,thresholds],axis=1)
            chunk_df["active"]=chunk_df["means"]>=chunk_df["thresh"]
            epochs[i]=np.intersect1d(epochs[i],list(chunk_df.loc[chunk_df["active"]==True]["genes"]))

            a=pd.DataFrame(chunk.index,columns=["gene"],index=chunk.index)
            b=pd.DataFrame([i]*len(chunk.index),columns=["epoch"],index=chunk.index)
            c=pd.DataFrame(np.mean(chunk,axis=1),columns=["mean_expression"],index=chunk.index)
            concatted=pd.concat([a,b,c],axis=1)
            mean_expression=pd.concat([mean_expression,concatted],axis=0)


    mean_expression['gene'] = mean_expression.index
    #epochs["mean_expression"]=mean_expression
    mean_expression = mean_expression.reset_index()
    mean_expression["mean_expression"] = pd.to_numeric(mean_expression.mean_expression)
    adata.uns['mean_expression']=mean_expression
    adata.uns['epochs']=epochs

    print("Epoch gene assignments stored in .uns['epochs'].")
    print("Mean expression per epoch stored in .uns['mean_expression'].")

    return adata


# In[5]:


def assign_epochs_simple(adata,num_epochs=2,pThresh_dyn=.01,toScale=False):
    #num_epochs=2
    #pThresh=.01
    #toScale=False


    exp = makeExpMat(adata)
    dgenes = adata.uns['dgenes'].copy()
    exp = exp.loc[dgenes]

    cells = adata.uns['cells'].copy()

    
    # genes=adata.var.index
    # expDat=pd.DataFrame(adata.X).T
    # expDat.columns=adata.obs.index
    # expDat.index=genes
    # expDat=expDat.loc[expDat.sum(axis=1)!=0]
    # #limit to dynamically expressed genes
    # exp=expDat.loc[adata.uns['genes'].loc[adata.uns['genes']["expression"]<pThresh_dyn].index,]

    if toScale==True:
        save_genes=exp.index
        save_cells=exp.columns
        exp=pd.DataFrame(preprocessing.scale(exp))
        exp.index=save_genes
        exp.columns=save_cells

    navg=math.ceil(len(exp.columns)*.05)

    # compute thresholds for each gene
    #for each gene, order expression --- compute top as average of top 5%, bottom as average of bottom 5% 
    # set threshold (active/inactive) as midpoint (or maybe 1/3 in case gene is expressed at different levels) between top and bottom

    thresholds=pd.DataFrame()
    thresholds["genes"]=exp.index
    #thresholds["thresh"]=[0]*len(exp.index)
    thresholds.index=thresholds["genes"]

    list_of_thresholds=[]
    for i in list(thresholds.index):
        #i=list(thresholds.index)[0]
        profile=exp.loc[i,].sort_values(ascending=True)
        bottom=np.mean(profile[0:navg])
        top=np.mean(profile[-(navg+1):])#mistake in the R code
        thresh=((top-bottom)*.33)+bottom
        list_of_thresholds.append(thresh)
    thresholds["thresh"]=list_of_thresholds

    # order cells in exp along pseudotime-- cells ordered in dynRes
    t1=pd.DataFrame()
    t1["pseudotime"]=cells["pseudotime"]
    t1.index=cells["cell_name"]

    t1=t1.sort_values(by=['pseudotime'])
    exp=exp[t1.index]

    mean_expression=pd.DataFrame(columns=["gene","epoch","mean_expression","peakTime"])

    #divide epochs by psuedotime
    epoch_names=list(cells["epochs"].unique())
    epochs={}
    for i in epoch_names:
        epochs[i]=None

    # determine activity based on average expression in each epoch
    ptmax=np.max(cells["pseudotime"])
    ptmin=np.min(cells["pseudotime"])
    chunk_size=(ptmax-ptmin)/num_epochs

    #cellsEps=pd.DataFrame()
    #cells=[]
    #for i in t1.index:
    #    cells.append(i)
    #cellsEps["cells"]=cells

    cellsEps={}
    for i in t1.index:
        cellsEps[i]=None

    for i in (np.arange(len(epochs))+1):

        lower_bound=ptmin+((i-1)*chunk_size)
        upper_bound=ptmin+((i)*chunk_size)
        #chunk_cells
        #chunk_cells<-rownames(dynRes$cells[dynRes$cells$pseudotime>=lower_bound & dynRes$cells$pseudotime<=upper_bound,])
        chunk_cells=list(cells.loc[(cells['pseudotime'] >= lower_bound) & (cells['pseudotime'] <= upper_bound)]["cell_name"])
        chunk=exp[chunk_cells]
        chunk_df=pd.DataFrame()
        chunk_df["means"]=np.mean(chunk,axis=1)
        chunk_df=pd.concat([chunk_df,thresholds],axis=1)
        chunk_df["active"]=chunk_df["means"]>=chunk_df["thresh"]
        epochs[list(epochs.keys())[i-1]]=list(chunk_df.loc[chunk_df["active"]==True]["genes"])
        genesPeakTimes=chunk.apply(np.argmax,axis=1)

        gpt=cells[cells.cell_name.isin(chunk_cells)].iloc[genesPeakTimes.values]["pseudotime"]



        a=pd.DataFrame(chunk.index,columns=["gene"],index=chunk.index)
        b=pd.DataFrame([i]*len(chunk.index),columns=["epoch"],index=chunk.index)
        c=pd.DataFrame(np.mean(chunk,axis=1),columns=["mean_expression"],index=chunk.index)
        d=pd.DataFrame(gpt.values,columns=["peakTime"],index=chunk.index)
        concatted=pd.concat([a,b,c,d],axis=1)
        mean_expression=pd.concat([mean_expression,concatted],axis=0)
        mean_expression.index=np.arange(mean_expression.shape[0])


        for j in chunk_cells:
            cellsEps[j]=epoch_names[i-1]


    #assign genes to epochs
    genes=np.unique(mean_expression["gene"])
    print("n genes:",len(genes))


    eps={}
    geneEpPT={}
    epMean={}

    for i in genes:
        eps[i]=None
        geneEpPT[i]=None
        epMean[i]=None

    for i in genes:

        x=mean_expression.loc[mean_expression["gene"]==i]
        x.index=np.arange(x.shape[0])
        xi=np.argmax(x["mean_expression"])
        eps[i]=x.loc[xi]["epoch"]
        geneEpPT[i]=x.loc[xi]["peakTime"]
        epMean[i]=np.max(x["mean_expression"])

    geneDF=pd.DataFrame()
    geneDF["genes"]=genes
    geneDF["epoch"]=eps.values()
    geneDF["peakTime"]=geneEpPT.values()
    geneDF["epMean"]=epMean.values()
    geneDF["pval"]=adata.uns['genes'].loc[genes]["expression"].values

    cells2=cells.loc[cells["cell_name"].isin(t1.index)]
    cells2["epoch"]=cellsEps.values()

    #epochs["mean_expression"]=mean_expression
    
    adata.uns['epochs']=epochs
    mean_expression = mean_expression.reset_index()
    mean_expression["mean_expression"] = pd.to_numeric(mean_expression.mean_expression)
    adata.uns['mean_expression']=mean_expression

    print("Epoch gene assignments stored in .uns['epochs'].")
    print("Mean expression per epoch stored in .uns['mean_expression'].")
    return adata

# In[ ]:


def epochGRN(adata, epoch_network=None):
    
    grnDF = adata.uns['grnDF'].copy()
    epochs = adata.uns['epochs'].copy()

    #epochs.pop('mean_expression', None)
    keys=epochs.keys()

    all_dyngenes=[]
    for i in epochs.keys():
        if i != "mean_expression":
            for j in epochs[i]:
                all_dyngenes.append(j)

    all_dyngenes=np.unique(all_dyngenes)
    len(all_dyngenes)

    if epoch_network==None:
        epoch_network=pd.DataFrame()
        epoch_network["from"]=None
        epoch_network["to"]=None
        for i in np.arange(len(epochs)-1):
            df=pd.DataFrame()
            df["from"]=[list(epochs.keys())[i]]
            df["to"]=[list(epochs.keys())[i+1]]
            epoch_network=pd.concat([epoch_network,df],axis=0)


    throwaway=pd.DataFrame()
    throwaway["from"]=epochs.keys()
    throwaway["to"]=epochs.keys()
    epoch_network=pd.concat([epoch_network,throwaway])
    epoch_network

    epoch_network=epoch_network.reset_index()
    epoch_network=epoch_network.drop(["index"],axis=1)

    name=[]
    for i in np.arange(epoch_network.shape[0]):
        name.append("..".join(epoch_network.iloc[i].values))
    epoch_network["name"]=name
    print(epoch_network)

    GRN={}
    for i in epoch_network["name"]:
        GRN[i]=None

    for i in np.arange(epoch_network.shape[0]):
        From=epoch_network.iloc[i][0]
        To=epoch_network.iloc[i][1]
        a=pd.Series(grnDF["TF"]).isin(epochs[From])
        temp=grnDF[a]
            # For transition network
        if From != To:
            # target turns on: target is not active in source epoch but is active in target epoch
            # target turns off: target is active in source epoch but is not active in target epoch
            # remove any other interaction (i.e. interactions that are constant -- target on in both epochs or off in both epochs)

            remove_tgs_in_both=list(set(epochs[To]) & set(epochs[From]))
            remove_tgs_in_neither=list(set(all_dyngenes).difference(set(epochs[From])) & set(all_dyngenes).difference(set(epochs[To])))

            b=pd.Series(temp["TG"]).isin(remove_tgs_in_both)
            temp=temp[~b]
            c=pd.Series(temp["TG"]).isin(remove_tgs_in_neither)
            temp=temp[~c]

        GRN[epoch_network["name"][i]]=temp
    adata.uns["dynamic_GRN"]=GRN
    print("Done. Dynamic GRN stored in .uns['dynamic_GRN'].")
    return adata


# In[ ]:


def compute_pagerank(adata,weight_column="weighted_score",directed_graph=False):
    dynnet = adata.uns['dynamic_GRN'].copy()

    keys=list(dynnet.keys())
    ranks={}
    for i in keys:
        ranks[i]=None

    #net=list(dynnet.keys())[0]
    for i in list(dynnet.keys()):
        df=dynnet[i]
        df=df[["TF","TG",weight_column]]
        df.columns=["TF","TG","weight"]
        df.index=np.arange(df.shape[0])


        pagerank_genes=[]
        for j in np.arange(df.shape[0]):
            a=df.iloc[j]
            pagerank_genes.append(a[0])
            pagerank_genes.append(a[1])
        indexes=np.unique(pagerank_genes, return_index=True)[1]
        pagerank_genes=[pagerank_genes[index]for index in sorted(indexes)]  

        tuples = [tuple(x) for x in df.values]
        Gm = igraph.Graph.TupleList(tuples, directed = False, edge_attrs = ['weight'])

        #page.rank(graph, nodes = V(graph), directed = is.directed(graph),
        #niter = 1000, eps = 0.001, damping = 0.85)
        #pagerank<-data.frame(page_rank(tfnet,directed=directed_graph)$vector)
        pagerank=pd.DataFrame()
        pagerank["page_rank"]=Gm.pagerank(directed=False,implementation="arpack",weights=list(df["weight"]))
        pagerank.index=pagerank_genes
        pagerank["gene"]=pagerank.index
        pagerank=pagerank[["gene","page_rank"]]
        pagerank=pagerank.sort_values(by='page_rank', ascending=False)
        #pagerank["is_regulator"]=[False]*pagerank.shape[0]
        #pagerank
        pagerank["is_regulator"]=list(pd.Series(pagerank["gene"]).isin(np.unique(df["TF"])))
        ranks[i]=pagerank
    adata.uns["pagerank"]=ranks

    print("PageRank stored in .uns['pagerank'].")
    return adata


def compute_betweenness_degree(adata,weight_column="weighted_score",directed_graph=False):

    dynnet = adata.uns['dynamic_GRN'].copy()
    
    keys=list(dynnet.keys())
    ranks={}
    for i in keys:
        ranks[i]=None
    for i in list(dynnet.keys()):
        df=dynnet[i]
        df=df[["TF","TG",weight_column]]
        df.columns=["TF","TG","weight"]
        df.index=np.arange(df.shape[0])
        pagerank_genes=[]
        
        for j in np.arange(df.shape[0]):
            a=df.iloc[j]
            pagerank_genes.append(a[0])
            pagerank_genes.append(a[1])
        indexes=np.unique(pagerank_genes, return_index=True)[1]
        pagerank_genes=[pagerank_genes[index]for index in sorted(indexes)]  
        pagerank_genes
        tuples = [tuple(x) for x in df.values]
        Gm = igraph.Graph.TupleList(tuples, directed = False, edge_attrs = ['weight'])
        betweenness_degree=pd.DataFrame()

        betweenness_degree["betweenness"]=Gm.betweenness(directed=False,weights=list(df["weight"]))
        n=Gm.vcount()
        #Bnorm=2*B/(n*n-3*n+2)
        betweenness_degree["betweenness"]=2*betweenness_degree/(n*n-3*n+2)
    
        betweenness_degree["degree"]=Gm.degree(mode='all',loops='true')
        betweenness_degree["degree"]=betweenness_degree["degree"]/(n-1)
    
        betweenness_degree["betweenness*degree"]=betweenness_degree["betweenness"]*betweenness_degree["degree"]
        betweenness_degree.index=pagerank_genes
        betweenness_degree["gene"]=betweenness_degree.index
        betweenness_degree=betweenness_degree[["gene","betweenness","degree","betweenness*degree"]]
    
        betweenness_degree=betweenness_degree.sort_values(by='betweenness*degree', ascending=False)
    
        betweenness_degree["is_regulator"]=list(pd.Series(betweenness_degree["gene"]).isin(np.unique(df["TF"])))
        ranks[i]=betweenness_degree
        
    adata.uns["betweenness_degree"]=ranks
    print("Betweenness-degree stored in .uns['betweenness_degree'].")    
    return adata

