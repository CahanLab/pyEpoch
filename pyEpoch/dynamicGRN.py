#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import scanpy as sc
import math
import sys
from scipy import stats
import statsmodels.stats.multitest as multi
import igraph as igraph


# In[1]:


def define_epochs(xdyn,dyn_expDat,method,num_epochs,group_assignments=None):
    if method=="cell_order":
        t1=xdyn[1]["pseudotime"]
        t1.index=list(xdyn[1]["cell_name"].values)
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
                    epoch_names.append("Epoch"+str(j))

        xdyn[1]["epochs"]=epoch_names
        return xdyn

    if method=="pseudotime":
        pseudotime_cuts=np.linspace(0, np.max(xdyn[1]["pseudotime"]),num_epochs+1)
        pseudotime_cuts=pseudotime_cuts[1:]
        pseudotime_cuts=pseudotime_cuts[:(len(pseudotime_cuts)-1)]
        pseudotime_cuts=pseudotime_cuts.tolist()
        a=split_epochs_by_pseudotime(xdyn,pseudotime_cuts)
        return a

    if method=="group":
        if group_assignments==None:
            sys.exit("Must provide group_assignments for group method.")
        a=split_epochs_by_group(xdyn,group_assignments)
        return a


# In[2]:


def split_epochs_by_pseudotime(xdyn,pseudotime_cuts,epoch_names2=None):
    sampTab=xdyn[1]
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
    xdyn[1]["epochs"]=epoch
    return xdyn


# In[3]:


#split_epochs_by_group<-function(dynRes,assignment)
def split_epochs_by_group(xdyn,group_assignments):
    xdyn[0]["epochs"]=group_assignments
    return xdyn


# In[4]:


def assign_epochs(expSmoothed, xdyn, method="active_expression",pThresh_dyn=.05,pThresh_DE=.05,toScale=False):
    #method="active_expresion"
    #pThresh_dyn=.05
    #pThresh_DE=.05
    #toScale=False
    exp=expSmoothed.loc[xdyn[0].loc[xdyn[0]["expression"]<pThresh_dyn].index,]

    if toScale==True:
        save_genes=exp.index
        save_cells=exp.columns
        exp=pd.DataFrame(preprocessing.scale(exp))
        exp.index=save_genes
        exp.columns=save_cells

    epoch_names=list(xdyn[1]["epochs"].unique())
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
            chunk_cells=list(xdyn[1].loc[xdyn[1]["epochs"]==i]["cell_name"])
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
            chunk_cells=list(xdyn[1].loc[xdyn[1]["epochs"]==i]["cell_name"])
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


    epochs["mean_expression"]=mean_expression
    return epochs


# In[5]:


def assign_epochs_simple(expSmoothed,xdyn,num_epochs=2,pThresh=.01,toScale=False):
    #num_epochs=2
    #pThresh=.01
    #toScale=False

    #limit to dynamically expressed genes
    exp=expSmoothed.loc[xdyn[0].loc[xdyn[0]["expression"]<pThresh].index,]

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
    t1["pseudotime"]=xdyn[1]["pseudotime"]
    t1.index=xdyn[1]["cell_name"]

    t1=t1.sort_values(by=['pseudotime'])
    exp=exp[t1.index]

    mean_expression=pd.DataFrame(columns=["gene","epoch","mean_expression","peakTime"])

    #divide epochs by psuedotime
    epoch_names=list(xdyn[1]["epochs"].unique())
    epochs={}
    for i in epoch_names:
        epochs[i]=None

    # determine activity based on average expression in each epoch
    ptmax=np.max(xdyn[1]["pseudotime"])
    ptmin=np.min(xdyn[1]["pseudotime"])
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
        chunk_cells=list(xdyn[1].loc[(xdyn[1]['pseudotime'] >= lower_bound) & (xdyn[1]['pseudotime'] <= upper_bound)]["cell_name"])
        chunk=exp[chunk_cells]
        chunk_df=pd.DataFrame()
        chunk_df["means"]=np.mean(chunk,axis=1)
        chunk_df=pd.concat([chunk_df,thresholds],axis=1)
        chunk_df["active"]=chunk_df["means"]>=chunk_df["thresh"]
        epochs[list(epochs.keys())[i-1]]=list(chunk_df.loc[chunk_df["active"]==True]["genes"])
        genesPeakTimes=chunk.apply(np.argmax,axis=1)

        gpt=xdyn[1][xdyn[1].cell_name.isin(chunk_cells)].iloc[genesPeakTimes.values]["pseudotime"]



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
    geneDF["pval"]=xdyn[0].loc[genes]["expression"].values

    cells2=xdyn[1].loc[xdyn[1]["cell_name"].isin(t1.index)]
    cells2["epoch"]=cellsEps.values()

    epochs["mean_expression"]=mean_expression
    return epochs


# In[ ]:


def epochGRN(grnDF, epochs, epoch_network=None):
    #epoch_network=None
    #epochs=epoch_assignments
    
    #keys=epochs.keys()
    #for i in keys:
    #    if i=="mean_expression":
    #        del epochs["mean_expression"]
    epochs.pop('mean_expression', None)
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
    return GRN


# In[ ]:


def compute_pagerank(dynnet,weight_column="weighted_score",directed_graph=False):
    #dynnet=dynamic_grn
    #weight_column="weighted_score"
    #directed_graph=False

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
    return ranks

