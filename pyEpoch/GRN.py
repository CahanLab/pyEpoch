#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.metrics import normalized_mutual_info_score
import sys
from scipy import stats
from pygam import GAM, s,l



# In[ ]:


#from .utils import *


# In[ ]:


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


# In[ ]:


def findDynGenes(adata, group_column="leiden", pseudotime_column="dpt_pseudotime"):
    
    sampTab=pd.DataFrame(adata.obs)
    sampTab.rename(columns={'psuedotime':'pseudotime'}, inplace=True)
    
    genes=adata.var.index
    expDat=pd.DataFrame(adata.X).T
    expDat.columns=sampTab.index
    expDat.index=genes
    expDat=expDat.loc[expDat.sum(axis=1)!=0]



    sampTab["dpt_groups"]=sampTab[group_column]
    sampTab["pseudotime"]=sampTab[pseudotime_column]
    sampTab["cell_name"]=sampTab.index
    path=np.unique(sampTab["leiden"])
    ids=[]
    for grp in path:
        a=sampTab.loc[sampTab["dpt_groups"]==grp]
        b=a["cell_name"]
        ids=np.append(ids,b)
    sampTab=sampTab.loc[ids,:]
    #print(sampTab)
    expDat=expDat[ids]
    t1=sampTab["pseudotime"]
    t1C=t1[ids]
    print("starting gamma...")
    #print(expDat[t1C.index])
    gpChr=pd.DataFrame(gamFit(expDat[t1C.index],expDat.index,t1))
    gpChr.columns=["dynamic_pval"]
    

    cells=pd.DataFrame()
    cells["cell_name"]=pd.DataFrame(t1).index
    cells["pseudotime"]=t1.values
    cells["group"]=sampTab["dpt_groups"].values
    cells.index=cells["cell_name"]
    cells=cells.sort_values(by="pseudotime")
    #ans=list([gpChr,cells])
    adata.uns["genes"]=gpChr
    adata.uns["cells"]=cells
    print("Done. Dynamic pvalues stored in .uns['genes']. Ordered cells and pseudotime stored in .uns['cells'].")
    return adata


# In[ ]:


#def build_mim(exp,estimator="pearson"):
def build_mim(exp,estimator="pearson"):
    genes=exp.index
    
    if estimator=="pearson":
        est=KBinsDiscretizer(n_bins=int(np.sqrt(exp.shape[0])),encode='ordinal',strategy="uniform")
        est.fit(exp.T)
        dataset=est.transform(exp.T)
        dataset=pd.DataFrame(dataset)
        #dataset=dataset+1
        mim=dataset.corr(method="pearson")**2
        np.fill_diagonal(mim.values, 0)
        maxi=0.999999
        mim[mim>maxi]=maxi
        mim=-.5*np.log(1-mim)
        mim[mim<0]=0
        mim.index=genes
        mim.columns=genes
        return mim
    elif estimator=="MI":
        est=KBinsDiscretizer(n_bins=int(np.sqrt(exp.shape[0])),encode='ordinal',strategy="uniform")
        est.fit(exp.T)
        dataset=est.transform(exp.T)
        dataset=pd.DataFrame(dataset)
        #dataset=dataset+1
        mim=pd.DataFrame()#create empty mim dataframe
        #create dataframe with mutual infrormation for each gene pair 
        #texp is the expression matrix with the genes as columns and the samples as rows
        for i in dataset.columns:
            throwaway=[]
            for j in dataset.columns:
                a=dataset[i] 
                b=dataset[j]
                mi_pair=normalized_mutual_info_score(a,b)
                throwaway.append(mi_pair)
            mim[i]=throwaway
        mim=mim.set_index(dataset.columns)
        mim=mim.replace(1, 0)  # for replacing 1 to 0
        mim.index=genes
        mim.columns=genes
        return mim
        #return dataset
        #dataset= dataset.astype('int32') 
        #a=drv.information_mutual_normalised(dataset)
        #dataset.dtype
        #mim=pd.DataFrame(a)
        #return a
    else:
        sys.exit("Must enter valid estimator: MI or pearson.")


# In[ ]:


def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)


# In[ ]:


def clr(mim):
    if len(mim.columns) != len(mim.index):
        sys.exit("Arguement matrix must be square")
    sym=check_symmetric(mim)
    #if sym == False:
    #    sys.exit("Error. Please Enter Symmetric Matrix.")
    #else:
    genes=mim.index
    a=pd.DataFrame(stats.zscore(mim,axis=0))
    a[a<0]=0
    b=pd.DataFrame(stats.zscore(mim,axis=1))
    b[b<0]=0
    xnet=np.sqrt(a**2+b**2)
    xnet.index=genes
    xnet.columns=genes
    return xnet



# In[ ]:


def reconstructGRN(adata,tfs,pThresh=0.05,zThresh=0,method="pearson"):
    # limit adata to dynamic genes
    DataFrameGenes=pd.DataFrame(adata.uns["genes"]["dynamic_pval"]<pThresh)
    dgenes=DataFrameGenes[DataFrameGenes["dynamic_pval"]==True].index.values
    adata = adata[:,dgenes]

    # reconstruct
    genes=adata.var.index
    expDat=pd.DataFrame(adata.X).T
    expDat.columns=adata.obs.index
    expDat.index=genes
    exp=expDat.loc[expDat.sum(axis=1)!=0]
    
    texp=exp.T
    mim=build_mim(exp,method)
    xnet=clr(mim)
    xcorr=texp.corr()
    tfsI= list(set(tfs) & set(texp.columns))
    #print(tfsI)
    xnet=xnet[tfsI]
    xcorr=xcorr[tfsI]
    grn=cn_extractRegsDF(xnet,xcorr,zThresh)
    adata.uns["grnDF"]=grn
    print("Done. Static GRN stored in .uns['grnDF']")
    return adata


# In[ ]:

def cn_extractRegsDF(xnet, xcorr,zThresh):
    targets=[None]*10000000
    regulators=[None]*10000000
    zscoresX=[0]*10000000
    correlations=[0]*10000000
    start=0
    stp=0
    genes=xnet.columns
    for i in genes:
        x=pd.DataFrame(xnet[i])
        regs=list(x.loc[(x[i]>zThresh)].index)#zthresh
        if len(regs)>0:
            zzs=x.loc[regs]
            zzs=zzs[i].to_list()
            corrs=xcorr.loc[regs,i].to_list()
            ncount=len(regs)
            stp=start+ncount
            for j in np.arange(start,stp):
                targets[j]=i
            for j in np.arange(len(regs)):
                regulators[start+j]=regs[j]
            for j in np.arange(len(zzs)):
                zscoresX[start+j]=zzs[j]
            for j in np.arange(len(corrs)):
                correlations[start+j]=corrs[j]
            start=stp
    targets=targets[0:stp]
    regulators=regulators[0:stp]
    zscoresX=zscoresX[0:stp]
    correlations=correlations[0:stp]
    grn=pd.DataFrame()
    grn["TG"]=targets
    grn["TF"]=regulators
    grn["zscore"]=zscoresX
    grn["Correlation"]=correlations
    return grn