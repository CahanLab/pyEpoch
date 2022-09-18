#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import scanpy as sc
import loess as lo
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.metrics import normalized_mutual_info_score
import sys
from scipy import stats
from pygam import GAM, s,l
from .utils import *



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


# def loessFit(expMat,genes,celltime):

#     genes2=(set(genes) & set(expMat.index))
#     def abcd(input_data):
#         z=pd.DataFrame()
#         z["z"]=input_data.values
#         z["t"]=celltime.values
#         z.index=expMat.columns
#         X=celltime.values.reshape((celltime.shape[0],1))
#         y=z["z"].values

#         xout, yout, wout = loess_1d(X, y, xnew=None, degree=1, frac=0.2, npoints=None, rotate=False, sigy=None)
        
#         # gam=GAM(l(0)).fit(X,y)
#         # p=gam.statistics_['p_values'][0]
#         return p
#     ans=expMat.loc[genes2][celltime.index].apply(abcd,axis=1)
#     return ans


# In[ ]:


def findDynGenes(adata, pseudotime_column="dpt_pseudotime", group_column=None, path=None,pThresh=0.05):
    
    sampTab = makeSampTab(adata)
    expDat = makeExpMat(adata)

    #sampTab=pd.DataFrame(adata.obs)
    #sampTab.rename(columns={'psuedotime':'pseudotime'}, inplace=True)
    
    #genes=adata.var.index
    #expDat=pd.DataFrame(adata.X).T
    #expDat.columns=sampTab.index
    #expDat.index=genes

    sampTab["pseudotime"]=sampTab[pseudotime_column]
    sampTab["cell_name"]=sampTab.index

    if group_column is not None:
        sampTab['group'] = sampTab[group_column]
        if path is not None:
            sampTab = sampTab.loc[sampTab['group'].isin(path)]
    else:
        sampTab = sampTab.assign(group='1')

    # sampTab["dpt_groups"]=sampTab[group_column]
    # sampTab["pseudotime"]=sampTab[pseudotime_column]
    # sampTab["cell_name"]=sampTab.index
    # path=np.unique(sampTab["dpt_groups"])
    # ids=[]
    # for grp in path:
    #     a=sampTab.loc[sampTab["dpt_groups"]==grp]
    #     b=a["cell_name"]
    #     ids=np.append(ids,b)
    # sampTab=sampTab.loc[ids,:]
    # #print(sampTab)

    expDat = expDat[sampTab.index]
    expDat = expDat.loc[expDat.sum(axis=1)!=0]

    #expDat=expDat[ids]
    t1=sampTab["pseudotime"]
    t1C=t1[sampTab.index]
    print("starting gamma...")
    #print(expDat[t1C.index])
    gpChr=pd.DataFrame(gamFit(expDat[t1C.index],expDat.index,t1))
    gpChr.columns=["dynamic_pval"]
    

    cells=pd.DataFrame()
    cells["cell_name"]=pd.DataFrame(t1).index
    cells["pseudotime"]=t1.values
    cells["group"]=sampTab["group"].values
    cells.index=cells["cell_name"]
    cells=cells.sort_values(by="pseudotime")
    #ans=list([gpChr,cells])
    adata.uns["genes"]=gpChr
    adata.uns["cells"]=cells
    print("Done. Dynamic pvalues stored in .uns['genes']. Ordered cells and pseudotime stored in .uns['cells'].")
    
    dgenes = gpChr.index[gpChr['dynamic_pval']<pThresh].tolist()
    adata.uns['dgenes']=dgenes

    print("Dynamic genes stored in .uns['dgenes'].")


    return adata


# In[ ]:

# def build_mim(exp,estimator="pearson"):
#     genes=exp.index
#     if estimator=="pearson":
#         est=KBinsDiscretizer(n_bins=int(np.sqrt(exp.shape[0])),encode='ordinal',strategy="uniform")
#         est.fit(exp.T)
#         dataset=est.transform(exp.T)
#         dataset=pd.DataFrame(dataset)
#         #dataset=dataset+1
#         mim=dataset.corr(method="pearson")**2
#         np.fill_diagonal(mim.values, 0)
#         maxi=0.999999
#         mim[mim>maxi]=maxi
#         mim=-.5*np.log(1-mim)
#         mim[mim<0]=0
#         mim.index=genes
#         mim.columns=genes
#         return mim
#     elif estimator=="MI":
#         print('discretize')
#         est=KBinsDiscretizer(n_bins=int(np.sqrt(exp.shape[0])),encode='ordinal',strategy="uniform")
#         est.fit(exp.T)
#         dataset=est.transform(exp.T)
#         print("MI")
#         nrows = dataset.shape[0]
#         ncols = dataset.shape[1]
#         print(nrows)
#         print(ncols)
#         mim = buildMIM(dataset, nrows, ncols)        
#         return mim
#     else:
#         sys.exit("Must enter valid estimator: MI or pearson.")

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
        print('discretize')
        est=KBinsDiscretizer(n_bins=int(np.sqrt(exp.shape[0])),encode='ordinal',strategy="uniform")
        est.fit(exp.T)
        dataset=est.transform(exp.T)
        #dataset=pd.DataFrame(dataset)
        #dataset=dataset+1
        #mim=pd.DataFrame()#create empty mim dataframe
        #create dataframe with mutual infrormation for each gene pair 
        #texp is the expression matrix with the genes as columns and the samples as rows
        print("MI")
        #mim = np.array([[normalized_mutual_info_score(dataset[:,i],dataset[:,j]) for i in range(dataset.shape[1])] for j in range(dataset.shape[1])])
        mi = [[normalized_mutual_info_score(dataset[:,i],dataset[:,j]) for i in range(0,j+1)] for j in range(dataset.shape[1])]
        mim = np.zeros((dataset.shape[1],dataset.shape[1]))
        mask = np.tril(np.ones((dataset.shape[1],dataset.shape[1])),0) != 0
        mim[mask]=np.concatenate(mi)

        i_upper = np.triu_indices(dataset.shape[1], 0)
        mim[i_upper] = mim.T[i_upper]
        np.fill_diagonal(mim,0)

        print("format")
        mim = pd.DataFrame(mim)
        #mim=mim.set_index(dataset.columns)
        #mim=mim.replace(1, 0)  # for replacing 1 to 0
        mim.index=genes
        mim.columns=genes
        return mim
    else:
        sys.exit("Must enter valid estimator: MI or pearson.")



# def buildMIM(data, nrows, ncols): ## choice is the choice of entropy
#     ##res = pd.DataFrame(index=range(nrows),columns=range(ncols))
#     # res = np.ndarray(shape=(nrows, ncols))
#     res = [None] * (ncols * ncols)
#     sel = [None]* ncols
#     for i in range(ncols):
#         res[i*(ncols) + i] = 0
#         sel[i] = False
    
    
#     for i in range(ncols):
#         sel[i] = True
#         res[i*(ncols) + i] = entropy(data, nrows, ncols, sel)
#         # print(res[i*(ncols) + i]) 
#         sel[i] = False
        
#     for i in range(1, ncols):
#         sel[i] = True
#         for j in range(i):
#             sel[j] = True
#             # print(res)
#             res[j*(ncols) + i] = res[i*(ncols) + i] + res[j*(ncols) + j] - entropy(data, nrows, ncols, sel)
#             res[i*(ncols) + j] = res[j*(ncols) + i]
#             # try:
#             #     res[j*(ncols) + i] = res[i*(ncols) + i] + res[j*(ncols) + j] - entropy(data, nrows, ncols, sel)
#             #     res[i*(ncols) + j] = res[j*(ncols) + i]
#             # except IndexError:
#             #     break
#             sel[j] = False
#         sel[i] = False
#     # print(res) 
#     res = np.nan_to_num(res)
#     res = np.reshape(res, (ncols, ncols))
#     res = pd.DataFrame(res)
#     print(res)
#     return res



# def entropy(data, nrows, ncols, v):
#     data = data.flatten()
#     freq = dict()
#     sel = []
#     ok = True
#     nsample_ok = 0
#     H = 0
#     for s in range(nrows):
#         ok = True
#         sel.clear()
#         for i in range(ncols):
#             if v[i] == True:
#                 if data[s+(i*(nrows))] is not None:
#                     sel.append(data[s+(i*nrows)])
#                 else:
#                     ok = False
#         if ok == True:
#             try:
#                 freq[tuple(sel)] += 1
#             except:
#                 freq[tuple(sel)] = 0
                
#             nsample_ok += 1
    
#     H = emp_entropy(freq,nsample_ok)
#     return H            
  
    
    

    
    
# def emp_entropy(freq: dict, nsample: int):
#     e = 0.0
#     for x in freq.values():
#         e -= x*np.log(x)
#     return np.log(nsample) + e/nsample



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


def reconstructGRN(adata,tfs,zThresh=0,method="pearson"):
    
    expDat = makeExpMat(adata)
    adata.uns['tfs'] = tfs
    print("tfs stored in .uns['tfs']")
    dgenes = adata.uns['dgenes']
    cells = adata.uns['cells'].index.tolist()

    # limit adata to dynamic genes
    # DataFrameGenes=pd.DataFrame(adata.uns["genes"]["dynamic_pval"]<pThresh)
    # dgenes=DataFrameGenes[DataFrameGenes["dynamic_pval"]==True].index.values
    # adata = adata[:,dgenes]

    expDat = expDat[cells]
    expDat = expDat.loc[dgenes]
    expDat = expDat.loc[expDat.sum(axis=1)!=0]

    # reconstruct
    # genes=adata.var.index
    # expDat=pd.DataFrame(adata.X).T
    # expDat.columns=adata.obs.index
    # expDat.index=genes
    # exp=expDat.loc[expDat.sum(axis=1)!=0]
    
    texp=expDat.T
    mim=build_mim(expDat,method)
    xnet=clr(mim)
    xcorr=texp.corr()
    tfsI= list(set(tfs) & set(texp.columns))
    xnet=xnet.loc[tfsI]
    xcorr=xcorr.loc[tfsI]
    grn=cn_extractRegsDF(xnet,xcorr,zThresh)
    grn = grn[grn['zscore']>zThresh]
    
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
    grn["corr"]=correlations
    return grn