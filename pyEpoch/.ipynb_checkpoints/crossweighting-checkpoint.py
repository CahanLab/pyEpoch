#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import scanpy as sc
import math
import scipy.signal as ss
from .utils import *

# In[ ]:


def crossweight(adata, lag, minimum, maximum, symmetric_filter=False, weightThresh=0):

    grnDF = adata.uns['grnDF'].copy()

    expDat = makeExpMat(adata)
    #dgenes = adata.uns['dgenes'].copy()
    cells = adata.uns['cells'].index.tolist()

    expDat = expDat[cells]

    # genes=adata.var.index
    # expDat=pd.DataFrame(adata.X).T
    # expDat.columns=adata.obs.index
    # expDat.index=genes
    # expDat=expDat.loc[expDat.sum(axis=1)!=0]
    # expSmoothed=expDat

    
    # lag=math.floor(len(expDat.columns)/5)
    # minimum=math.ceil(len(expDat.columns)/50)
    # maximum=math.floor(len(expDat.columns)/12)
    offset=grnDF.apply(cross_corr,axis=1,expSmoothed=expDat,lag=lag)
    #print(offset)
    grnDF["offset"]=offset
    weighted_scores=[]
    for i in np.arange(grnDF.shape[0]):
        new=score_offset(grnDF["zscore"][i],grnDF["offset"][i],expDat)
        weighted_scores.append(new)
    grnDF["weighted_score"]=weighted_scores
    grnDF = grnDF[grnDF["weighted_score"]>weightThresh]
    adata.uns["grnDF"]=grnDF
    print("Done. Cross-weighted and updated GRN stored in .uns['grnDF'].")
    return adata


# In[ ]:


def cross_corr(grnDF,expSmoothed,lag):
    #grn_row=grnTab
    tg=str(grnDF["TG"])
    tf=str(grnDF["TF"])
    targets=expSmoothed.loc[tg].values
    transcription=expSmoothed.loc[tf].values
    x = ccf(targets,transcription,lag=lag)
    return x


# In[ ]:


def ccf(x, y, lag):
    result = ss.correlate(y - np.mean(y), x - np.mean(x), method='direct') / (np.std(y) * np.std(x) * len(y))
    length = (len(result) - 1) // 2
    lo = length - lag
    hi = length + (lag + 1)
    correlation=result[lo:hi]
    lags=np.arange(lo,hi)-length
    #print(lags)
    df=pd.DataFrame()
    df["lag"]=lags
    df["correlation"]=correlation
    df= df.sort_values('correlation',ascending=False)
    result= np.mean(df["lag"][0:(math.ceil(2/3)*lag)])
    return result


# In[ ]:


def score_offset(score,offset,expSmoothed):
    mini=math.ceil(len(expSmoothed.columns)/50)
    maxi=math.floor(len(expSmoothed.columns)/12)
    if offset<=mini:
        res=score
    elif offset>maxi:
        res=0
    else:
        #linear weighting scheme according to y=(-1/max-min)+1
        weight=-(-offset/(maxi-mini))+1
        res=score*weight
    return res

