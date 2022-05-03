#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import scanpy as sc
import skfda
from pygam import GAM, s,l
from skfda import FDataGrid
import skfda.preprocessing.smoothing.kernel_smoothers as ks


# In[ ]:


def makeExpMat(adata):
    expMat = pd.DataFrame(adata.X.T, index = adata.var_names, columns = adata.obs_names).copy()
    return expMat


# In[ ]:


def makeSampTab(adata):
    sampTab = adata.obs.copy()
    return sampTab


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

