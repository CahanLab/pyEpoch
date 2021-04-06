#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import skfda
from pygam import GAM, s,l
from skfda import FDataGrid
import skfda.preprocessing.smoothing.kernel_smoothers as ks


# In[ ]:


#gamFit<-functin(expMat,genes, # genes to testcelltime){
def gamFit(expMat,genes,celltime):
    #expMat=expDat[t1C.index]
    #genes=expDat.index
    #celltime=t1C

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


def grnKsmooth(expDat,cells,BW=.25):
    #cells=ccells
    #BW=.1
    BW=min(BW, max(cells["pseudotime"])-min(cells["pseudotime"])/10)
    t1=pd.DataFrame(cells["pseudotime"])
    t1.index=cells["cell_name"]
    t1=t1.sort_values(by='pseudotime', ascending=True)
    #expDat.iloc[t1.index]
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
    return ans

