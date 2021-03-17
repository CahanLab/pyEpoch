#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def findDynGenes(expDat, sampTab, group_column="leiden", pseudotime_column="dpt_pseudotime"):


    sampTab["dpt_groups"]=sampTab[group_column]
    sampTab["psuedotime"]=sampTab[pseudotime_column]
    sampTab["cell_name"]=sampTab.index
    path=np.unique(sampTab["leiden"])
    ids=[]
    for grp in path:
        a=sampTab.loc[sampTab["dpt_groups"]==grp]
        b=a["cell_name"]
        ids=np.append(ids,b)
    sampTab=sampTab.loc[ids,:]
    expDat=expDat[ids]
    t1=sampTab["psuedotime"]
    t1C=t1[ids]
    print("starting gamma...")

    gpChr=gamFit(expMat=expDat[t1C.index],genes=expDat.index,celltime=t1C)

    cells=pd.DataFrame()
    cells["cell_name"]=pd.DataFrame(t1).index
    cells["psuedotime"]=t1.values
    cells["group"]=sampTab["dpt_groups"].values
    cells.index=cells["cell_name"]
    cells=cells.sort_values(by="psuedotime")
    ans=list([gpChr,cells])
    return ans


# In[ ]:


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
    if sym == False:
        sys.exit("Error. Please Enter Symmetric Matrix.")
    else:
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


def reconstructGRN(exp,tfs,zThresh,method="pearson"):
    texp=exp.T
    mim=build_mim(exp,method)
    xnet=clr(mim)
    xcorr=texp.corr()
    tfsI= list(set(tfs) & set(texp.columns))
    #print(tfsI)
    xnet=xnet[tfsI]
    xcorr=xcorr[tfsI]
    grn=cn_extractRegsDF(xnet,xcorr,zThresh)
    return grn


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
    grn["Zscores"]=zscoresX
    grn["Correlation"]=correlations
    return grn

