#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def crossweight(grnDF, expSmoothed):
    lag=math.floor(len(expSmoothed.columns)/5)
    minimum=math.ceil(len(expSmoothed.columns)/50)
    maximum=math.floor(len(expSmoothed.columns)/12)
    offset=grnDF.apply(cross_corr,axis=1,expSmoothed=expSmoothed,lag=lag)
    #print(offset)
    grnDF["offset"]=offset
    weighted_scores=[]
    for i in np.arange(grnDF.shape[0]):
        new=score_offset(grnDF["zscore"][i],grnDF["offset"][i],expSmoothed)
        weighted_scores.append(new)
    grnDF["weighted_score"]=weighted_scores
    return grnDF


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

