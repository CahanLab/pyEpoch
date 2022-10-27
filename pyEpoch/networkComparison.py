import numpy as np
import pandas as pd
import sys
import igraph as igraph
from collections import defaultdict
from .utils import *



def edge_uniqueness(adata_list,grn='grnDF',tfs=None,weight_column="weighted_score"):
    grnDFs=[]
    for a in adata_list:
        grnDFs.append(a.uns[grn])

    if tfs is None:
        tfs = adata_list[0].uns['tfs']

    # grn_DFs = 
    genes=[]
    for g in grnDFs:
        values=list(g["TG"])+list(g["TF"])
        genes=genes+values
    genes=np.unique(genes)
    tfs=list(set(genes) & set(tfs))

    adj_list=[]

    for g in grnDFs:
        df = g[['TF','TG',weight_column]]

        ig = igraph.Graph.TupleList(df.itertuples(index=False), directed=True, weights=True)
        vtcs=igraph.VertexSeq(ig)
        to_add = list(set(genes)-set(vtcs['name']))

        ig.add_vertices(to_add)

        adjmat = pd.DataFrame(ig.get_adjacency(attribute="weight",default=0),index=ig.vs['name'],columns=ig.vs['name'])
        adjmat = adjmat.loc[tfs]
        # adjmat = adjmat.T

        adj_list.append(adjmat)

    res=edge_rank(adj_list)

    return res



def edge_rank(adj_list):
    #grnMats=adj_list

    nums = [*range(1,len(adj_list)+1)]
    names = ["net" + str(x) for x in nums]

    adj_dfs = []
    for i in range(len(names)):
        adj = adj_list[i]
        adj['TF']=adj.index
        df = adj.melt(id_vars=['TF'])
        df.columns=['TF','TG',names[i]]
        adj_dfs.append(df)

    res = adj_dfs[0][['TF','TG']]
    for df in adj_dfs:
        res = res.merge(df,on=['TF','TG'],how="left")

    # get max and mins
    res.fillna(0)
    res['max'] = res.max(axis=1,numeric_only=True)
    res['min'] = res.min(axis=1,numeric_only=True)
    res['diff'] = res['max'] - res['min']

    res = res.sort_values(by='diff',ascending=False)

    return res




# Condition is the column in edgeDF corresponding to the network of interest
def dynamic_difference_network(edgeDF, adata_list, condition, type="on", diff_thresh=3, condition_thresh=6):
    #edgeDF=res
    #epochs=[adata1,adata2]
    #condition="network1"
    #types="on"
    #diff_thresh=7.5
    #condition_thresh=10
    
    edgeDF = edgeDF.loc[edgeDF['diff']!=0]
    conditions = list(set(list(edgeDF.columns))-set(['TG', 'TF','min', 'max', 'diff']))
    if condition not in conditions:
        sys.exit("condition not represented.")

    epochs_combined = adata_list[0].uns['epochs']
    for ad in adata_list.pop(0):
        epochs_combined={key:np.hstack([epochs_combined[key],ad.uns['epochs'][key]]) for key in epochs_combined.keys()}

    #https://stackoverflow.com/questions/54040858/python-3-x-merge-two-dictionaries-with-same-keys-and-values-being-array
    #d3 = {key:np.hstack([d1[key],d2[key]]) for key in d1.keys()}

    # epochs=[]
    # for ad in adata_list:
    #     epochs.append(ad.uns['epochs'])

    # epochs_combined = defaultdict(list)
    # for d in epochs:
    #     for key,value in d.items():
    #         epochs_combined[key].append(value)
    for d in epochs_combined.keys():
        epochs_combined[d] = list(set(epochs_combined[d]))

    dynamic_edges = dict.fromkeys(epochs_combined.keys())
    for d in epochs_combined.keys():
        dynamic_edges[d] = edgeDF.loc[edgeDF['TF'].isin(epochs_combined[d])]

    diffnet = dict.fromkeys(dynamic_edges.keys())
    for d in diffnet.keys():
        diffnet[d] = dynamic_edges[d].drop(['min','max'],axis=1)
        diffnet[d] = diffnet[d].loc[diffnet[d]['diff']>diff_thresh]

    if type=="on":
        # Get edges specifically on in condition (is the max)
        for d in diffnet.keys():
            diffnet[d]['win'] = diffnet[d].drop(['diff'],axis=1).select_dtypes(np.number).idxmax(axis=1)
            diffnet[d] = diffnet[d].loc[diffnet[d]['win']==condition]

            diffnet[d] = diffnet[d].drop(['win'],axis=1)
            diffnet[d] = diffnet[d].loc[diffnet[d][condition]>=condition_thresh]
    elif type=="off":
        # Get edges specifically not in condition (is the min)
        for d in diffnet.keys():
            diffnet[d]['win'] = diffnet[d].drop(['diff'],axis=1).select_dtypes(np.number).idxmin(axis=1)
            diffnet[d] = diffnet[d].loc[diffnet[d]['win']==condition]

            diffnet[d] = diffnet[d].drop(['win'],axis=1)
            diffnet[d] = diffnet[d].loc[diffnet[d][condition]<condition_thresh]
    else:
        sys.exit("type should be 'on' or 'off'.")

    return diffnet


# diffnet is the result of runnin dynamic_difference_network
# adata contains the expression data that will be used to assign the interaction type
#       If looking for edges specific to network 1, use adata corresponding to network 1...
def add_diffnet_interaction_type(diffnet,adata):
    if type(diffnet)==pd.core.frame.DataFrame:
        GRN=diffnet
        GRN = add_corr(GRN,adata)

        true_false=GRN["corr"]>0
        activation=[]
        for i in true_false:
            if i==True:
                activation.append("activation")
            else:
                activation.append("repression")
        GRN["interaction"]=activation

    elif type(diffnet)==dict:
        GRN=diffnet
        for i in GRN:
            GRN[i] = add_corr(GRN[i],adata)

            true_false=GRN[i]["corr"]>0
            activation=[]
            for j in true_false:
                if j==True:
                    activation.append("activation")
                else:
                    activation.append("repression")
            GRN[i]["interaction"]=activation

    return GRN





