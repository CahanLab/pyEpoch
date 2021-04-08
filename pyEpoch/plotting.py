#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import matplotlib.pyplot as plt
import sys
import seaborn as sns
from matplotlib.lines import Line2D


# In[ ]:


def plot_dynamic_network(grn,tfs,only_TFs=True,order=None,thresh=None):
    if type(grn)==dict:
        if order != None:
            grn = {k: grn[k] for k in order}


        for i in list(grn.keys()):
            print(i)
            df=grn[i]

            if only_TFs==True:
                df=df.loc[df["TG"].isin(tfs)]

            if thresh != None:
                df=df.loc[df["zscore"]>thresh]

            interactions=[]
            for k in df["Correlation"]:
                if k>0:
                    interactions.append("activation")
                else:
                    interactions.append("repression")
            df["interactions"]=interactions


            G=nx.Graph()
            for j in df.index:
                a=df.loc[j]
                if a["interactions"]=="activation":
                    G.add_edge(a["TF"],a["TG"],color="b")
                else:
                    G.add_edge(a["TF"],a["TG"],color="r")


            fig=plt.figure(figsize=(10,10))
            edges = G.edges()
            colors = [G[u][v]['color'] for u,v in edges]
            #nx.draw(G, edge_color=colors, with_labels=True, font_weight='bold')

            #legend stuff
            _c = 'rb' 
            clrs = [c for c in _c[:2]]
            pos=nx.spring_layout(G)


            h1 = nx.draw_networkx_nodes(G, pos=pos, node_color = 'grey',
                                    alpha = 0.9, node_size = 300, linewidths=1)


            h2 = nx.draw_networkx_edges(G, pos=pos, width=1, edge_color=colors)


            h3 = nx.draw_networkx_labels(G, pos=pos, font_size=8, font_color='k',font_weight='bold')


            def make_proxy(clr, mappable, **kwargs):
                return Line2D([0, 1], [0, 1], color=clr, **kwargs)
            proxies = [make_proxy(clr, h2, lw=5) for clr in clrs]
            labels=["repression","activation"]
            #end legend stuff



            plt.title(i,fontsize=20)
            plt.legend(["a plot"])
            plt.legend(proxies,labels,fontsize=20)
            plt.show()
    else:
        df=grn
        if only_TFs==True:
            df=df.loc[df["TG"].isin(tfs)]

        if thresh != None:
            df=df.loc[df["zscore"]>thresh]

        interactions=[]
        for k in df["Correlation"]:
            if k>0:
                interactions.append("activation")
            else:
                interactions.append("repression")
        df["interactions"]=interactions


        G=nx.Graph()
        for j in df.index:
            a=df.loc[j]
            if a["interactions"]=="activation":
                G.add_edge(a["TF"],a["TG"],color="b")
            else:
                G.add_edge(a["TF"],a["TG"],color="r")


        fig=plt.figure(figsize=(10,10))
        edges = G.edges()
        colors = [G[u][v]['color'] for u,v in edges]
        #nx.draw(G, edge_color=colors, with_labels=True, font_weight='bold')

        #legend stuff
        _c = 'rb' 
        clrs = [c for c in _c[:2]]
        pos=nx.spring_layout(G)


        h1 = nx.draw_networkx_nodes(G, pos=pos, node_color = 'grey',
                                alpha = 0.9, node_size = 300, linewidths=1)


        h2 = nx.draw_networkx_edges(G, pos=pos, width=1, edge_color=colors)


        h3 = nx.draw_networkx_labels(G, pos=pos, font_size=8, font_color='k',font_weight='bold')


        def make_proxy(clr, mappable, **kwargs):
            return Line2D([0, 1], [0, 1], color=clr, **kwargs)
        proxies = [make_proxy(clr, h2, lw=5) for clr in clrs]
        labels=["repression","activation"]
        #end legend stuff



        plt.title("Network",fontsize=20)
        plt.legend(["a plot"])
        plt.legend(proxies,labels,fontsize=20)
        plt.show()

# In[ ]:


#plot_top_regulators<-function(grn,gene_ranks,tfs,numTopTFs=5, numTargets=5, only_TFs=TRUE,order=NULL)
#plot_top_regulators(dynamic_grn, gene_rank, mmTFs, only_TFs=FALSE)

def plot_top_regulators(grn,gene_ranks,tfs,numTopTFs=5,numTargets=5, only_TFs=True, order=None):
    #grn=dynamic_grn
    #gene_ranks=gene_rank
    #tfs=list(mmTFs["mmTFs"].values)
    #numTopTFs=5
    #numTargets=5
    #only_TFs=True
    #order=None

    if order != None:
            grn = {k: grn[k] for k in order}

    for i in np.arange(len(list(grn.keys()))):

        epoch=list(grn.keys())[i]
        print(epoch,i)
        df=grn[epoch]

        if only_TFs==True:
            df=df.loc[df["TG"].isin(tfs)]


        interactions=[]
        for k in df["corr"]:
            if k>0:
                interactions.append("activation")
            else:
                interactions.append("repression")
        df["interactions"]=interactions

        rank=gene_ranks[epoch]
        topregs=rank


        topregs=list(rank.loc[rank['is_regulator'] == True].index[0:numTopTFs]) 
        df=df.loc[df['TF'].isin(topregs)]


        topdf=pd.DataFrame(columns=["TG","TF","interactions"])
        for reg in topregs:
            targets=list(df.loc[df["TF"]==reg]["TG"].values)
            rank_targets=rank.loc[targets,]
            rank_targets=rank_targets.sort_values(by=['page_rank'],ascending=False)

            num=numTargets
            if numTargets>len(targets):
                num=len(targets)

            toptargets=list(rank_targets.index)[0:num]


            add=df.loc[df["TF"]==reg][["TG","TF","interactions"]]

            add=add.loc[add['TG'].isin(toptargets)]
            topdf=pd.concat([topdf,add])


        G=nx.Graph()
        for j in topdf.index:
            a=topdf.loc[j]
            if a["interactions"]=="activation":
                G.add_edge(a["TF"],a["TG"],color="b")
            else:
                G.add_edge(a["TF"],a["TG"],color="r")


        fig=plt.figure(figsize=(10,10))
        edges = G.edges()
        colors = [G[u][v]['color'] for u,v in edges]
        #nx.draw(G, edge_color=colors, with_labels=True, font_weight='bold')

        #legend stuff
        _c = 'rb' 
        clrs = [c for c in _c[:2]]
        pos=nx.spring_layout(G)


        h1 = nx.draw_networkx_nodes(G, pos=pos, node_color = 'grey', #nodelist=list(topdf["TF"].values),
                                alpha = 0.9, node_size = 300, linewidths=1)
        #h11 = nx.draw_networkx_nodes(G, pos=pos, node_color = 'grey', nodelist=list(topdf["TG"].values),
                                #alpha = 0.9, node_size = 300, linewidths=1)


        h2 = nx.draw_networkx_edges(G, pos=pos, width=1, edge_color=colors)


        h3 = nx.draw_networkx_labels(G, pos=pos, font_size=8, font_color='k',font_weight='bold')


        def make_proxy(clr, mappable, **kwargs):
            return Line2D([0, 1], [0, 1], color=clr, **kwargs)
        proxies = [make_proxy(clr, h2, lw=5) for clr in clrs]
        labels=["repression","activation"]
        #end legend stuff

        plt.title(epoch,fontsize=20)
        plt.legend(["a plot"])
        plt.legend(proxies,labels,fontsize=20)
        plt.show()
    


# In[ ]:


#plot_targets_with_top_regulators<-function(grn,targets,weight_column="zscore",gene_ranks=NULL,numTopRegulators=5,order=NULL){
#plot_targets_with_top_regulators(dynamic_grn,interesting_targets,weight_column="zscore")


def plot_targets_with_top_regulators(grn,targets,weight_column="zscore",gene_ranks=None,numTopRegulators=5,order=None):
    #grn=dynamic_grn
    #targets=interesting_targets
    #weight_column="zscore"
    #gene_ranks=gene_rank
    #numTopRegulators=5
    #order=None

    if order != None:
        grn = {k: grn[k] for k in order}

    for i in np.arange(len(list(grn.keys()))):
        #i=np.arange(len(list(grn.keys())))[0]
        epoch=list(grn.keys())[i]
        print(epoch,i)
        df=grn[epoch]

        #look for targets in epoch GRN
        tgs=list(np.unique(list(df.loc[df["TG"].isin(targets)]["TG"].values)))

        interactions=[]
        #colors=[]
        for k in df["corr"]:
            if k>0:
                interactions.append("activation")
                #colors.append('b')
            else:
                interactions.append("repression")
                #colors.append('r')
        df["interactions"]=interactions


        #find top regulators for each target
        edges_to_keep=pd.DataFrame(columns=["TF","TG","interactions"])

        if weight_column=="page_rank":
            print("pagerank")
            for tg in tgs:
                #tg=tgs[0]

                if gene_ranks==None:
                    sys.exit("Need to supply gene ranks")

                rank=gene_ranks[epoch]
                regs_of_target=list(df.loc[df["TG"]==tg]["TF"])
                rank_regs=rank.loc[regs_of_target]

                rank_regs=rank_regs.sort_values(by=['page_rank'],ascending=False)

                top_regs=list(rank_regs.index)[0:numTopRegulators]

                edges=df.loc[df["TG"]==tg]
                edges=edges.loc[edges["TF"].isin(top_regs)][["TF","TG","interactions"]]

                edges_to_keep=pd.concat([edges_to_keep,edges])
        else:
            for tg in tgs:

                edges=df.loc[df["TG"]==tg]
                edges=edges.sort_values(by=["zscore"],ascending=False)
                edges=edges.iloc[0:numTopRegulators][["TG","TF","interactions"]]
                edges_to_keep=pd.concat([edges_to_keep,edges])

        G=nx.Graph()
        for j in edges_to_keep.index:
            a=edges_to_keep.loc[j]
            if a["interactions"]=="activation":
                G.add_edge(a["TF"],a["TG"],color="b")
            else:
                G.add_edge(a["TF"],a["TG"],color="r")

        G.add_nodes_from(edges_to_keep["TF"],s="^")
        G.add_nodes_from(edges_to_keep["TG"],s="o")

        fig,ax=plt.subplots(figsize=(10,10))
        edges = G.edges()
        colors = [G[u][v]['color'] for u,v in edges]
        #nx.draw(G, edge_color=colors, with_labels=True, font_weight='bold')

        #legend stuff
        _c = 'rb' 
        clrs = [c for c in _c[:2]]

        pos=nx.spring_layout(G)

        nodeShapes = set((aShape[1]["s"] for aShape in G.nodes(data = True)))

        target_labels=["regulator","target"]
        #For each node class...
        count=0
        for aShape in nodeShapes:
            #...filter and draw the subset of nodes with the same symbol in the positions that are now known through the use of the layout.
            h1 = nx.draw_networkx_nodes(G, pos=pos, node_shape=aShape,nodelist = [sNode[0] for sNode in filter(lambda x: x[1]["s"]==aShape,G.nodes(data = True))], node_color = 'grey',
                                    label=target_labels[count],alpha = 0.9, node_size = 300, linewidths=1)
            count=count+1

        h2 = nx.draw_networkx_edges(G, pos=pos, width=1, edge_color=colors)

        h3 = nx.draw_networkx_labels(G, pos=pos, font_size=8, font_color='k',font_weight='bold')


        def make_proxy(clr, mappable, **kwargs):
            return Line2D([0, 1], [0, 1], color=clr, **kwargs)
        proxies = [make_proxy(clr, h2, lw=5) for clr in clrs]
        labels=["repression","activation"]
        #end legend stuff

        plt.title(epoch,fontsize=20)
        leg= plt.legend(scatterpoints = 1,fontsize=20, loc=(1.03,0))
        ax.add_artist(leg)
        plt.legend(["a plot"])
        plt.legend(proxies,labels,fontsize=20, loc=(1.03,0.5))
        plt.show()


# In[ ]:


def hm_dyn(expDat,dynRes,topX=25,cRow=False,cCol=False,limits=[0,10],toScale=False,fontsize_row=4,geneAnn=False):
    #expDat=expSmoothed
    #dynRes=dynTFs
    #topX=25
    #cRow=False
    #cCol=False
    #limits=[0,10]
    #toScale=False
    #fontsize_row=4
    #geneAnn=False

    #topX=100
    allgenes=expDat.index

    sampTab=dynRes[1]
    t1=pd.DataFrame(sampTab["pseudotime"])
    t1.index=sampTab["cell_name"]
    grps=pd.DataFrame(sampTab["group"])
    grps.index=sampTab["cell_name"]

    ord1=t1.sort_values(by='pseudotime', ascending=True)

    expDat=expDat[ord1.index]
    grps=grps.loc[ord1.index]

    genes=dynRes[0]
    genes=list(genes.sort_values(by="expression",ascending=True).iloc[0:topX].index)


    missingGenes=set(genes).difference(set(allgenes))
    if len(missingGenes)>0:
        print("Missing genes:"+str(missingGenes))
        genes=list(set(genes).intersection(set(allgenes)))


    peakTime=pd.DataFrame(expDat.loc[genes].apply(np.argmax,axis=1))
    #print(peakTime)
    peakTime.columns=["peakTime"]
    genesOrdered=peakTime.sort_values(by="peakTime",ascending=True).index

    value=expDat.loc[genesOrdered]

    if toScale==True:
        save_genes=value.index
        save_cells=value.columns
        value=pd.DataFrame(preprocessing.scale(value))
        value.index=save_genes
        value.columns=save_cells
    value=value[value>limits[0]]
    value=value.fillna(limits[0])
    value=value[value<limits[1]]
    value=value.fillna(limits[1])
    groupNames=np.unique(list(grps.index))
    cells=list(groupNames)

    sns.set_theme()
    sns.set(font_scale=.5)
    fig, ax = plt.subplots(figsize=(10,10))         
    ax = sns.heatmap(value,yticklabels=True,xticklabels=False)
    return ax