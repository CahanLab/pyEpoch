##### Functions to integrate signaling pathways with reconstructed GRNs

# required modules
import os
import glob
import pandas as pd
import numpy as np
import igraph
import itertools
from .utils import *

# Function to load in signaling effectors
# path the path to the folder containing binding data tables
# extension the data table file types
def load_SP_effectors(path,extension='.tsv'):
    files = glob.glob(os.path.join(path,'*{}'.format(extension)))
    binding_dict = {}
    for f in files:
        if (extension=='.tsv'):
            df = pd.read_csv(f,sep='\t',header=0)
        else:
            df = pd.read_csv(f,sep=',',header=0)
        name=f.split("/")[-1]
        print(name)
        name=name.split(extension)[0]
        binding_dict[name] = df
    return binding_dict

# Function to score targets
def score_targets(binding_dict,gene_column='Target_genes',threshold=50):
    for effector in binding_dict:
        df = binding_dict[effector]
        # remove string score
        df = df.drop(['STRING'],axis=1)
        # set index names
        df = df.set_index(gene_column,drop=False)
        # add mean score column (this is different from existing "average" column)
        df = df.drop(df.filter(like="Average").columns,axis=1)
        df['mean_score'] = df.mean(numeric_only=True,axis=1)
        # add max score column
        df['max_score'] = df.max(numeric_only=True,axis=1)
        # add percent freq
        df['percent_freq'] = (df.drop(['Target_genes','mean_score','max_score'],axis=1).gt(threshold).sum(axis=1))/len(df.drop(['Target_genes','mean_score','max_score'],axis=1).columns)
        # up-date binding_dict
        binding_dict[effector] = df
    return binding_dict

# Function to find targets
def find_targets(binding_dict,column="max_score",by_rank=True,n_targets=2000,threshold=50):
    effector_targets={}
    for effector in binding_dict:
        df = binding_dict[effector]
        if (by_rank):
            df = df.sort_values(by=[column],ascending=False)
            targets = df.head(n_targets).index.values
        else:
            targets = df[df[column]>threshold].index.values
        effector_targets[effector] = targets
    return effector_targets

# refine effector targets to only include those in adata  这个看起来有点不对，adata在哪里用的？
def target_overlap(adata,effectortargets_dict):
    effector_targets={}
    for effector in effectortargets_dict:
        overlap = set(adata.var_names).intersection(effectortargets_dict[effector])
        effector_targets[effector] = list(overlap)
    return effector_targets

# Function to compute mean signaling activity (or any module) over pseudotime
# result_name -- what it will be stored as in adata.uns
def mean_module_expression(adata,module_dict,condition=None,condition_column=None,result_name="signaling_activity"):
    if condition is not None:
        adata_sub = adata[adata.uns['cells'][condition_column]== condition,:]
        res = adata.uns['cells'][adata.uns['cells'][condition_column]==condition,:]
    else:
        adata_sub = adata
        res = adata.uns['cells']
    # mean_dict={'cell':adata_sub.obs_names}

    if scipy.sparse.issparse(adata_sub.X):
        expX = pd.DataFrame(adata_sub.X.todense())
    else:
        expX=pd.DataFrame(adata_sub.X)
    
    expX.index=adata_sub.obs_names
    expX.columns=adata_sub.var_names

    # order the rows to match
    expX = expX.reindex(res.index.tolist())

    for module in module_dict.keys():
        genes=adata_sub.var_names.intersection(module_dict[module])
        res[module] = expX[genes].mean(axis=1)

    adata.uns[result_name]=res
    print("Mean module expression stored in adta.uns['"+result_name+"']")

    return(adata)


# function to fine shortest path from a TF to a target in a static network
# compare_to_average whether or not to normalize by avreage path length
# grn_name -- name of static GRN stored in .uns. Defaults to grnDF.
def static_shortest_path(adata,TF,target,grn_name='grnDF',weight_column="weighted_score",compare_to_average=True,quiet=False):
	grnDF = adata.uns[grn_name]
	# compute relative edge lengths
	grnDF['normalized_score'] = grnDF[weight_column]/grnDF[weight_column].max()
	grnDF['edge_length'] = 1-grnDF['normalized_score']
	ig=igraph.Graph.DataFrame(grnDF[['TF','TG','edge_length','corr']],directed=True)
	if (TF not in list(grnDF['TF'])):
		vpath=[[]]
		if not quiet:
			print(TF + " not a TF in the network.")
	elif (target not in list(grnDF['TG'])):
		vpath=[[]]
		if not quiet:
			print(target + " not a TG in the network.")
	else:
		vpath = ig.get_shortest_paths(TF,target,weights='edge_length',mode='out',output='vpath')
	# compute distance
	vtcs = vpath[0]
	dists=[]
	if (len(vtcs)>1):
		for i in range(0,len(vtcs)-1):
			dist = ig.es.select(_source=vtcs[i],_target=vtcs[i+1])['edge_length']
			dists.append(dist[0])

	# compute average distance
	if (compare_to_average==True):
		# distances=np.array(ig.shortest_paths(weights='edge_length',mode='out'))
		# # ignore all 0's and infinite lengths
		# distances[distances==np.inf]=0
		# average_path_length = np.mean(distances,where=(distances>0))
		# # return vtx path, distance, normalized distance
		average_path_length = ig.average_path_length()
		return [ig.vs(vpath[0])['name'], sum(dists), sum(dists)/average_path_length]
	else:
		return [ig.vs(vpath[0])['name'], sum(dists)]


# dynamic shortest path
def dynamic_shortest_path(adata,TF,target,grn_name='dynamic_GRN',weight_column="weighted_score",compare_to_average=True):
	# merge. No, this is not the same as using the static network (some edges are excluded based on dynamic activity)
	dynamic_GRN = adata.uns[grn_name]
	grnDF = pd.concat(dynamic_GRN.values(),ignore_index=True)
	grnDF = grnDF.drop_duplicates()
	new_adata = adata
	new_adata.uns['grn_for_shortest_path'] = grnDF
	return static_shortest_path(new_adata,TF,target,grn_name='grn_for_shortest_path',weight_column=weight_column,compare_to_average=compare_to_average,quiet=True)


# dynamic shortest path for multiple TFs and targets
# this can be sped up by computing average path length once rather than every time calling dynamic_shortest_path/static_shortest_path
def dynamic_shortest_path_multiple(adata, TFs, targets, grn_name='dynamic_GRN', weight_column="weighted_score"):
    if "tfs" in adata.uns.keys():
        TFs = list(set(TFs) & set(adata.uns['tfs']))

    if type(adata.uns[grn_name]) is dict:
        conc = pd.concat(adata.uns[grn_name].values(), ignore_index=True)
        TFs = list(set(TFs) & set(conc["TF"]))

        missing_targets = list(set(targets).difference(set(conc["TG"])))
        print(str(missing_targets) + " not in network.")
    else:
        TFs = list(set(TFs) & set(adata.uns[grn_name]["TF"]))

        missing_targets = list(set(targets).difference(set(adata.uns[grn_name]["TG"])))
        print(missing_targets + " not in network.")

    res = pd.DataFrame(list(itertools.product(*[TFs,targets])),columns=['source','target'])
    add = res.apply(lambda row: pd.Series(dynamic_shortest_path(adata,row.source,row.target,grn_name=grn_name,weight_column=weight_column,compare_to_average=True)), axis=1)
    res = pd.concat([res,add],axis=1)
    res.columns = ['from','to','path','distance','distance_over_average']

    # rm rows with empty paths
    res = res[res['path'].map(lambda d: len(d)) > 0]

    # cor_and_add_action
    genes = list(set(res['from']).union(set(res['to'])))
    adata_sub = adata[:,genes]
    if scipy.sparse.issparse(adata_sub.X):
        expX = pd.DataFrame(adata_sub.X.todense())
    else:
        expX = pd.DataFrame(adata_sub.X)
    expX.index=adata_sub.obs_names
    expX.columns=adata_sub.var_names

    corrmat=expX.corr()
    corrmat['from'] = corrmat.index
    corrdf=pd.melt(corrmat, id_vars=["from"])
    corrdf.columns=["from","to","corr"]

    newdf = res.merge(corrdf,on=['from','to'],how="left")

    newdf['action_by_corr']=1
    newdf.loc[newdf['corr'] < 0, 'action_by_corr'] = -1
    newdf = newdf.drop(['corr'],axis=1)

    return newdf







    






