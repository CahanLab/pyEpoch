##### Functions to integrate signaling pathways with reconstructed GRNs

# required modules
import os
import glob
import pandas as pd
import numpy as np
import igraph
import itertools

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

# refine effector targets to only include those in adata
def target_overlap(adata,effectortargets_dict):
	effector_targets={}
	for effector in effectortargets_dict:
		overlap = set(adata.var_names).intersection(effectortargets_dict[effector])
		effector_targets[effector] = list(overlap)
	return effector_targets

# Function to compute mean signaling activity (or any module) over pseudotime
# result_name -- what it will be stored as in adata.uns
def mean_module_expression(adata,module_dict,result_name):
	mean_dict={'cell':adata.obs_names}
	for module in module_dict:
		subset=adata[:,module_dict[module]]
		df=pd.DataFrame(subset.X,index=subset.obs_names,columns=subset.var_names)
		mean_activity=df.mean(axis=1).to_frame(name=module).to_dict('list')
		mean_dict.update(mean_activity)
	mean_df = pd.DataFrame.from_dict(mean_dict)
	mean_df = mean_df.set_index('cell')
	adata.uns[result_name]=mean_df
	print("Mean module expression stored in adta.uns['"+result_name+"']")
	return(adata)


# function to fine shortest path from a TF to a target in a static network
# compare_to_average whether or not to normalize by avreage path length
# grn_name -- name of static GRN stored in .uns. Defaults to grnDF.
def static_shortest_path(adata,TF,target,grn_name='grnDF',weight_column="weighted_score",compare_to_average=True):
	grnDF = adata.uns[grn_name]
	# compute relative edge lengths
	grnDF['normalized_score'] = grnDF[weight_column]/grnDF[weight_column].max()
	grnDF['edge_length'] = 1-grnDF['normalized_score']
	ig=igraph.Graph.DataFrame(grnDF[['TF','TG','edge_length','corr']],directed=True)
	if ((TF in list(grnDF['TF']))==False):
		print(TF + " not a TF in the network.")
		vpath=[[]]
	elif ((target in list(grnDF['TG']))==False):
		print(target + " not a TG in the network.")
		vpath=[[]]
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
		distances=np.array(ig.shortest_paths(weights='edge_length',mode='out'))
		# ignore all 0's and infinite lengths
		distances[distances==np.inf]=0
		average_path_length = np.mean(distances,where=(distances>0))
		# return vtx path, distance, normalized distance
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
	return static_shortest_path(new_adata,TF,target,grn_name='grn_for_shortest_path',weight_column=weight_column,compare_to_average=compare_to_average)


# dynamic shortest path for multiple TFs and targets
# this can be sped up by computing average path length once rather than every time calling dynamic_shortest_path/static_shortest_path
def dynamic_shortest_path_multiple(adata, TFs, targets, grn_name='dynamic_GRN', weight_column="weighted_score"):
	res = pd.DataFrame(list(itertools.product(*[TFs,targets])),columns=['source','target'])
	add = res.apply(lambda row: pd.Series(dynamic_shortest_path(adata,row.source,row.target,grn_name=grn_name,weight_column=weight_column,compare_to_average=True)), axis=1)
	res = pd.concat([res,add],axis=1)
	res.columns = ['from','to','path','distance','distance_over_average']
	return res













