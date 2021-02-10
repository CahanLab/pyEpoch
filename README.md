# pyEpoch

Gene regulatory network reconstruction from scRNA-seq data. This program was translated from R to Python to be compatible with SCAN-PY


## Introduction
Epoch leverages single-cell transcriptomic data, single-cell analysis methods, and graph theoretic approaches to elucidate GRN structure and dynamic activity. 

## Example Walk Thru 0: The Basics

### Set up
```Python
import numpy as np
import pandas as pd
import scanpy as sc
from pygam import GAM, s,l
from scipy import stats
from sklearn.metrics import normalized_mutual_info_score
import scipy.signal as ss
import math

import igraph as igraph
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import statsmodels.stats.multitest as multi;
import researchpy as rp
from bioinfokit.analys import get_data, stat
from scipy import stats
from sklearn import preprocessing
import sys
import skfda
from skfda import FDataGrid
from skfda.misc import kernels
import skfda.preprocessing.smoothing.kernel_smoothers as ks
import seaborn as sns
```
### Load Data

```Python
#create and name data frames
genes=adata.var.index
sampTab=pd.DataFrame(adata.obs)
cells=list(sampTab.index.values)
expDat=pd.DataFrame(adata.X).T
expDat.columns=sampTab.index
expDat.index=genes
expDat=expDat.loc[expDat.sum(axis=1)!=0]

```
### Static Network Reconstruction
Reconstruction occurs in three steps: 

1. Find dynamically expressed genes
2. Infer edges across dynamic genes using CLR (or other supported method)
3. Perform optional cross-weighting to refine network structure
``` Python
#Find Dynamically Expressed Genes
xdyn=findDynGenes(expDat, sampTab, group_column="leiden",pseudotime_column="dpt_pseudotime")
pThresh=0.05
DataFrameGenes=pd.DataFrame(xdyn[0]<pThresh)
dgenes=DataFrameGenes[DataFrameGenes[0]==True].index.values

# Reconstruct and perform optional crossweighting
grnDF=reconstructGRN(expDat.iloc[dgenes,],zThresh=3)
grnDF=crossweight(grnDF,expSmoothed=expDat)
```
The object grnDF contains the reconstructed network. TG and TF refer to target gene and transcription factor respectively. The column "zscore" is the network prior to crossweighting. The column "weighted_score" is the network after crossweighting:

```Python
print(grnDF.iloc[0:5,:])

#     TG      TF    zscore      corr    offset  weighted_score
#0  Eya1    Myog  4.178345 -0.261096  2.365385        4.178345
#1  Eya1   Dmrt2  4.772928  0.213328 -2.346154        4.772928
#2  Eya1    Lbx1  3.556854  0.227854  1.365385        3.556854
#3   Msc    Myog  5.340096 -0.482617  3.500000        5.340096
#4   Msc  Cited1  5.910916  0.274095 -1.750000        5.910916
```

### Dynamic Network Extraction
We can further explore changes in the network across time by defining "epochs" or time periods in our trajectory, assigning genes to these epochs, and extracting a dynamic network across time.  

Defining epochs can be done in a number of ways. Here we show an example with method="pseudotime". This will partition cells based on pseudotime (pseudotime will be divided evenly, unless specified with parameter psuedotime_cuts). Althernatively, we can define epochs by "cell_order", in which cells are partitioned based on raw cell order rather than pseudotime, or "group", in which partitions are pre-defined.  

For a simpler approach, assign_epoch_simple() will define and assign epochs based on maximum mean expression of a gene. This approach assumes genes cannot belong to more than one epoch.


