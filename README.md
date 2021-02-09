# pyEpoch

Gene regulatory network reconstruction from scRNA-seq data. This program was translated from R to Python to be compatible with SCAN-PY


## Introduction
Epoch leverages single-cell transcriptomic data, single-cell analysis methods, and graph theoretic approaches to elucidate GRN structure and dynamic activity. 

## Example Walk Thru 0: The Basics

### Set up

```R
list12<-loadDataFromLoom("data/adMuscle_E12_DPT_071919.loom")
expDat<-list12[['expDat']]
sampTab<-list12[['sampTab']] 
expDat<-expDat[rowSums(expDat)!=0,]

mmTFs<-utils_loadObject("data/mmTFs_123019.rda")
mmTFs<-intersect(rownames(expDat),mmTFs)
```
