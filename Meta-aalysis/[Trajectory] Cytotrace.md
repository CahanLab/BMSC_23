
```R
# Cytotrace analysis in R
R
library(singleCellNet)
library(dplyr)
library(loomR)
library(CytoTRACE)

lfile <- connect(filename = "/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/ALL_BMSC/All_BMSC_bone_clustered_bc_20230329.loom", skip.validate = TRUE)
geneNames<-lfile[["row_attrs"]][["var_names"]][]
cellNames<-lfile[["col_attrs"]][["obs_names"]][]
expPark<- t(lfile[["matrix"]][1:length(cellNames),])
rownames(expPark)<-geneNames
colnames(expPark)<-cellNames

leiden<-lfile[["col_attrs"]][["type"]][]
stPark<- data.frame(row.names = cellNames, cell_name=cellNames, cluster=leiden)
stPark2 <- stPark$cluster
stPark2 <- as.character(stPark2)
names(stPark2) <- stPark$cell_name

expPark <- as.data.frame(expPark)
results <- CytoTRACE(expPark)

results_df <- as.data.frame(results$CytoTRACE)
write.csv(results_df, file = "BMSC_bone_cytotrace_20230330.csv")
```

```python
# Visualization in Python
python3
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib import rcParams
import loompy
import scvelo as scv

cyto = pd.read_csv("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/ALL_BMSC/Cytotrace/BMSC_bone_cytotrace_20230329.csv")

cyto.index = adata_bone.obs.index
adata_bone.obs["CytoTRACE"] = ""
adata_bone.obs["CytoTRACE"] = cyto["results$CytoTRACE"]

sc.pl.umap(adata_bone, color = 'CytoTRACE')
```

