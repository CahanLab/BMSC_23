python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy as sc
from matplotlib import rcParams
import seaborn as sns
import loompy
import scvelo as scv
import gseapy as gp
from gseapy.plot import barplot, dotplot

# Import data
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/All_BMSC_marrow_clustered_bc_20230329.loom", "All_BMSC_marrow_clustered_bc_20230329.loom")
adata_all = sc.read("All_BMSC_clustered_20230310.loom")

adata_all.obs_names = adata_all.obs['obs_names']
adata_all.var_names = adata_all.var['var_names']

sc.pp.normalize_per_cell(adata_all, counts_per_cell_after=1e4)
sc.pp.log1p(adata_all)

# Check the list of available libraries and select the libraries
names = gp.get_library_name()
gene_sets= ["GO_Biological_Process_2021", "Reactome_2022", "WikiPathway_2021_Human"]

# Comparing Cluster Features
sc.tl.rank_genes_groups(adata_all,'type', n_genes=adata_all.shape[1], use_raw=False)
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(20)

### BMSC1
glsc0 = adata_all.uns['rank_genes_groups']['names']['BMSC1'].tolist()[0:100]
er_Sig_0 = gp.enrichr(gene_list=glsc0, gene_sets=gene_sets, organism = 'mouse', no_plot=True)
er_Sig_0.results.to_csv("BMSC1_Sig_20230310.csv")

### BMSC2
glsc1 = adata_all.uns['rank_genes_groups']['names']['BMSC2'].tolist()[0:100]
er_Sig_1 = gp.enrichr(gene_list=glsc1, gene_sets=gene_sets, organism = 'mouse', no_plot=True)
er_Sig_1.results.to_csv("BMSC2_Sig_20230310.csv")

### BMSC3
glsc2 = adata_all.uns['rank_genes_groups']['names']['BMSC3'].tolist()[0:100]
er_Sig_2 = gp.enrichr(gene_list=glsc2, gene_sets=gene_sets, organism = 'mouse', no_plot=True)
er_Sig_2.results.to_csv("BMSC3_Sig_20230310.csv")

### BMSC4
glsc3 = adata_all.uns['rank_genes_groups']['names']['BMSC4'].tolist()[0:100]
er_Sig_3 = gp.enrichr(gene_list=glsc3, gene_sets=gene_sets, organism = 'mouse', no_plot=True)
er_Sig_3.results.to_csv("BMSC4_Sig_20230310.csv")

### BMSC5
glsc4 = adata_all.uns['rank_genes_groups']['names']['BMSC5'].tolist()[0:100]
er_Sig_4 = gp.enrichr(gene_list=glsc4, gene_sets=gene_sets, organism = 'mouse', no_plot=True)
er_Sig_4.results.to_csv("BMSC5_Sig_20230310.csv")


# Import data
adata_all = sc.read("All_BMSC_marrow_clustered_bc_20230329.loom")

adata_all.obs_names = adata_all.obs['obs_names']
adata_all.var_names = adata_all.var['var_names']

sc.pp.normalize_per_cell(adata_all, counts_per_cell_after=1e4)
sc.pp.log1p(adata_all)

# Check the list of available libraries and select the libraries
names = gp.get_library_name()
gene_sets= ["GO_Biological_Process_2021", "Reactome_2022", "WikiPathway_2021_Human"]

# Comparing Cluster Features
sc.tl.rank_genes_groups(adata_all,'type1', n_genes=adata_all.shape[1], use_raw=False)
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(20)

### BMSC1
glsc0 = adata_all.uns['rank_genes_groups']['names']['BMSC1'].tolist()[0:100]
er_Sig_0 = gp.enrichr(gene_list=glsc0, gene_sets=gene_sets, organism = 'mouse', no_plot=True)
er_Sig_0.results.to_csv("BMSC1_Sig_20230331.csv")

### BMSC2
glsc1 = adata_all.uns['rank_genes_groups']['names']['BMSC2'].tolist()[0:100]
er_Sig_1 = gp.enrichr(gene_list=glsc1, gene_sets=gene_sets, organism = 'mouse', no_plot=True)
er_Sig_1.results.to_csv("BMSC2_Sig_20230331.csv")

### BMSC3
glsc2 = adata_all.uns['rank_genes_groups']['names']['BMSC3'].tolist()[0:100]
er_Sig_2 = gp.enrichr(gene_list=glsc2, gene_sets=gene_sets, organism = 'mouse', no_plot=True)
er_Sig_2.results.to_csv("BMSC3_Sig_20230331.csv")

### BMSC4
glsc3 = adata_all.uns['rank_genes_groups']['names']['BMSC4'].tolist()[0:100]
er_Sig_3 = gp.enrichr(gene_list=glsc3, gene_sets=gene_sets, organism = 'mouse', no_plot=True)
er_Sig_3.results.to_csv("BMSC4_Sig_20230331.csv")

### BMSC5
glsc4 = adata_all.uns['rank_genes_groups']['names']['BMSC5'].tolist()[0:100]
er_Sig_4 = gp.enrichr(gene_list=glsc4, gene_sets=gene_sets, organism = 'mouse', no_plot=True)
er_Sig_4.results.to_csv("BMSC5_Sig_20230331.csv")

### BMSC6
glsc5 = adata_all.uns['rank_genes_groups']['names']['BMSC6'].tolist()[0:100]
er_Sig_5 = gp.enrichr(gene_list=glsc5, gene_sets=gene_sets, organism = 'mouse', no_plot=True)
er_Sig_5.results.to_csv("BMSC6_Sig_20230331.csv")

### BMSC7
glsc6 = adata_all.uns['rank_genes_groups']['names']['BMSC7'].tolist()[0:100]
er_Sig_6 = gp.enrichr(gene_list=glsc6, gene_sets=gene_sets, organism = 'mouse', no_plot=True)
er_Sig_6.results.to_csv("BMSC7_Sig_20230331.csv")

