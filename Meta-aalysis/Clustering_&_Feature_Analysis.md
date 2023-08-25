# Prepare python environment
```python
python3
import numpy as np
import pandas as pd
import scanpy as sc
import loompy
```

# Import Data
```python
adata = sc.read("the_loom_file_of_sample_of_your_interest", cache=True)
adata.var_names_make_unique()
adata
```

# Define criteria to filter out bad quality cell/potential doublets
```python
def filter_10x(mdata, topX = 95, mCells=3, mGenes=200):
    sc.pp.filter_genes(mdata, min_cells=mCells)
    sc.pp.filter_cells(mdata, min_genes=mGenes)
    mdata.obs['n_counts'] = np.sum(mdata.X, axis=1).A1
    thresh = np.percentile(mdata.obs['n_counts'],topX)
    mdata = mdata[mdata.obs['n_counts'] < thresh, :]
    return mdata

def percent_cal(MSC): 
	MSC= filter_10x(MSC, mCells=3)
	mito_genes = [name for name in MSC.var_names if name.startswith('mt-')]
	MSC.obs['percent_mito'] = np.sum(MSC[:, mito_genes].X, axis=1) / np.sum(MSC.X,axis=1)
	ribo_genes = [name for name in MSC.var_names if name.startswith('Rps') or name.startswith('Rpl')]
	MSC.obs['percent_ribo'] = np.sum(MSC[:,ribo_genes].X,axis=1) / np.sum(MSC.X,axis=1)
	MSC.obs['n_counts'] = MSC.X.sum(axis=1)
	return MSC

adata = percent_cal(adata)
sc.pl.violin(adata , ['n_genes', 'n_counts', 'percent_mito', 'percent_ribo'], jitter=0.4, multi_panel=True)
```

# Filter out bad cells
```python
def filter_prep(adata, ngenes, percentmito, percentribo):
	adata= adata[adata.obs['n_genes'] < ngenes, :]
	adata= adata[adata.obs['percent_mito'] < percentmito, :]
	adata= adata[adata.obs['percent_ribo'] < percentribo, :]
	sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
	sc.pp.log1p(adata)
	return (adata)

adata= filter_prep(adata, 3000, 0.15, 0.2) # the numbers are defined based on the characteristics of each sample
adata

adata_raw = adata.copy()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

def variable_filter(adata):
	adata= adata[:, adata.var['highly_variable']]
	sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
	sc.pp.scale(adata, max_value=10)
	return(adata)

adata = variable_filter(adata)
```

# Regress out cell cycle effect

```python
def cycle_score(adata, species):
	cell_cycle_genes = [x.strip() for x in open("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Reference/regev_lab_cell_cycle_genes.txt")] 
	if species == "mouse":
		cell_cycle_genes = [x.lower() for x in cell_cycle_genes]
		cell_cycle_genes = [x.title() for x in cell_cycle_genes]
		s_genes = cell_cycle_genes[:43]
		g2m_genes = cell_cycle_genes[43:]
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		if cell_cycle_genes == []:
			print("There is no cell cylce genes in the dataset")
		s_genes = [x for x in s_genes if x in adata.var_names]
		g2m_genes = [x for x in g2m_genes if x in adata.var_names]
		sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
		return(adata)
	elif species == "human":
		s_genes = cell_cycle_genes[:43]
		g2m_genes = cell_cycle_genes[43:]
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		if cell_cycle_genes == []:
			print("There is no cell cylce genes in the dataset")
		s_genes = [x for x in s_genes if x in adata.var_names]
		g2m_genes = [x for x in g2m_genes if x in adata.var_names]
		sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
		return(adata)
	else:
		print("please specify the species to be mouse or human")

def oricycle_plot(adata, species):
	cell_cycle_genes = [x.strip() for x in open("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Reference/regev_lab_cell_cycle_genes.txt")] 
	if species == "mouse":
		cell_cycle_genes = [x.lower() for x in cell_cycle_genes]
		cell_cycle_genes = [x.title() for x in cell_cycle_genes]
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		adata_cc_genes = adata[:, cell_cycle_genes]
		sc.tl.pca(adata_cc_genes,  svd_solver='arpack')
		sc.pl.pca_scatter(adata_cc_genes, color='phase')
	elif species == "human":
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		adata_cc_genes = adata[:, cell_cycle_genes]
		sc.tl.pca(adata_cc_genes,  svd_solver='arpack')
		sc.pl.pca_scatter(adata_cc_genes, color='phase')
	else:
		print("please specify the species to be mouse or human")

def cycle_regress(adata):
	sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
	sc.pp.scale(adata)
	return(adata)

adata = variable_filter(adata)
adata = cycle_score(adata, "mouse")
oricycle_plot(adata, "mouse")

adata = cycle_regress(adata)
oricycle_plot(adata, "mouse")
```

# Unsupervised clustering
```python
PCA
sc.tl.pca(adata,use_highly_variable=True, n_comps=60, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, 60, log=True)

adata_pre = adata.copy()
adata= adata_pre.copy()
sc.pp.neighbors(adata, n_neighbors=100,n_pcs=20)
sc.tl.louvain(adata, resolution = 0.4)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['louvain'])

adata_all= adata_raw.copy()
adata_all.obsm['X_pca'] = adata.obsm['X_pca']
adata_all.obsm['X_umap'] = adata.obsm['X_umap']
adata_all.obs['louvain'] = adata.obs['louvain']
sc.pl.umap(adata_all, color=['louvain'])
```

# Differentially expressed gene analysis
```python
sc.tl.rank_genes_groups(adata_all, 'louvain')
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(25)
```
