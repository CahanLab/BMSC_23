# Prepare the envrionment
python3
import numpy as np
import pandas as pd
import scanpy as sc
import episcanpy as epi

# Import the file
# Please download the GSE234540_BMSC_CH data from GEO before start
adata = sc.read_10x_mtx('filtered_feature_bc_matrix', gex_only = False)
gex_rows = list(map(lambda x: x == 'Gene Expression', adata.var['feature_types']))
atac_rows = list(map(lambda x: x == 'Peaks', adata.var['feature_types']))
adata_gem = adata[:, gex_rows].copy()
adata_atac = adata[:, atac_rows].copy()

#Make variables unique
adata_gem.var_names_make_unique()

#Mito_genes and Ribo_genes filter
sc.pp.filter_genes(adata_gem, min_cells=3)
sc.pp.filter_cells(adata_gem, min_genes=100)
adata_gem.obs['n_counts'] = np.sum(adata_gem.X, axis=1).A1

mito_genes = [name for name in adata_gem.var_names if name.startswith('mt-')]
adata_gem.obs['percent_mito'] = np.sum(adata_gem[:, mito_genes].X, axis=1) / np.sum(adata_gem.X,axis=1)
ribo_genes = [name for name in adata_gem.var_names if name.startswith('Rps') or name.startswith('Rpl')]
adata_gem.obs['percent_ribo'] = np.sum(adata_gem[:,ribo_genes].X,axis=1) / np.sum(adata_gem.X,axis=1)
adata_gem.obs['n_counts'] = adata_gem.X.sum(axis=1)

sc.pl.violin(adata_gem, ['n_genes', 'n_counts'], jitter=0.4, multi_panel=True)

sc.pl.violin(adata_gem, ['percent_mito', 'percent_ribo'], jitter=0.4, multi_panel=True)

# Filter out low quality cells
adata_gem= adata_gem[adata_gem.obs['n_counts'] < 30000, :]
adata_gem= adata_gem[adata_gem.obs['n_counts'] > 1000, :]
# adata_gem= adata_gem[adata_gem.obs['percent_mito'] < 0.85, :]

# Remove mitochondria genes and ribosomal genes
non_mito_genes_list = [name for name in adata_gem.var_names if not name.startswith('mt-')]
adata_gem= adata_gem[:, non_mito_genes_list]
non_ribo_genes_list = [name for name in adata_gem.var_names if not name.startswith('Rps')]
adata_gem= adata_gem[:, non_ribo_genes_list]
non_ribo_genes_list2 = [name for name in adata_gem.var_names if not name.startswith('Rpl')]
adata_gem= adata_gem[:, non_ribo_genes_list2]

# Normalization
sc.pp.normalize_total(adata_gem, target_sum=1e4)
sc.pp.log1p(adata_gem)

# Save the original count
adata_raw = adata_gem.copy()

# Find and keep only the most variable genes
sc.pp.highly_variable_genes(adata_gem, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_gem)

adata_gem= adata_gem[:, adata_gem.var['highly_variable']]

#Find PCA features and compute embedding
sc.tl.pca(adata_gem,use_highly_variable=True, n_comps=60, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_gem, 60, log=True)

adata_pre = adata_gem.copy()

adata_gem = adata_pre.copy()
sc.pp.neighbors(adata_gem, n_neighbors=4,n_pcs=30)
sc.tl.leiden(adata_gem, resolution = 0.4)
sc.tl.umap(adata_gem)
adata_all = adata_raw.copy()
adata_all.obsm = adata_gem.obsm
adata_all.obs= adata_gem.obs
sc.pl.umap(adata_all, color=['leiden'])

#assignment
#adata_all.obs['barcode2'] = [x.split('-')[0] for x in adata_all.obs_names.tolist()]
#adata_all.obs.index = adata_all.obs['barcode2']
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/tag_assignment_final_040423.txt", "tag_assignment_final_040423.txt")
assignment = pd.read_csv('tag_assignment_manual_040423.txt', sep = '\t', header = 0)
assignment.index = assignment.index.astype(str) + '-1'
adata_all.obs['multitag'] = assignment['multitag']
sc.pl.umap(adata_all, color=['multitag'])

sc.pl.umap(adata_all, color=['Ptprc', 'Cdh5'])

sc.pl.umap(adata_all, color=['Sox9', 'Col2a1', 'Acan', 'Cxcl12', 'Col1a1', 'Spp1'])

# Macrophage markers
sc.pl.umap(adata_all, color=['Itgam', 'Cd14', 'Cd68'])

sc.tl.rank_genes_groups(adata_all, 'leiden', n_genes=50, method='t-test', use_raw=False)
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(20)

#Keep cluster 0, 1, 2, 3, 5

# Save the File
adata_gem.write_loom('CH_pre_clustered_GEX_042423.loom', write_obsm_varm=True)
adata_all.write_loom('CH_clustered_GEX_042423.loom', write_obsm_varm=True)

a= adata_all[adata_all.obs['leiden'].isin(['0', '1', '2', '3', '5'])].obs.index

# Import the file
adata_gem = adata[:, gex_rows].copy()

#Make variables unique
adata_gem.var_names_make_unique()
adata_gem= adata_gem[adata_gem.obs.index.isin(a)]
adata_gem

#Mito_genes and Ribo_genes filter
sc.pp.filter_genes(adata_gem, min_cells=3)
sc.pp.filter_cells(adata_gem, min_genes=0)
adata_gem.obs['n_counts'] = np.sum(adata_gem.X, axis=1).A1

mito_genes = [name for name in adata_gem.var_names if name.startswith('mt-')]
adata_gem.obs['percent_mito'] = np.sum(adata_gem[:, mito_genes].X, axis=1) / np.sum(adata_gem.X,axis=1)
ribo_genes = [name for name in adata_gem.var_names if name.startswith('Rps') or name.startswith('Rpl')]
adata_gem.obs['percent_ribo'] = np.sum(adata_gem[:,ribo_genes].X,axis=1) / np.sum(adata_gem.X,axis=1)

# Normalization
sc.pp.normalize_per_cell(adata_gem, counts_per_cell_after=1e4)
sc.pp.log1p(adata_gem)

adata_raw = adata_gem.copy()

sc.pp.highly_variable_genes(adata_gem, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_gem)

adata_gem= adata_gem[:, adata_gem.var['highly_variable']]

# Find PCA features and compute embedding
sc.tl.pca(adata_gem,use_highly_variable=True, n_comps=60, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_gem, 60, log=True)

adata_pre = adata_gem.copy()

adata_gem= adata_pre.copy()
sc.pp.neighbors(adata_gem, n_neighbors=50, n_pcs=30)
sc.tl.leiden(adata_gem, resolution = 0.5)
sc.tl.umap(adata_gem)
adata_all = adata_raw.copy()
adata_all.obsm = adata_gem.obsm
adata_all.obs= adata_gem.obs
sc.pl.umap(adata_all, color=['leiden'])


sc.tl.rank_genes_groups(adata_all, 'leiden', n_genes=50, method='t-test_overestim_var', use_raw=False)
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(20)

non_mito_genes_list = [name for name in adata_all.var_names if not name.startswith('mt-')]
adata_all = adata_all[:, non_mito_genes_list]
non_ribo_genes_list = [name for name in adata_all.var_names if not name.startswith('Rps')]
adata_all = adata_all[:, non_ribo_genes_list]
non_ribo_genes_list2 = [name for name in adata_all.var_names if not name.startswith('Rpl')]
adata_all = adata_all[:, non_ribo_genes_list2]
 
sc.tl.rank_genes_groups(adata_all, 'leiden', n_genes=200, use_raw=False)
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(20)

#assignment
assignment = pd.read_csv('tag_assignment_final_042423.txt', sep = '\t', header = 0)
assignment.index = assignment.index.astype(str) + '-1'
adata_all.obs['multitag'] = assignment['multitag']

sc.pl.umap(adata_all, color=['multitag'])

#Diffusion Map and Pseudo Time
adata_diff = adata_all.copy()
sc.pp.neighbors(adata_diff, n_neighbors=25,n_pcs=25)
sc.tl.diffmap(adata_diff , n_comps=25)
sc.pl.diffmap(adata_diff , color = ['multitag'], cmap='OrRd')

adata_diff .uns['iroot'] = np.argmin(adata_diff .obsm['X_umap'][:,1])
sc.pp.neighbors(adata_diff , n_neighbors=25,n_pcs=25)
sc.tl.dpt(adata_diff , n_dcs=4, n_branchings=1)
sc.pl.umap(adata_diff , color = ['multitag', 'dpt_groups', 'dpt_pseudotime'], cmap='OrRd')

adata_diff.obs['dpt_groups_ori'] = adata_diff.obs['dpt_groups']
adata_diff.obs['dpt_groups'] = adata_all.obs['leiden2']

# Save the File
adata_gem.write_loom('CH_pre_clustered_GEX_cleaned_042423.loom', write_obsm_varm=True)
adata_all.write_loom('CH_clustered_GEX_cleaned_042423.loom', write_obsm_varm=True)
