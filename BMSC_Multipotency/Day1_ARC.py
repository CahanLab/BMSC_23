# Prepare the envrionment
python3
import numpy as np
import pandas as pd
import scanpy as sc
import episcanpy.api as epi

# Import the file
# Please download the GSE234540_BMSC_day1 data from GEO before start
adata = sc.read_10x_mtx('filtered_feature_bc_matrix', gex_only = False)
gex_rows = list(map(lambda x: x == 'Gene Expression', adata.var['feature_types']))
atac_rows = list(map(lambda x: x == 'Peaks', adata.var['feature_types']))
adata_gem = adata[:, gex_rows].copy()
adata_atac = adata[:, atac_rows].copy()

# AnnData object with n_obs × n_vars = 5988 × 141320

#Make variables unique
adata_gem.var_names_make_unique()

# AnnData object with n_obs × n_vars = 5988 × 32285

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
adata_gem= adata_gem[adata_gem.obs['n_counts'] < 20000, :]
adata_gem= adata_gem[adata_gem.obs['percent_mito'] < 0.6, :]

# View of AnnData object with n_obs × n_vars = 3881 × 16049

# Remove mitochondria genes and ribosomal genes
non_mito_genes_list = [name for name in adata_gem.var_names if not name.startswith('mt-')]
adata_gem= adata_gem[:, non_mito_genes_list]
non_ribo_genes_list = [name for name in adata_gem.var_names if not name.startswith('Rps')]
adata_gem= adata_gem[:, non_ribo_genes_list]
non_ribo_genes_list2 = [name for name in adata_gem.var_names if not name.startswith('Rpl')]
adata_gem= adata_gem[:, non_ribo_genes_list2]

# View of AnnData object with n_obs × n_vars = 3881 × 15940

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
sc.pp.neighbors(adata_gem, n_neighbors=20,n_pcs=20)
sc.tl.leiden(adata_gem, resolution = 0.5)
sc.tl.umap(adata_gem)
adata_all = adata_raw.copy()
adata_all.obsm = adata_gem.obsm
adata_all.obs= adata_gem.obs
sc.pl.umap(adata_all, color=['leiden'])



sc.pl.umap(adata_all, color=['Hba-a1', 'Hba-a2', 'Ptprc', 'Cdh5'])


sc.pl.umap(adata_all, color=['Sox6', 'Kit', 'Myh11', 'Cxcl12'])

sc.pl.umap(adata_all, color=['Bglap', 'Spp1'])


sc.tl.rank_genes_groups(adata_all, 'leiden', n_genes=50, method='t-test', use_raw=False)
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(20)

#Keep cluster 0, 2, 3, and 6

# Save the File
adata_gem.write_loom('Day1_pre_clustered_GEX_112022.loom', write_obsm_varm=True)
adata_all.write_loom('Day1_clustered_GEX_112022.loom', write_obsm_varm=True)

a= adata_all[adata_all.obs['leiden'].isin(['0', '2', '3', '6'])].obs.index
# 2592 cells

# Import the file
adata_gem = adata[:, gex_rows].copy()

#Make variables unique
adata_gem.var_names_make_unique()
adata_gem= adata_gem[adata_gem.obs.index.isin(a)]
adata_gem

# View of AnnData object with n_obs × n_vars = 2592 × 32285

#Mito_genes and Ribo_genes filter
sc.pp.filter_genes(adata_gem, min_cells=3)
sc.pp.filter_cells(adata_gem, min_genes=0)
adata_gem.obs['n_counts'] = np.sum(adata_gem.X, axis=1).A1

mito_genes = [name for name in adata_gem.var_names if name.startswith('mt-')]
adata_gem.obs['percent_mito'] = np.sum(adata_gem[:, mito_genes].X, axis=1) / np.sum(adata_gem.X,axis=1)
ribo_genes = [name for name in adata_gem.var_names if name.startswith('Rps') or name.startswith('Rpl')]
adata_gem.obs['percent_ribo'] = np.sum(adata_gem[:,ribo_genes].X,axis=1) / np.sum(adata_gem.X,axis=1)

# Remove mitochondria genes and ribosomal genes
non_mito_genes_list = [name for name in adata_gem.var_names if not name.startswith('mt-')]
adata_gem= adata_gem[:, non_mito_genes_list]
non_ribo_genes_list = [name for name in adata_gem.var_names if not name.startswith('Rps')]
adata_gem= adata_gem[:, non_ribo_genes_list]
non_ribo_genes_list2 = [name for name in adata_gem.var_names if not name.startswith('Rpl')]
adata_gem= adata_gem[:, non_ribo_genes_list2]

# View of AnnData object with n_obs × n_vars = 2592 × 14793

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
sc.pp.neighbors(adata_gem, n_neighbors=10, n_pcs=20)
sc.tl.leiden(adata_gem, resolution = 0.2)
sc.tl.umap(adata_gem)
adata_all = adata_raw.copy()
adata_all.obsm = adata_gem.obsm
adata_all.obs= adata_gem.obs
sc.pl.umap(adata_all, color=['leiden'])

sc.tl.rank_genes_groups(adata_all, 'leiden', n_genes=50, method='t-test_overestim_var', use_raw=False)
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(20)

# Save the File
adata_gem.write_loom('Day1_pre_clustered_GEX_cleaned_112222.loom', write_obsm_varm=True)
adata_all.write_loom('Day1_clustered_GEX_cleaned_112222.loom', write_obsm_varm=True)

adata_all.obs_names= adata_all.obs.index.str[:-2]
adata_all.obs_names.to_frame().to_csv('Day1_filtered_cells_20221122.txt', header=False, index=False)


# ATAC Part:
adata_atac = adata_atac[adata_all.obs_names, :]

# View of AnnData object with n_obs × n_vars = 2805 × 109035
# Add annotation
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/GSE234540_BMSC_day1_atac_peak_annotation.tsv", "atac_peak_annotation.tsv")
ann = pd.read_csv("atac_peak_annotation.tsv", sep='\t', header = 0)

adata_atac.var_names_make_unique()

ann['start'] = ann['start'].astype(str)
ann['end'] = ann['end'].astype(str)
ann['Name1'] = ann['chrom'].str.cat(ann['start'],sep=":")
ann['Name2'] = ann['Name1'].str.cat(ann['end'],sep="-")
ann.index = ann['Name2']

# adata_atac.var['gene'] = ann['gene']
# error will show if there is a duplicate

ann.index.unique()
is_duplicate = ann.index.duplicated(keep="first")
not_duplicate = ~is_duplicate
no_duplicate_indeces = ann[not_duplicate]
adata_atac.var['gene'] = no_duplicate_indeces['gene']

adata_atac.var['peak_type'] = no_duplicate_indeces['peak_type']

# Preprocessing
# Binarize the data
epi.pp.binarize(adata_atac)
print(np.max(adata_atac.X))

epi.pp.filter_cells(adata_atac, min_features=1)
epi.pp.filter_features(adata_atac, min_cells=1)

# AnnData object with n_obs × n_vars = 2805 × 108774

# Quality controls
adata_atac.obs['log_nb_features'] = [np.log10(x) for x in adata_atac.obs['nb_features']]
epi.pl.violin(adata_atac, ['nb_features'])



epi.pl.violin(adata_atac, ['log_nb_features'])



min_features = 500
epi.pp.coverage_cells(adata_atac, binary=True, log=False, bins=50,
               threshold=min_features, save='BMSCMM_bulk_peaks_coverage_cells.png')


epi.pp.coverage_cells(adata_atac, binary=True, log=10, bins=50,
               threshold=min_features, save='BMSCMM_bulk_peaks_coverage_cells_log10.png')



min_cells = 5
epi.pp.coverage_features(adata_atac, binary=True, log=False, 
                        threshold=min_cells, save='BMSCMM_bulk_peaks_coverage_peaks.png')



epi.pp.coverage_features(adata_atac, binary=True, log=True, 
                        threshold=min_cells, save='BMSCMM_bulk_peaks_coverage_peaks_log10.png')




epi.pp.filter_cells(adata_atac, min_features=min_features)
epi.pp.filter_features(adata_atac, min_cells=min_cells)
adata_atac

# AnnData object with n_obs × n_vars = 2303 × 104159

# Filter GEX based on ATAC selection
adata_all = adata_all[adata_atac.obs_names, :]

sc.tl.rank_genes_groups(adata_all, 'leiden', n_genes=50, method='t-test_overestim_var', use_raw=False)
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(20)


epi.pp.coverage_cells(adata_atac, binary=True, log='log10', bins=50, threshold=min_features)

epi.pp.coverage_features(adata_atac, binary=True, log='log10', bins=50, threshold=min_cells)

# Identify most variable features
epi.pp.cal_var(adata_atac)

# Pick the values after the elbow to include as much as variation possible without much noise
min_score_value = 0.515
nb_feature_selected = 60000
epi.pl.variability_features(adata_atac,log=None,
                     min_score=min_score_value, nb_features=nb_feature_selected)

epi.pl.variability_features(adata_atac,log='log10',
                     min_score=min_score_value, nb_features=nb_feature_selected)

adata_atac.raw = adata_atac

adata_atac = epi.pp.select_var_feature(adata_atac,
                                  nb_features=nb_feature_selected,
                                  show=False,
                                  copy=True)

epi.pl.violin(adata_atac, ['nb_features'])

epi.pl.violin(adata_atac, ['log_nb_features'])

# Normalization
adata_atac_raw = adata_atac.copy()

adata_atac  = adata_atac_raw.copy()
epi.pp.lazy(adata_atac)
sc.pl.umap(adata_atac, color=['nb_features'], wspace=0.3)

adata_atac.layers['binary'] = adata_atac.X.copy()

epi.tl.leiden(adata_atac, resolution = 0.4)
epi.pl.umap(adata_atac, color=['leiden'], wspace=0.4)

adata_atac.obs['GEX'] = adata_all.obs['leiden']
epi.pl.umap(adata_atac, color=['GEX'], wspace=0.4)

adata_atac_all = adata_atac_raw.copy()
adata_atac_all.obsm = adata_atac.obsm
adata_atac_all.obs= adata_atac.obs
epi.pl.umap(adata_atac_all, color=['leiden'])

adata_all.obs['ATAC'] = adata_atac.obs['leiden']
sc.pl.umap(adata_all, color = 'ATAC')


epi.tl.rank_features(adata_atac_all, 'leiden', omic='ATAC', use_raw=False)
epi.pl.rank_feat_groups(adata_atac_all, feature_symbols='gene')


adata_all.write_loom('Day1_clustered_GEX_ATAC_filtered_112222.loom', write_obsm_varm=True)
adata_atac_raw.write_loom("Day1_ATAC_112222.loom", write_obsm_varm=True)
adata_atac_all.write_loom("Day1_clustered_ATAC_112222.loom", write_obsm_varm=True)
