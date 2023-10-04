# Prepare the envrionment
python3
import numpy as np
import pandas as pd
import scanpy as sc
import episcanpy.api as epi

# RNA Part:
# Import Data
# Day0 BMSC
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day0_clustered_GEX_ATAC_filtered_112022.loom", "Day0_clustered_GEX_ATAC_filtered_112022.loom")
adata_all1 = sc.read("Day0_clustered_GEX_ATAC_filtered_112022.loom")
adata_all1.obs.index = adata_all1.obs['obs_names']
adata_all1.var.index = adata_all1.var['var_names']

# Please download the GSE234540_BMSC_day0_re1 data from GEO
adata1 = sc.read_10x_mtx('filtered_feature_bc_matrix')
adata1.var_names_make_unique()
a= adata_all1.obs.index
adata1= adata1[adata1.obs.index.isin(a)]
adata1.obs = adata_all1.obs
adata1.obs['time'] = 'Day0'

sc.pp.normalize_per_cell(adata1, counts_per_cell_after=1e4)
sc.pp.log1p(adata1)

# Day0 BMSC - Re
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day0_clustered_GEX_ATAC_filtered_030723.loom", "Day0_clustered_GEX_ATAC_filtered_030723.loom")
adata_all4 = sc.read("Day0_clustered_GEX_ATAC_filtered_030723.loom")
adata_all4.obs.index = adata_all4.obs['obs_names']
adata_all4.var.index = adata_all4.var['var_names']
# Please download the GSE234540_BMSC_day0_re2 data from GEO
adata4 = sc.read_10x_mtx('filtered_feature_bc_matrix')

adata4.var_names_make_unique()
d= adata_all4.obs.index
adata4= adata4[adata4.obs.index.isin(d)]
adata4.obs = adata_all4.obs
adata4.obs['time'] = 'Day0'

sc.pp.normalize_per_cell(adata4, counts_per_cell_after=1e4)
sc.pp.log1p(adata4)

# Day1 BMSC
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day1_clustered_GEX_ATAC_filtered_112222.loom", "Day1_clustered_GEX_ATAC_filtered_112222.loom")
adata_all2 = sc.read("Day1_clustered_GEX_ATAC_filtered_112222.loom")
adata_all2.obs.index = adata_all2.obs['obs_names']
adata_all2.var.index = adata_all2.var['var_names']

# Please download the GSE234540_BMSC_day1 data from GEO
adata2 = sc.read_10x_mtx('filtered_feature_bc_matrix')

adata2.var_names_make_unique()
b= adata_all2.obs.index
adata2= adata2[adata2.obs.index.isin(b)]
adata2.obs = adata_all2.obs
adata2.obs['time'] = 'Day1'

sc.pp.normalize_per_cell(adata2, counts_per_cell_after=1e4)
sc.pp.log1p(adata2)

# Day3 BMSC
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day3_clustered_GEX_ATAC_filtered_112122.loom", "Day3_clustered_GEX_ATAC_filtered_112122.loom")
adata_all3 = sc.read("Day3_clustered_GEX_ATAC_filtered_112122.loom")

adata_all3.obs.index = adata_all3.obs['obs_names']
adata_all3.var.index = adata_all3.var['var_names']

# Please download the GSE234540_BMSC_day3 data from GEO
adata3= sc.read_10x_mtx('filtered_feature_bc_matrix')

adata3.var_names_make_unique()
c= adata_all3.obs.index
adata3= adata3[adata3.obs.index.isin(c)]
adata3.obs = adata_all3.obs

# adata3.obs['time'] = 'Day3'
adata_all3.obs['time'] = ''
conditions = [(adata_all3.obs['multitag'] == 'Bar1'), (adata_all3.obs['multitag'] == 'Bar2'), (adata_all3.obs['multitag'] == 'Bar3')]
choices = ['Day0', 'Day1', 'Day3']
adata_all3.obs['time'] = np.select(conditions, choices, default='0')

sc.pp.normalize_per_cell(adata3, counts_per_cell_after=1e4)
sc.pp.log1p(adata3)

adata = adata1.concatenate([adata2, adata3, adata4])
adata= adata[adata.obs['n_counts'] > 1000, :]

# Save the original count
adata.raw = adata
adata_raw = adata.copy()

# Find and keep only the most variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

adata= adata[:, adata.var['highly_variable']]

#Find PCA features and compute embedding
sc.tl.pca(adata,use_highly_variable=True, n_comps=60, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, 60, log=True)

adata_pre = adata.copy()

adata= adata_pre.copy()
sc.pp.neighbors(adata, n_neighbors=300,n_pcs=40)
# sc.tl.leiden(adata, resolution = 0.2)
sc.tl.umap(adata)
adata_all= adata_raw.copy()
adata_all.obsm['X_pca'] = adata.obsm['X_pca']
adata_all.obsm['X_umap'] = adata.obsm['X_umap']
adata_all.obsm['X_umap'] = adata.obsm['X_umap']
adata_all.obs['leiden'] = adata.obs['leiden']
sc.pl.umap(adata_all, color=['leiden', 'time', 'Top2a'])

sc.pl.umap(adata_all, color=['multitag'])

# Check number of cells in each cluster
adata_all[adata_all.obs['leiden']=='0'].shape[0] 
adata_all[adata_all.obs['leiden']=='1'].shape[0] 
adata_all[adata_all.obs['leiden']=='2'].shape[0] 
adata_all[adata_all.obs['leiden']=='3'].shape[0]

# Batch corrected embeding
# scanorama
import scanpy.external as sce
sce.pp.scanorama_integrate(adata_all, 'batch')
sc.pp.neighbors(adata_all, use_rep="X_scanorama", n_neighbors=100,n_pcs=40)
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color = ['leiden0', 'multitag', 'time'])

# harmony
sce.pp.harmony_integrate(adata_all, 'batch')
sc.pp.neighbors(adata_all, use_rep="X_pca_harmony", n_neighbors=300,n_pcs=40)
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color = ['leiden0', 'multitag', 'time'])

sc.tl.leiden(adata_all, key_added="clusters", resolution = 0.5)
sc.pl.umap(adata_all, color = ['clusters', 'multitag', 'time'])

# Identify subclusters
sc.pp.neighbors(adata_all, n_neighbors=25,n_pcs=25)
# sc.tl.leiden(adata_all, resolution = 0.05, restrict_to=('leiden', ['3']), key_added='leiden0')
sc.tl.leiden(adata_all, resolution = 0.35, restrict_to=('clusters', ['2']), key_added='leiden0')
sc.tl.leiden(adata_all, resolution = 0.1, restrict_to=('leiden0', ['1']), key_added='leiden1')
sc.pl.umap(adata_all, color = ['leiden1'])
sc.pl.umap(adata_all, color = ['type'])

sc.pl.umap(adata_all, color = ['time'], palette = 'cool')

# Plot MSC marker
sc.pl.umap(adata_all, color = ['Lepr', 'Ly6a', 'Nes', 'Cspg4'], ncols = 2, color_map = 'Wistia')

sc.pl.umap(adata_all, color = ['Lepr', 'Ly6a', 'Nes', 'Cspg4'], ncols = 2, color_map = 'BuGn')

# PLot stem finder
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/stemFinder_All_BMSC_DAY_20230330.csv", "stemFinder_All_BMSC_DAY_20230330.csv")
new_meta = pd.read_csv("stemFinder_All_BMSC_DAY_20230330.csv")
new_meta.index = adata_all.obs.index
adata_all.obs['stemFinder_invert'] = new_meta['stemFinder_invert'].copy()
sc.pl.umap(adata_all, color= 'stemFinder_invert', color_map = 'RdBu')

adata_all.obs['stemFinder_comp'] = new_meta['stemFinder_comp'].copy()

# Rename the clusters
adata_all.obs['type'] = ''
conditions = [(adata_all.obs['leiden0'] == '0'), (adata_all.obs['leiden0'] == '2,0'), (adata_all.obs['leiden0'] == '2,1'), (adata_all.obs['leiden0'] == '3'), (adata_all.obs['leiden0'] == '1')]
choices = ['BMSC1', 'BMSC2', 'BMSC3', 'BMSC4', 'BMSC5']
adata_all.obs['type'] = np.select(conditions, choices, default='0')

adata_all.obs['type2'] = ''
conditions = [(adata_all.obs['leiden0'] == '0'), (adata_all.obs['leiden0'] == '2,0'), (adata_all.obs['leiden0'] == '2,1'), (adata_all.obs['leiden0'] == '3'), (adata_all.obs['leiden0'] == '1')]
choices = ['1', '2', '3', '4', '5']
adata_all.obs['type2'] = np.select(conditions, choices, default='0')

adata_all.obs['type3'] = ''
conditions = [(adata_all.obs['leiden0'] == '0'), (adata_all.obs['leiden0'] == '2,0'), (adata_all.obs['leiden0'] == '2,1'), (adata_all.obs['leiden0'] == '3'), (adata_all.obs['leiden0'] == '1')]
choices = ['BMSC4', 'BMSC2', 'BMSC3', 'BMSC1', 'BMSC5']
adata_all.obs['type3'] = np.select(conditions, choices, default='0')

adata_all.obs['type4'] = ''
conditions = [
(np.logical_and(adata_all.obs['leiden0'] == '0', adata_all.obs['time'] == 'Day0')),
(np.logical_and(adata_all.obs['leiden0'] == '2,0', adata_all.obs['time'] == 'Day0')), 
(np.logical_and(adata_all.obs['leiden0'] == '2,1' , adata_all.obs['time'] == 'Day0')), 
(np.logical_and(adata_all.obs['leiden0'] == '3' , adata_all.obs['time'] == 'Day0')), 
(np.logical_and(adata_all.obs['leiden0'] == '1' , adata_all.obs['time'] == 'Day0')),
(np.logical_and(adata_all.obs['leiden0'] == '0' , adata_all.obs['time'] == 'Day1')), 
(np.logical_and(adata_all.obs['leiden0'] == '2,0' , adata_all.obs['time'] == 'Day1')), 
(np.logical_and(adata_all.obs['leiden0'] == '2,1' , adata_all.obs['time'] == 'Day1')), 
(np.logical_and(adata_all.obs['leiden0'] == '3' , adata_all.obs['time'] == 'Day1')), 
(np.logical_and(adata_all.obs['leiden0'] == '1' , adata_all.obs['time'] == 'Day1')),
(np.logical_and(adata_all.obs['leiden0'] == '0' , adata_all.obs['time'] == 'Day3')), 
(np.logical_and(adata_all.obs['leiden0'] == '2,0' , adata_all.obs['time'] == 'Day3')), 
(np.logical_and(adata_all.obs['leiden0'] == '2,1' , adata_all.obs['time'] == 'Day3')), 
(np.logical_and(adata_all.obs['leiden0'] == '3' , adata_all.obs['time'] == 'Day3')), 
(np.logical_and(adata_all.obs['leiden0'] == '1' ,  adata_all.obs['time'] == 'Day3'))
]
choices = ['BMSCa', 'BMSCb', 'BMSCc', 'BMSCd', 'BMSCe', 'BMSCf', 'BMSCg', 'BMSCh', 'BMSCi', 'BMSCj', 'BMSCk', 'BMSCl', 'BMSCm', 'BMSCn', 'BMSCo']
adata_all.obs['type4'] = np.select(conditions, choices, default='0')

a= adata_all[adata_all.obs['leiden'].isin(['0', '1'])].obs.index
adata_cluster1 = adata_all[adata_all.obs.index.isin(a)]

adata_cluster1 = adata_all[adata_all.obs['leiden']=='0']
sc.tl.rank_genes_groups(adata_cluster1, 'leiden0', n_genes=adata_cluster1.shape[1], use_raw=False)

non_mito_genes_list = [name for name in adata_all.var_names if not name.startswith('mt-')]
adata_all = adata_all[:, non_mito_genes_list]
non_ribo_genes_list = [name for name in adata_all.var_names if not name.startswith('Rps')]
adata_all = adata_all[:, non_ribo_genes_list]
non_ribo_genes_list2 = [name for name in adata_all.var_names if not name.startswith('Rpl')]
adata_all = adata_all[:, non_ribo_genes_list2]

# gm_genes_list = [name for name in adata_all.var_names if name.startswith('Gm')]

non_gm_genes_list = [name for name in adata_all.var_names if not name.startswith('Gm')]
adata_all = adata_all[:, non_gm_genes_list]

sc.tl.rank_genes_groups(adata_all, 'type', n_genes=200, use_raw=False)
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(20)

a= adata_all[adata_all.obs['leiden'].isin(['3', '4'])].obs.index
adata_day0 = adata_all[adata_all.obs.index.isin(a)]
b= adata_day0 [adata_day0 .obs['batch'].isin(['0', '3'])].obs.index
adata_day0 = adata_day0 [adata_day0 .obs.index.isin(b)]

sc.tl.rank_genes_groups(adata_day0 , 'type', n_genes=adata_day0 .shape[1], use_raw=False)
pd.DataFrame(adata_day0 .uns['rank_genes_groups']['names']).head(20)

sc.pl.rank_genes_groups_dotplot(adata_day0 , n_genes=10, values_to_plot='logfoldchanges', min_logfoldchange=1.5, vmax=7, vmin=-7, cmap='bwr', dendrogram = False)

# Plot DE genes
sc.tl.rank_genes_groups(adata_all, 'type', n_genes=adata_all.shape[1], use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_all, n_genes=10, values_to_plot='logfoldchanges', min_logfoldchange=1.5, vmax=7, vmin=-7, cmap='bwr', dendrogram = False)

genes = ['Cxcl12', 'Ebf3', 'Kitl', 'Lepr', 'Bmp6', 'Tnc', 'Limch1', 'Lifr', 'Ncam1', 'Runx2', 'Col1a1', 'Col1a2', 'Col11a2', 'Col11a1', 'Mef2c', 'Mki67', 'Top2a', 'Hmgb2', 'Sptb', 'Kel', 'Vegfa', 'Pvt1', 'Ngf', 'Ncl', 'Thbd', 'Timp1', 'Actb', 'Spp1', 'Lgals1', 'Pkm']
sc.pl.dotplot(adata_all, genes, groupby='type', vmax=7, vmin=-7, cmap='bwr', dendrogram = False)
sc.pl.rank_genes_groups_dotplot(adata_all, var_names = genes, values_to_plot='logfoldchanges', min_logfoldchange=1.5, vmax=7, vmin=-7, cmap='bwr', dendrogram = False)

genes = ['Cxcl12', 'Ebf3', 'Kitl', 'Lepr', 'Igfbp5', 'Igfbp4', 'Vcam1', 'Limch1', 'Alpl', 'Lifr', 'Ncam1', 'Runx2', 'Col11a1', 'Col11a2', 'Col5a2', 'Col1a1', 'Mef2c', 'Top2a', 'Hmgb2', 'Mki67', 'Spp1', 'Inhba', 'Actb',  'Timp1', 'Pvt1', 'Ngf', 'Glis3']

sc.tl.rank_genes_groups(adata_all, 'type', n_genes=adata_all.shape[1], use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_all, var_names = genes, values_to_plot='logfoldchanges', min_logfoldchange=1.5, vmax=7, vmin=-7, cmap='bwr', dendrogram = False)

genes = ['Tgfbr3', 'Ebf3', 'Bmp6', 'Egfr', 'Kitl', 'Ebf1', 'Limch1', 'Ncam1', 'Alpl', 'Itga9', 'Gli2', 'Foxn3', 'Prex1', 'Col13a1', 'Col11a1', 'Col11a1', 'Mef2c', 'Col5a1', 'Cdk8', 'Lars2', 'Sp4', 'Klf11', 'Brd2', 'Ccne1', 'Glis3', 'Glis1', 'Hivep3',  'Dst', 'Ank3', 'Arhgef3']

sc.tl.rank_genes_groups(adata_all, 'type', n_genes=adata_all.shape[1], use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_all, var_names = genes, values_to_plot='logfoldchanges', min_logfoldchange=1.5, vmax=7, vmin=-7, cmap='bwr', dendrogram = False)

sc.pl.umap(adata_all, color = ['Gdpd2', 'Ftl1', 'Malat1', 'Myl12a', 'Gnas'])

df1 = pd.DataFrame(adata_all.uns['rank_genes_groups']['names'])
df2 = pd.DataFrame(adata_all.uns['rank_genes_groups']['pvals'])
df3 = pd.DataFrame(adata_all.uns['rank_genes_groups']['logfoldchanges'])
df4 = pd.concat([df1, df2], axis=1)
df5 = pd.concat([df4, df3], axis=1)
df6 = df5[['BMSC1', 'BMSC2', 'BMSC3', 'BMSC4', 'BMSC5', 'BMSC6']]
df6.to_csv("All_days_DE_gene_list_20230309.csv")

# Identify expressed genes
res = pd.DataFrame(columns=adata_all.var_names, index=adata_all.obs['leiden0'].cat.categories)                                                                                                 
for clust in adata_all.obs['leiden0'].cat.categories: 
    res.loc[clust] = adata_all[adata_all.obs['leiden0'].isin([clust]),:].X.mean(0)

res = res.transpose()
res.to_csv('Leiden_expression_20230309.csv', header=True, index=True)co

# Create new names for fragment mapping
adata_all.obs['new_name'] = adata_all.obs.index.str[:-2]

# Save the File
adata_all.obs_names= adata_all.obs.index.str[:-2]
adata_all.obs_names.to_frame().to_csv('BMSC_All_days_filtered_cells_20221122.txt', header=False, index=False)

adata.write_loom('BMSC_All_day_pre_clustered_030823.loom', write_obsm_varm=True)
adata_all.write_loom('BMSC_All_day_clustered_bc_092123.loom', write_obsm_varm=True)

# Day0
a= adata_all[adata_all.obs['leiden0'].isin(['0', '3', '2,0', '2,1'])].obs.index
adata_day0 = adata_all[adata_all.obs.index.isin(a)]
b= adata_day0 [adata_day0 .obs['time'].isin(['Day0'])].obs.index
adata_day0 = adata_day0 [adata_day0 .obs.index.isin(b)]

a= adata_all[adata_all.obs['time'].isin(['Day0'])].obs.index
adata_day0 = adata_all[adata_all.obs.index.isin(a)]

del adata_day0.uns
adata_day0.write_loom('BMSC_day0_clustered_051223.loom', write_obsm_varm=True)

# Day1
a= adata_all[adata_all.obs['leiden'].isin(['1', '3'])].obs.index
adata_day0 = adata_all[adata_all.obs.index.isin(a)]
b= adata_day0 [adata_day0 .obs['time'].isin(['Day1'])].obs.index
adata_day0 = adata_day0 [adata_day0 .obs.index.isin(b)]

del adata_day0.uns
adata_day0.write_loom('BMSC_day1_clustered_032923.loom', write_obsm_varm=True)

# Day3
a= adata_all[adata_all.obs['leiden'].isin(['1', '3'])].obs.index
adata_day0 = adata_all[adata_all.obs.index.isin(a)]
b= adata_day0 [adata_day0 .obs['time'].isin(['Day3'])].obs.index
adata_day0 = adata_day0 [adata_day0 .obs.index.isin(b)]

del adata_day0.uns
adata_day0.write_loom('BMSC_day3_clustered_032923.loom', write_obsm_varm=True)

adata_day0 = adata_all[adata_all.obs['batch'] == '0']
adata_day1 = adata_all[adata_all.obs['batch'] == '1']
adata_day3 = adata_all[adata_all.obs['batch'] == '2']
adata_day0re = adata_all[adata_all.obs['batch'] == '3']

adata_day0.obs_names = adata_day0.obs['obs_names']
adata_day1.obs_names = adata_day1.obs['obs_names']
adata_day3.obs_names = adata_day3.obs['obs_names']
adata_day0re.obs_names = adata_day0re.obs['obs_names']

adata_day0.obs.to_csv("Day0_metatable_20230308.csv")
adata_day1.obs.to_csv("Day1_metatable_20230308.csv")
adata_day3.obs.to_csv("Day3_metatable_20230308.csv")
adata_day0re.obs.to_csv("Day0re_metatable_20230308.csv")

