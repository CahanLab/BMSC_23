# Prepare the envrionment
python3
import numpy as np
import pandas as pd
import scanpy as sc

# RNA Part:
# Import Data
# AD BMSC
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/AD_clustered_GEX_cleaned_021723.loom", "AD_clustered_GEX_cleaned_021723.loom")
adata_all1 = sc.read("AD_clustered_GEX_cleaned_021723.loom")
adata_all1.obs.index = adata_all1.obs['obs_names']
adata_all1.var.index = adata_all1.var['var_names']

# Please download the GSE234540_BMSC_AD data from GEO
adata1 = sc.read_10x_mtx('filtered_feature_bc_matrix')

adata1.var_names_make_unique()
a= adata_all1.obs.index
adata1= adata1[adata1.obs.index.isin(a)]
adata1.obs = adata_all1.obs
adata1.obs['lineage'] = 'AD'

sc.pp.normalize_per_cell(adata1, counts_per_cell_after=1e4)
sc.pp.log1p(adata1)

# Label time
adata1.obs['time'] = ''
conditions = [(adata1.obs['multitag'] == 'Bar1'), (adata1.obs['multitag'] == 'Bar2'), (adata1.obs['multitag'] == 'Bar3'), (adata1.obs['multitag'] == 'Bar4')]
choices = ['AD_Day0', 'AD_Day1','AD_Day3','AD_Day7']
adata1.obs['time'] = np.select(conditions, choices)

# OS BMSC
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/OS_clustered_GEX_cleaned_032823.loom", "OS_clustered_GEX_cleaned_032823.loom")
adata_all2 = sc.read("OS_clustered_GEX_cleaned_032823.loom")
adata_all2.obs.index = adata_all2.obs['obs_names']
adata_all2.var.index = adata_all2.var['var_names']

# Please download the GSE234540_BMSC_OS data from GEO
adata2 = sc.read_10x_mtx('filtered_feature_bc_matrix')

adata2.var_names_make_unique()
b= adata_all2.obs.index
adata2= adata2[adata2.obs.index.isin(b)]
adata2.obs = adata_all2.obs
adata2.obs['lineage'] = 'OS'

sc.pp.normalize_per_cell(adata2, counts_per_cell_after=1e4)
sc.pp.log1p(adata2)

# Label time
adata2.obs['time'] = ''
conditions = [(adata2.obs['multitag'] == 'Bar1'), (adata2.obs['multitag'] == 'Bar2'), (adata2.obs['multitag'] == 'Bar3'), (adata2.obs['multitag'] == 'Bar4')]
choices = ['OS_Day0', 'OS_Day1','OS_Day3','OS_Day7']
adata2.obs['time'] = np.select(conditions, choices)

# CH BMSC
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/CH_clustered_GEX_cleaned_042423.loom", "CH_clustered_GEX_cleaned_042423.loom")
adata_all3 = sc.read("CH_clustered_GEX_cleaned_042423.loom")
adata_all3.obs.index = adata_all3.obs['obs_names']
adata_all3.var.index = adata_all3.var['var_names']

# Please download the GSE234540_BMSC_CH data from GEO
adata3= sc.read_10x_mtx('filtered_feature_bc_matrix')

adata3.var_names_make_unique()
c= adata_all3.obs.index
adata3= adata3[adata3.obs.index.isin(c)]
adata3.obs = adata_all3.obs

adata3.obs['lineage'] = 'CH'

sc.pp.normalize_per_cell(adata3, counts_per_cell_after=1e4)
sc.pp.log1p(adata3)

# Label time
adata3.obs['time'] = ''
conditions = [(adata3.obs['multitag'] == 'Bar1'), (adata3.obs['multitag'] == 'Bar2'), (adata3.obs['multitag'] == 'Bar3'), (adata3.obs['multitag'] == 'Bar4')]
choices = ['CH_Day0', 'CH_Day1','CH_Day3','CH_Day7']
adata3.obs['time'] = np.select(conditions, choices)

adata = adata1.concatenate([adata2, adata3])
#adata= adata[adata.obs['n_counts'] > 1000, :]

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
sc.pp.neighbors(adata, n_neighbors=25,n_pcs=25)
sc.tl.leiden(adata, resolution = 0.6)
sc.tl.umap(adata)
adata_all= adata_raw.copy()
adata_all.obsm['X_pca'] = adata.obsm['X_pca']
adata_all.obsm['X_umap'] = adata.obsm['X_umap']
adata_all.obsm['X_umap'] = adata.obsm['X_umap']
adata_all.obs['leiden'] = adata.obs['leiden']
sc.pl.umap(adata_all, color=['leiden', 'multitag', 'lineage'])


# Check number of cells in each cluster
adata_all[adata_all.obs['leiden']=='0'].shape[0] 
adata_all[adata_all.obs['leiden']=='1'].shape[0] 
adata_all[adata_all.obs['leiden']=='2'].shape[0] 
adata_all[adata_all.obs['leiden']=='3'].shape[0]

adata_AD = adata_all[adata_all.obs['lineage'] == 'AD']
adata_OS = adata_all[adata_all.obs['lineage'] == 'OS']
adata_CH = adata_all[adata_all.obs['lineage'] == 'CH']

# Batch corrected embeding
import scanpy.external as sce
# scanorama
import scanorama
sce.pp.scanorama_integrate(adata_all, 'batch')
sc.pp.neighbors(adata_all, use_rep="X_scanorama", n_neighbors=100,n_pcs=30)
sc.tl.umap(adata_all)
sc.tl.leiden(adata_all, key_added="clusters", resolution = 0.3)
sc.pl.umap(adata_all, color = ['clusters', 'multitag', 'lineage'])

# harmony
sce.pp.harmony_integrate(adata_all, 'batch')
sc.pp.neighbors(adata_all, use_rep="X_pca_harmony", n_neighbors=10,n_pcs=30)
sc.tl.umap(adata_all)
sc.tl.leiden(adata_all, key_added="clusters", resolution = 0.3)
sc.pl.umap(adata_all, color = ['clusters', 'multitag', 'lineage'])



# Identify subclusters
#sc.pp.neighbors(adata_all, n_neighbors=25,n_pcs=25)
sc.pp.neighbors(adata_all, use_rep="X_pca_harmony", n_neighbors=25,n_pcs=25)
sc.tl.leiden(adata_all, resolution = 0.4, restrict_to=('clusters', ['1']), key_added='leiden0')
sc.tl.leiden(adata_all, resolution = 0.3, restrict_to=('leiden0', ['2']), key_added='leiden1')
sc.tl.leiden(adata_all, resolution = 0.3, restrict_to=('leiden1', ['2,0']), key_added='leiden2')
sc.pl.umap(adata_all, color = ['leiden2', 'multitag', 'lineage'])

# PLot stem finder
new_meta = pd.read_csv("/Users/raycheng/Dropbox (CahanLab)/BMSC_identity/Data/Sequencing/BMSC_All_Diff/StemFinder/stemFinder_All_Diff_20230330.csv")
new_meta.index = adata_all.obs.index
adata_all.obs['stemFinder_invert'] = new_meta['stemFinder_invert'].copy()
sc.pl.umap(adata_all, color= 'stemFinder_invert', color_map = 'RdBu')

# Rename the clusters
adata_all.obs['type'] = ''
conditions = [(adata_all.obs['leiden2'] == '0'), (adata_all.obs['leiden2'] == '2,1'), (adata_all.obs['leiden2'] == '2,0,0'), (adata_all.obs['leiden2'] == '2,0,1'), (adata_all.obs['leiden2'] == '2,2'), (adata_all.obs['leiden2'] == '1,0'), (adata_all.obs['leiden2'] == '1,1'), (adata_all.obs['leiden2'] == '1,2')]
choices = ['BMSC1', 'BMSC2', 'BMSC3', 'BMSC4', 'BMSC5', 'BMSC6', 'BMSC7', 'BMSC8']
adata_all.obs['type'] = np.select(conditions, choices, default='0')

adata_all.obs['type2'] = ''
conditions = [(adata_all.obs['leiden2'] == '0'), (adata_all.obs['leiden2'] == '2,1'), (adata_all.obs['leiden2'] == '2,0,0'), (adata_all.obs['leiden2'] == '2,0,1'), (adata_all.obs['leiden2'] == '2,2'), (adata_all.obs['leiden2'] == '1,0'), (adata_all.obs['leiden2'] == '1,1'), (adata_all.obs['leiden2'] == '1,2')]
choices = ['1', '2', '3', '4', '5', '6', '7', '8']
adata_all.obs['type2'] = np.select(conditions, choices, default='0')

non_mito_genes_list = [name for name in adata_all.var_names if not name.startswith('mt-')]
adata_all = adata_all[:, non_mito_genes_list]
non_ribo_genes_list = [name for name in adata_all.var_names if not name.startswith('Rps')]
adata_all = adata_all[:, non_ribo_genes_list]
non_ribo_genes_list2 = [name for name in adata_all.var_names if not name.startswith('Rpl')]
adata_all = adata_all[:, non_ribo_genes_list2]

sc.tl.rank_genes_groups(adata_all, 'type', n_genes=200, use_raw=False)
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(40)


df = pd.DataFrame(adata_all.uns['rank_genes_groups']['names'])
df.to_csv("AD_OS_CH_DE_list_bc_20230425.csv")

genes = ['Acta2', 'Actb', 'Ccn2', 'Tagln', 'Inhba', 'Actg1', 'Actg2',  'Ibsp', 'Mmp13', 'Ogn', 'Gpx3', 'Gas6', 'Ptx3', 'Igfbp5', 'Alpl', 'Cxcl12', 'Lpl', 'Tgfbr3', 'Kitl', 'Negr1', 'Vcam1', 'Vcan', 'Lepr',  'Fabp4', 'Adipoq', 'Lpl', 'Cd36', 'Fabp5', 'Adipor2', 'Pparg', 'Ebf1', 'Malat1', 'Zeb2', 'Runx1', 'Tnc', 'Cald1']

a= adata_all[adata_all.obs['type'].isin(['BMSC1', 'BMSC2', 'BMSC3', 'BMSC4', 'BMSC5'])].obs.index
adata_ad = adata_all[adata_all.obs.index.isin(a)]

sc.tl.rank_genes_groups(adata_ad, 'type', n_genes=adata_ad.shape[1], use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_ad, var_names = genes, values_to_plot='logfoldchanges', min_logfoldchange=1.5, vmax=5, vmin=-5, cmap='bwr', dendrogram = False)

a= adata_all[adata_all.obs['type'].isin(['BMSC6', 'BMSC7', 'BMSC8'])].obs.index
adata_os = adata_all[adata_all.obs.index.isin(a)]

sc.tl.rank_genes_groups(adata_os, 'type', n_genes=200, use_raw=False)
pd.DataFrame(adata_os.uns['rank_genes_groups']['names']).head(40)

genes = ['Acta2', 'Ccn2', 'S100a6', 'Tagln', 'Tpm2', 'Lgals1', 'Actg1', 'Actg2', 'Actb',  'Col11a1', 'Col12a1',  'Fgfr1', 'Edil3', 'Plcb1', 'Fgfr2', 'Exoc4', 'Fndc3b', 'Runx2', 'Tnc', 'Nfia', 'Runx1', 'Zeb2', 'Lifr', 'Cxcl12', 'Mmp13', 'Kitl', 'Negr1', 'Limch1']

adata_os.obs['type2'] = ''
conditions = [(adata_os.obs['type'] == 'BMSC6'), (adata_os.obs['type'] == 'BMSC7'), (adata_os.obs['type'] == 'BMSC8')]
choices = ['BMSC1', 'BMSC3', 'BMSC2']
adata_os.obs['type2'] = np.select(conditions, choices, default='0')


sc.tl.rank_genes_groups(adata_os, 'type', n_genes=adata_os.shape[1], use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_os, var_names = genes, values_to_plot='logfoldchanges', min_logfoldchange=1.5, vmax=5, vmin=-5, cmap='bwr', dendrogram = False)


adata_all.obs['gene_counts'] = adata_all.obs['n_genes']/adata_all.obs['n_counts']
adata_all.obs['count_genes'] = adata_all.obs['n_counts']/adata_all.obs['n_genes']

# Plot for publication
sc.pl.umap(adata_all, color = ['type'], palette = ['#E97E73', '#CA9734', '#9EAF3C', '#7DC874', '#73C9B0', '#64BEE1', '#7FA5F9', '#CE78F4'])
sc.pl.umap(adata_all, color = ['multitag'], palette = 'cool')
sc.pl.umap(adata_all, color = ['lineage'], palette = 'Accent')

Save the File
adata.write_loom('BMSC_AD_OS_CH_pre_clustered_20231204.loom', write_obsm_varm=True)
adata_all.write_loom('BMSC_AD_OS_CH_clustered_bc_20231204.loom', write_obsm_varm=True)


