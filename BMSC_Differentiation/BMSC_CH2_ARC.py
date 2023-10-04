# Prepare the envrionment
python3
import numpy as np
import pandas as pd
import scanpy as sc
import episcanpy as epi

# Import the file
adata = sc.read_10x_mtx('/Users/raycheng/Dropbox (CahanLab)/BMSC_identity/Data/Sequencing/BMSC_MM_CH_20230127/CH_0127_ARC/outs/filtered_feature_bc_matrix', gex_only = False)
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

# View of AnnData object with n_obs × n_vars = 7876 × 17351

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

AnnData object with n_obs × n_vars = 3141 × 17295

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
assignment = pd.read_csv('/Users/raycheng/Dropbox (CahanLab)/BMSC_identity/Data/Sequencing/BMSC_MM_CH_20230127/CellRanger/tag_assignment_manual_040423.txt', sep = '\t', header = 0)
assignment.index = assignment.index.astype(str) + '-1'
adata_all.obs['multitag'] = assignment['multitag']
sc.pl.umap(adata_all, color=['multitag'])

sc.pl.umap(adata_all, color=['Ptprc', 'Cdh5'])



sc.pl.umap(adata_all, color=['Sox9', 'Col2a1', 'Acan', 'Cxcl12', 'Col1a1', 'Spp1'])


# Macrophage markers
sc.pl.umap(adata_all, color=['Itgam', 'Cd14', 'Cd68'])



sc.tl.rank_genes_groups(adata_all, 'leiden', n_genes=50, method='t-test', use_raw=False)
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(20)
          0              1         2       3        4        5
0    Igfbp5           Cdk8    Igfbp7    Ccn2    Ptprj   Tmsb4x
1      Ptx3  E330009J07Rik      Ccn2   Tagln     Lyz2     Tpt1
2     Dclk1          Cmss1     Aebp1   Acta2     Apoe  Col18a1
3      Gas6         Nkain3      Ctsl   Timp1     Ctsb      Vim
4      Gpx3          Mutyh     Mmp14  Malat1     C1qc     Nme2
5    Tgfbr3        Gm15952       Bgn   Postn   Dock10      Dbi
6      Alpl           Aox2      Sdc2    Tpm2  Gm42418    Actg1
7       Fst           Acpp      Ecm1   Actg2    Hmox1    Postn
8      Ibsp        Gm48742      Psap  Tmsb10     C1qb   Col5a2
9    Plxna2        Slc35f1      Ftl1   Thbs1     C1qa   Eef1a1
10  Adamts1         Cfap53    Itga11    Flna   Tyrobp   Pabpc1
11   Igfbp4        Gm15834    S100a6    Actb   Wfdc17     Cdh2
12     Ccn1          Vash1      Spp1   Cemip     Cd14     Ppia
13       Hp          Cdh19    Rbfox2  Rbfox2   Cyfip1   Igfbp7
14   Olfml3          Egln2  Serpinh1   Cald1    Stab1    Itgb1
15   Chrdl1          Itih4       Lox    Myh9      Lpl   Tuba1a
16     Ndnf          Wdr63    Malat1  Igfbp7     Ccr5    Anxa2
17   Col1a2        Gm28710    Col4a1  Mical2   Fcer1g     Rbm3
18     Comp        Ccdc117     Timp3   Inhba   Frmd4b    Hspa8
19     Mest     Atxn7l1os1    Lgals1  Col1a1     Lgmn   Tmsb10

Keep cluster 0, 1, 2, 3, 5

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
          0       1              2
0      Ptx3    Ccn2          Csmd1
1      Gas6  Igfbp7        Slc35f1
2      Gpx3   Aebp1           Gas7
3      Alpl   Acta2          Uchl1
4       Fst    Ctsl           Cdk8
5      Ibsp    Ecm1          Cdh19
6    Igfbp5  Malat1          Sox10
7    Tgfbr3   Thbs1         Nkain3
8     Dclk1     Bgn           Ngfr
9    Olfml3   Tagln           Sox6
10  Adamts1   Timp1          L1cam
11     Ccn1  Itga11        Fam178b
12     Cst3  Col4a1          Eif5b
13   Plxna2   Timp3           Nebl
14   Col8a1   Loxl3          Moxd1
15   Igfbp4     Fn1        Shroom4
16   Col1a2    Spp1           Klk8
17    Itm2a   Inhba        Cfap299
18   Chrdl1    Sdc2  1700025G04Rik
19       Hp     Lox         Gm2115

#assignment
assignment = pd.read_csv('/Users/raycheng/Dropbox (CahanLab)/BMSC_identity/Data/Sequencing/BMSC_MM_CH_20230127/CH_0127_ARC/tag_assignment_manual_042423.txt', sep = '\t', header = 0)
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
