python3
import scvelo as scv
import numpy as np
import pandas as pd
import scanpy as sc
import loompy

# Import Data
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/Addo_ZYBAF.loom", "Addo_ZYBAF.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/GSM3466901_YRHA5.loom", "GSM3466901_YRHA5.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/GSM3466900_BE8CC.loom", "GSM3466900_BE8CC.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/GSM3937217_VRRW0.loom", "GSM3937217_VRRW0.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/GSM3937216_D6GM2.loom", "GSM3937216_D6GM2.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/bm1_5VFON.loom", "bm1_5VFON.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/bm2_QMWB7.loom", "bm2_QMWB7.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/bm3_XT1AR.loom", "bm3_XT1AR.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/bm4_R4W09.loom", "bm4_R4W09.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/regev_lab_cell_cycle_genes.txt", "regev_lab_cell_cycle_genes.txt")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/b1_4E06Z.loom", "b1_4E06Z.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/b2_FGP7Y.loom", "b2_FGP7Y.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/b3_FJEW2.loom", "b3_FJEW2.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/b4_QJS0R.loom", "b4_QJS0R.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/regev_lab_cell_cycle_genes.txt", "regev_lab_cell_cycle_genes.txt")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/possorted_genome_bam_3QATU.loom", "possorted_genome_bam_3QATU.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/possorted_genome_bam_FCR1O.loom", "possorted_genome_bam_FCR1O.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/possorted_genome_bam_B18BA.loom", "possorted_genome_bam_B18BA.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/possorted_genome_bam_FPPRE.loom", "possorted_genome_bam_FPPRE.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/Tikhonova_MSC_12wk_SKB70.loom", "Tikhonova_MSC_12wk_SKB70.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/Tikhonova_MSC_7wk_0BK8W.loom", "Tikhonova_MSC_7wk_0BK8W.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/Wolock_MSC_QCTZ0.loom", "Wolock_MSC_QCTZ0.loom")

# Prepare SCN classification gene list for embedding calculation
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/SCN_genes_marrow_20210817.csv", "SCN_genes_marrow_20210817.csv")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/SCN_genes.csv", "SCN_genes.csv")
genes = pd.read_csv("SCN_genes_marrow_20210817.csv", header = 0)
genes2 = pd.read_csv("SCN_genes.csv", header = 0)
genes_final = pd.merge(genes['x'], genes2['x'], how='outer')
genes_final = genes_final.drop_duplicates()

# Baryawno_BMSC-MSC
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/Baryawno_MSC_OsP_clustered.loom", "Baryawno_MSC_OsP_clustered.loom")
adata_all1 = sc.read("Baryawno_MSC_OsP_clustered.loom")

BMSC_1 = sc.read("bm1_5VFON.loom", cache=True)
BMSC_2 = sc.read("bm2_QMWB7.loom", cache=True)
BMSC_3 = sc.read("bm3_XT1AR.loom", cache=True)
BMSC_4 = sc.read("bm4_R4W09.loom", cache=True)

BMSC_1.var_names_make_unique()
BMSC_2.var_names_make_unique()
BMSC_3.var_names_make_unique()
BMSC_4.var_names_make_unique()

adata1= BMSC_1.concatenate([BMSC_2, BMSC_3, BMSC_4])
a= adata_all1.obs.index
adata1= adata1[adata1.obs.index.isin(a)]
adata1.obs = adata_all1.obs

adata1.obs['Cell Type'] = ''
conditions = [(adata1.obs['louvain'] == '0'), (adata1.obs['louvain'] == '1'), (adata1.obs['louvain'] == '2'), (adata1.obs['louvain'] == '3'), (adata1.obs['louvain'] == '4'), (adata1.obs['louvain'] == '5'), (adata1.obs['louvain'] == '6')]
choices = ['CAR cell', 'Preadipocyte', 'Osteoprogenitor', 'Osteoprogenitor', 'Schwann cell ', 'Fibroblast', 'IFN responsive cell']
adata1.obs['Cell Type'] = np.select(conditions, choices, default='0')
adata1= adata1[:, genes_final['x']]
adata1.obs['study'] = 'Baryawno et al'
adata1.obs['tissue'] = 'marrow'

sc.pp.normalize_per_cell(adata1, counts_per_cell_after=1e4)
sc.pp.log1p(adata1)

# Baccin2
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/Baccin_BMSC4_extracted_clustered.loom", "Baccin_BMSC4_extracted_clustered.loom")
adata3=  sc.read("GSM3937217_VRRW0.loom", cache=True)
adata_all3 = sc.read("Baccin_BMSC4_extracted_clustered.loom", cache=True)

c= adata_all3.obs.index
adata3 = adata3[adata3.obs.index.isin(c)]
adata3.obs = adata_all3.obs

#sc.pl.umap(adata_all3, color = 'louvain')
adata3.obs['Cell Type'] = ''
conditions = [(adata3.obs['louvain'] == '0'), (adata3.obs['louvain'] == '1'), (adata3.obs['louvain'] == '2'), (adata3.obs['louvain'] == '3'), (adata3.obs['louvain'] == '4'), (adata3.obs['louvain'] == '5'), (adata3.obs['louvain'] == '6')]
choices = ['CAR cell', 'Preadipocyte', 'Osteoprogenitor', 'Osteoprogenitor', 'IFN responsive cell', 'Smooth muscle cell', 'Fibroblast']
adata3.obs['Cell Type'] = np.select(conditions, choices, default='0')
adata3.var_names_make_unique()
adata3= adata3[:, genes_final['x']]
adata3.obs['study'] = 'Baccin et al'
adata3.obs['tissue'] = 'marrow'

sc.pp.normalize_per_cell(adata3, counts_per_cell_after=1e4)
sc.pp.log1p(adata3)

#Addo_MSC
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/Addo_MSC_clustered.loom", "Addo_MSC_clustered.loom")
adata4 = scv.read("Addo_ZYBAF.loom", cache=True)
adata_all4 = sc.read("Addo_MSC_clustered.loom", cache=True)
d= adata_all4.obs.index
adata4= adata4[adata4.obs.index.isin(d)]
adata4.obs = adata_all4.obs

#sc.pl.umap(adata_all4, color = 'louvain')
adata4.obs['Cell Type'] = ''
conditions = [(adata4.obs['louvain'] == '0'), (adata4.obs['louvain'] == '1'), (adata4.obs['louvain'] == '2')]
choices = ['CAR cell', 'Osteoprogenitor', 'Preadipocyte']
adata4.obs['Cell Type'] = np.select(conditions, choices, default='0')
adata4.var_names_make_unique()
adata4= adata4[:, genes_final['x']]
adata4.obs['study'] = 'Addo et al'
adata4.obs['tissue'] = 'marrow'

sc.pp.normalize_per_cell(adata4, counts_per_cell_after=1e4)
sc.pp.log1p(adata4)

#Tikhonova_7wk
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/Tikhonova_MSC_7wk_clustered.loom", "Tikhonova_MSC_7wk_clustered.loom")
adata5= sc.read("Tikhonova_MSC_7wk_0BK8W.loom", cache=True)
adata_all5 = sc.read("Tikhonova_MSC_7wk_clustered.loom", cache=True)

e= adata_all5.obs.index
adata5= adata5[adata5.obs.index.isin(e)]
adata5.obs = adata_all5.obs

#sc.pl.umap(adata_all5, color = 'louvain')
adata5.obs['Cell Type'] = ''
conditions = [(adata5.obs['louvain'] == '0'), (adata5.obs['louvain'] == '1'), (adata5.obs['louvain'] == '2')]
choices = ['CAR cell', 'Preadipocyte', 'Osteoprogenitor']
adata5.obs['Cell Type'] = np.select(conditions, choices, default='0')
adata5.var_names_make_unique()
adata5= adata5[:, genes_final['x']]
adata5.obs['study'] = 'Tikhonova et al'
adata5.obs['tissue'] = 'marrow'

sc.pp.normalize_per_cell(adata5, counts_per_cell_after=1e4)
sc.pp.log1p(adata5)

#Tikhonova_12wk
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/Tikhonova_MSC_12wk_clustered.loom", "Tikhonova_MSC_12wk_clustered.loom")
adata6= sc.read("Tikhonova_MSC_12wk_SKB70.loom", cache=True)
adata_all6 = sc.read("Tikhonova_MSC_12wk_clustered.loom", cache=True)

f= adata_all6.obs.index
adata6= adata6[adata6.obs.index.isin(f)]
adata6.obs = adata_all6.obs

#sc.pl.umap(adata_all6, color = 'louvain')
adata6.obs['Cell Type'] = ''
conditions = [(adata6.obs['louvain'] == '0'), (adata6.obs['louvain'] == '1'), (adata6.obs['louvain'] == '2')]
choices = ['CAR cell', 'Osteoprogenitor', 'Preadipocyte']
adata6.obs['Cell Type'] = np.select(conditions, choices, default='0')
adata6.var_names_make_unique()
adata6= adata6[:, genes_final['x']]
adata6.obs['study'] = 'Tikhonova et al'
adata6.obs['tissue'] = 'marrow'

sc.pp.normalize_per_cell(adata6, counts_per_cell_after=1e4)
sc.pp.log1p(adata6)

# Baryawno_BMC-MSC
B_1 = sc.read("Baryawno_Cell_2019/b1_4E06Z.loom", cache=True)
B_2 = sc.read("Baryawno_Cell_2019/b2_FGP7Y.loom", cache=True)
B_3 = sc.read("Baryawno_Cell_2019/b3_FJEW2.loom", cache=True)
B_4 = sc.read("Baryawno_Cell_2019/b4_QJS0R.loom", cache=True)

download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/Baryawno_BMC_MSC_OsP_122719_clustered.loom", "Baryawno_BMC_MSC_OsP_122719_clustered.loom")
adata_all7 = sc.read("Baryawno_BMC_MSC_OsP_122719_clustered.loom", cache=True)
B_1.var_names_make_unique()
B_2.var_names_make_unique()
B_3.var_names_make_unique()
B_4.var_names_make_unique()

adata7= B_1.concatenate([B_2, B_3, B_4])
g= adata_all7.obs.index
adata7= adata7[adata7.obs.index.isin(g)]
adata7.obs = adata_all7.obs

#sc.pl.umap(adata_all7, color = 'louvain')

adata7.obs['Cell Type'] = ''
conditions = [(adata7.obs['louvain'] == '0'), (adata7.obs['louvain'] == '1'), (adata7.obs['louvain'] == '2'), (adata7.obs['louvain'] == '3'), (adata7.obs['louvain'] == '4'), (adata7.obs['louvain'] == '5'), (adata7.obs['louvain'] == '6'), (adata7.obs['louvain'] == '7')]
choices = ['Fibrochondroprogenitor', 'Chondrocyte', 'Smooth muscle cell', 'Fibroblast', 'CAR cell', 'Tenocyte', 'Osteoprogenitor', 'Synoviocyte']
adata7.obs['Cell Type'] = np.select(conditions, choices, default='0')
adata7= adata7[:, genes_final['x']]
adata7.obs['study'] = 'Baryawno et al'
adata7.obs['tissue'] = 'bone'

sc.pp.normalize_per_cell(adata7, counts_per_cell_after=1e4)
sc.pp.log1p(adata7)

# Baccin
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/Baccin_BMSC2_extracted_clustered.loom", "Baccin_BMSC2_extracted_clustered.loom")
adata8=  sc.read("GSM3466901_YRHA5.loom", cache=True)
adata_all8 = sc.read("Baccin_BMSC2_extracted_clustered.loom", cache=True)

h= adata_all8.obs.index
adata8= adata8[adata8.obs.index.isin(h)]
adata8.obs = adata_all8.obs

#sc.pl.umap(adata_all8, color = 'louvain')

adata8.obs['Cell Type'] = ''
conditions = [(adata8.obs['louvain'] == '0'), (adata8.obs['louvain'] == '1'), (adata8.obs['louvain'] == '2'), (adata8.obs['louvain'] == '3'), (adata8.obs['louvain'] == '4'), (adata8.obs['louvain'] == '5'), (adata8.obs['louvain'] == '6'), (adata8.obs['louvain'] == '7'), (adata8.obs['louvain'] == '8'), (adata8.obs['louvain'] == '9')]
choices = ['Fibroblast', 'Fibroblast', 'CAR cell', 'Fibroblast', 'Fibroblast', 'Fibroblast', 'Chondrocyte', 'Smooth muscle cell', 'Osteoprogenitor', 'Synoviocyte']
adata8.obs['Cell Type'] = np.select(conditions, choices, default='0')
adata8.var_names_make_unique()
adata8= adata8[:, genes_final['x']]
adata8.obs['study'] = 'Baccin et al'
adata8.obs['tissue'] = 'bone'

sc.pp.normalize_per_cell(adata8, counts_per_cell_after=1e4)
sc.pp.log1p(adata8)


# Sivaraj
MSC1 = sc.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Sivaraj_BMSC_Cell_Reports/GSE156636/GSM4735393/outs/velocyto/possorted_genome_bam_3QATU.loom", cache=True)
MSC2 = sc.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Sivaraj_BMSC_Cell_Reports/GSE156636/GSM4735394/outs/velocyto/possorted_genome_bam_FCR1O.loom", cache=True)
MSC3 = sc.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Sivaraj_BMSC_Cell_Reports/GSE156636/GSM4735395/outs/velocyto/possorted_genome_bam_B18BA.loom", cache=True)
MSC4 = sc.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Sivaraj_BMSC_Cell_Reports/GSE156636/GSM4735396/outs/velocyto/possorted_genome_bam_FPPRE.loom", cache=True)

MSC1.var_names_make_unique()
MSC2.var_names_make_unique()
MSC3.var_names_make_unique()
MSC4.var_names_make_unique()

adata9= MSC1.concatenate([MSC2, MSC3, MSC4])
adata9.var_names_make_unique()

adata_all9 = sc.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Sivaraj_BMSC_Cell_Reports/Sivaraj_BMSC_clustered_20210818-3.loom", cache=True)

i= adata_all9.obs.index
adata9= adata9[adata9.obs.index.isin(i)]
adata9.obs = adata_all9.obs

#sc.pl.umap(adata_all9, color = 'louvain')

adata9.obs['Cell Type'] = ''
conditions = [(adata9.obs['louvain'] == '0'), (adata9.obs['louvain'] == '1'), (adata9.obs['louvain'] == '2'), (adata9.obs['louvain'] == '3'), (adata9.obs['louvain'] == '4'), (adata9.obs['louvain'] == '5'), (adata9.obs['louvain'] == '6')]
choices = ['CAR cell', 'Chondrocyte', 'Chondrocyte', 'Osteoprogenitor',  'Osteoprogenitor', 'Smooth muscle cell', 'Fibroblast']
adata9.obs['Cell Type'] = np.select(conditions, choices, default='0')
adata9.var_names_make_unique()
adata9= adata9[:, genes_final['x']]
adata9.obs['study'] = 'Sivaraj et al'
adata9.obs['tissue'] = 'bone'

sc.pp.normalize_per_cell(adata9, counts_per_cell_after=1e4)
sc.pp.log1p(adata9)

# Wolock

adata10=  sc.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Wolock_CR_2019/Alignment/Wolock_final/Dropest/velocyto/Wolock_MSC_QCTZ0.loom", cache=True)

adata_all10 = sc.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Wolock_CR_2019/Wolock_BMSC_MSC_clustered.loom", cache=True)

j= adata_all10.obs.index
adata10= adata10[adata10.obs.index.isin(j)]
adata10.obs = adata_all10.obs

#sc.pl.umap(adata_all10, color = 'louvain')

adata10.obs['Cell Type'] = ''
conditions = [(adata10.obs['louvain'] == '0'), (adata10.obs['louvain'] == '1'), (adata10.obs['louvain'] == '2'), (adata10.obs['louvain'] == '3'), (adata10.obs['louvain'] == '4'), (adata10.obs['louvain'] == '5'), (adata10.obs['louvain'] == '6')]
choices = ['CAR cell', 'Preadipocyte', 'Osteoprogenitor', 'Osteoprogenitor',  'Osteoprogenitor', 'Smooth muscle cell', 'Osteoprogenitor']
adata10.obs['Cell Type'] = np.select(conditions, choices, default='0')
adata10.var_names_make_unique()
adata10= adata10[:, genes_final['x']]
adata10.obs['study'] = 'Wolock et al'
adata10.obs['tissue'] = 'marrow'

sc.pp.normalize_per_cell(adata10, counts_per_cell_after=1e4)
sc.pp.log1p(adata10)

# Compute embedding
adata1.var_names_make_unique()
adata3.var_names_make_unique()
adata4.var_names_make_unique()
adata5.var_names_make_unique()
adata6.var_names_make_unique()
adata7.var_names_make_unique()
adata8.var_names_make_unique()
adata9.var_names_make_unique()
adata10.var_names_make_unique()
adata = adata1.concatenate([adata3, adata4, adata5, adata6, adata7, adata8, adata9, adata10])

sc.tl.pca(adata,use_highly_variable=False, n_comps=60, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, 60, log=True)



adata_pre = adata.copy()
adata= adata_pre.copy()
sc.pp.neighbors(adata, n_neighbors=300,n_pcs=50)
sc.tl.umap(adata)
sc.pl.umap(adata, color = 'Cell Type')
sc.pl.umap(adata, color = 'study')
sc.pl.umap(adata, color = 'tissue')



#Overlap all the genes
adata1_raw= BMSC_1.concatenate([BMSC_2, BMSC_3, BMSC_4])
adata3_raw=  sc.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Baccin_NCB_2020/GSM3937217/outs/velocyto/GSM3937217_VRRW0.loom", cache=True)
adata4_raw = scv.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Addo_EJI_2019/Addo_ZYBAF.loom", cache=True)
adata5_raw= sc.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Tikhonova_Nature_2019/Tikhonova_MSC_7wk_0BK8W.loom", cache=True)
adata6_raw= sc.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Tikhonova_Nature_2019/Tikhonova_MSC_12wk_SKB70.loom", cache=True)
adata7_raw= B_1.concatenate([B_2, B_3, B_4])
adata8_raw=  sc.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Baccin_NCB_2020/GSM3466901/GSM3466901_YRHA5.loom", cache=True)
adata9_raw= MSC1.concatenate([MSC2, MSC3, MSC4])
adata10_raw= sc.read("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Wolock_CR_2019/Alignment/Wolock_final/Dropest/velocyto/Wolock_MSC_QCTZ0.loom", cache=True)

adata1_raw.var_names_make_unique()
adata3_raw.var_names_make_unique()
adata4_raw.var_names_make_unique()
adata5_raw.var_names_make_unique()
adata6_raw.var_names_make_unique()
adata7_raw.var_names_make_unique()
adata8_raw.var_names_make_unique()
adata9_raw.var_names_make_unique()
adata10_raw.var_names_make_unique()

adata_raw = adata1_raw.concatenate([adata3_raw, adata4_raw, adata5_raw, adata6_raw, adata7_raw, adata8_raw, adata9_raw, adata10_raw])
adata_raw = adata_raw[adata.obs_names]
adata_raw.obs = adata.obs
adata_raw.obsm['X_pca'] = adata.obsm['X_pca']
adata_raw.obsm['X_umap'] = adata.obsm['X_umap']

sc.pp.normalize_per_cell(adata_raw, counts_per_cell_after=1e4)
sc.pp.log1p(adata_raw)

sc.pl.umap(adata_all, color = ['Cxcl12', 'Fos', 'Bglap', 'Ly6a' ,'Col2a1', 'Anxa8', 'Prg4', 'Mbp', 'Myh11', 'Tnmd', 'Cxcl10'])

adata_raw.write_loom('All_BMSC_clustered_20210928.loom', write_obsm_varm=True)

def cluster_small_multiples(adata, clust_key, size=60, frameon=False, legend_loc=None, **kwargs):
    tmp = adata.copy()
    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]
    sc.pl.umap(tmp, groups=tmp.obs[clust].cat.categories[1:].values, color=adata.obs[clust_key].cat.categories.tolist(), size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)

test = adata_all.copy()
cluster_small_multiples(test, 'Cell Type')

# Identify clusters
sc.pp.neighbors(test, n_neighbors=100,n_pcs=25)
sc.tl.leiden(test, resolution = 0.1)
sc.pl.umap(adata_all, color = ['leiden'])


# Batch corrected embeding
import scanpy.external as sce
sce.pp.harmony_integrate(adata_all, 'batch')
sc.pp.neighbors(adata_all, use_rep="X_pca_harmony", n_neighbors=300,n_pcs=40)
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color = ['leiden5', 'type1', 'tissue'])


sc.tl.leiden(adata_all, key_added="clusters", resolution = 0.5)
sc.pl.umap(adata_all, color = ['clusters', 'tissue', 'study'])

BMSC1 CAR
BMSC2 OSP
BMSC3 CH
BMSC4 FIB
BMSC5 SMC

# Identify subclusters
adata_marrow = adata_all[adata_all.obs['tissue'] == 'marrow']
adata_bone = adata_all[adata_all.obs['tissue'] == 'bone']

sc.pp.neighbors(adata_marrow , use_rep="X_pca_harmony", n_neighbors=25,n_pcs=25)
sc.tl.leiden(adata_marrow , resolution = 0.1, restrict_to=('type', ['BMSC2']), key_added='leiden1')
sc.pl.umap(adata_marrow , color = ['leiden1'], palette = 'Set2')

sc.tl.rank_genes_groups(adata_marrow, 'leiden1', n_genes=adata_marrow.shape[1], use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_marrow, values_to_plot='logfoldchanges', min_logfoldchange=1.5, vmax=7, vmin=-7, cmap='bwr', dendrogram = False)


# Rename the clusters
adata_all.obs['type'] = ''
conditions = [(adata_all.obs['clusters'] == '0'), (adata_all.obs['clusters'] == '3'), (adata_all.obs['clusters'] == '2'), (adata_all.obs['clusters'] == '1'), (adata_all.obs['clusters'] == '4')]
choices = ['BMSC1', 'BMSC2', 'BMSC3', 'BMSC4', 'BMSC5']
adata_all.obs['type'] = np.select(conditions, choices, default='0')

adata_marrow.obs['type1'] = ''
conditions = [(adata_marrow.obs['leiden1'] == 'BMSC1'), (adata_marrow.obs['leiden1'] == 'BMSC2,0'), (adata_marrow.obs['leiden1'] == 'BMSC2,1'), (adata_marrow.obs['leiden1'] == 'BMSC2,2'), (adata_marrow.obs['leiden1'] == 'BMSC3'), (adata_marrow.obs['leiden1'] == 'BMSC4'), (adata_marrow.obs['leiden1'] == 'BMSC5')]
choices = ['BMSC1', 'BMSC2', 'BMSC3', 'BMSC4', 'BMSC5', 'BMSC6', 'BMSC7']
adata_marrow.obs['type1'] = np.select(conditions, choices, default='0')


# differential gene heatmap
non_mito_genes_list = [name for name in adata_all.var_names if not name.startswith('mt-')]
adata_all = adata_all[:, non_mito_genes_list]
non_ribo_genes_list = [name for name in adata_all.var_names if not name.startswith('Rps')]
adata_all = adata_all[:, non_ribo_genes_list]
non_ribo_genes_list2 = [name for name in adata_all.var_names if not name.startswith('Rpl')]
adata_all = adata_all[:, non_ribo_genes_list2]

sc.tl.rank_genes_groups(adata_all, 'type', n_genes=200, use_raw=False)
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(20)


