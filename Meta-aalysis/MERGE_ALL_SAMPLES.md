
# merge all the individual sample together
```python
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

# Identify subclusters
adata_marrow = adata_all[adata_all.obs['tissue'] == 'marrow']
adata_bone = adata_all[adata_all.obs['tissue'] == 'bone']

sc.pp.neighbors(adata_marrow , use_rep="X_pca_harmony", n_neighbors=25,n_pcs=25)
sc.tl.leiden(adata_marrow , resolution = 0.1, restrict_to=('type', ['BMSC2']), key_added='leiden1')
sc.pl.umap(adata_marrow , color = ['leiden1'], palette = 'Set2')

sc.tl.rank_genes_groups(adata_marrow, 'leiden1', n_genes=adata_marrow.shape[1], use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata_marrow, values_to_plot='logfoldchanges', min_logfoldchange=1.5, vmax=7, vmin=-7, cmap='bwr', dendrogram = False)
```
