# RNA velocity
# Stochastic model
import scvelo as scv
import scanpy as sc

download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/BMSC_AD_OS_CH_clustered_bc_20231204.loom", "BMSC_AD_OS_CH_clustered_bc_20230425.loom")
adata_all = sc.read("BMSC_AD_OS_CH_clustered_bc_20230425.loom")
adata_all.obs_names = adata_all.obs['obs_names']
adata_all.var_names = adata_all.var['var_names']

adata_AD= adata_all[adata_all.obs['lineage'] == 'AD']
adata_OS= adata_all[adata_all.obs['lineage'] == 'OS']
adata_CH= adata_all[adata_all.obs['lineage'] == 'CH']

# Reading in unspliced and spliced counts
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/gex_possorted_bam_8ZO19.loom", "gex_possorted_bam_8ZO19.loom")
adata1 = scv.read("gex_possorted_bam_8ZO19.loom")
adata1.obs.index = adata1.obs.index.str[:-1]
adata1.obs.index = adata1.obs.index.str[24:]
adata1.obs.index = adata1.obs.index + "-1"
adata1.var_names_make_unique()

download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/gex_possorted_bam_AWZ0K.loom", "gex_possorted_bam_AWZ0K.loom")
adata2 = scv.read("gex_possorted_bam_AWZ0K.loom")
adata2.obs.index = adata2.obs.index.str[:-1]
adata2.obs.index = adata2.obs.index.str[24:]
adata2.obs.index = adata2.obs.index + "-1"
adata2.var_names_make_unique()

download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/gex_possorted_bam_316SD.loom", "gex_possorted_bam_316SD.loom")
adata3 = scv.read("gex_possorted_bam_316SD.loom")
adata3.obs.index = adata3.obs.index.str[:-1]
adata3.obs.index = adata3.obs.index.str[24:]
adata3.obs.index = adata3.obs.index + "-1"
adata3.var_names_make_unique()

adata = adata1.concatenate([adata2, adata3])
adata

adVary = adata_all.copy()
adM = adata.copy()
adM = adM[  adVary.obs_names,:]
adM = adM[:,  adVary.var_names]
adM.obsm = adVary.obsm.copy()
adM.obs = adVary.obs.copy()

# Stoch model
scv.tl.velocity(adM , mode='stochastic')
scv.tl.velocity_graph(adM)

scv.tl.velocity_embedding(adM , basis='umap')
scv.pl.velocity_embedding(adM , basis='umap', alpha=.25, arrow_size =1.5, arrow_length=2)
adM_stoch = adM.copy()
scv.pl.velocity_embedding_stream(adM_stoch, basis='umap', color = 'type', title = '', legend_loc = 'right', palette = ['#E97E73', '#CA9734', '#9EAF3C', '#7DC874', '#73C9B0', '#64BEE1', '#7FA5F9', '#CE78F4', '#F292D0'])

scv.pl.velocity_embedding_stream(adM , basis='umap', color = 'multitag', title = '', legend_loc = 'right')

scv.tl.recover_dynamics(adM, var_names='velocity_genes')
scv.tl.latent_time(adM)
scv.pl.scatter(adM, color='latent_time', color_map='gnuplot', size=80, colorbar=True)

# Dynamic model
####
scv.pp.moments(adM, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adM, var_names='velocity_genes')
scv.tl.velocity(adM, mode='dynamical')
scv.tl.velocity_graph(adM)

scv.tl.velocity_embedding(adM, basis='umap')
# scv.pl.velocity_embedding(adM, basis='umap', alpha=.25, arrow_size =2.5, arrow_length=2)
scv.pl.velocity_embedding_stream(adM, basis='umap', title = '', legend_loc = 'right', color = 'type', palette = ['#E97E73', '#CA9734', '#9EAF3C', '#7DC874', '#73C9B0', '#64BEE1', '#7FA5F9', '#CE78F4', '#F292D0'])

