# RNA velocity
# Stochastic model
import scvelo as scv
import scanpy as sc
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/BMSC_All_day_clustered_bc_041923.loom", "BMSC_All_day_clustered_bc_041923.loom")
adata_all = sc.read("BMSC_All_day_clustered_bc_041923.loom")

adata_all.obs_names = adata_all.obs['obs_names']
adata_all.var_names = adata_all.var['var_names']

# Reading in unspliced and spliced counts
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/d0_gex_possorted_bam_NH4B3.loom", "d0_gex_possorted_bam_NH4B3.loom")
adata1 = scv.read("d0_gex_possorted_bam_NH4B3.loom")
adata1.obs.index = adata1.obs.index.str[:-1]
adata1.obs.index = adata1.obs.index.str[24:]
adata1.obs.index = adata1.obs.index + "-1"
adata1.var_names_make_unique()

download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/gex_possorted_bam_E8J1G.loom", "gex_possorted_bam_E8J1G.loom")
adata4 = scv.read("gex_possorted_bam_E8J1G.loom")
adata4.obs.index = adata4.obs.index.str[:-1]
adata4.obs.index = adata4.obs.index.str[24:]
adata4.obs.index = adata4.obs.index + "-1"
adata4.var_names_make_unique()

download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/d1_gex_possorted_bam_A5CIC.loom", "d1_gex_possorted_bam_A5CIC.loom")
adata2 = scv.read("d1_gex_possorted_bam_A5CIC.loom")
adata2.obs.index = adata2.obs.index.str[:-1]
adata2.obs.index = adata2.obs.index.str[24:]
adata2.obs.index = adata2.obs.index + "-1"
adata2.var_names_make_unique()

download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/d3_gex_possorted_bam_GL9C5.loom", "d3_gex_possorted_bam_GL9C5.loom")
adata3 = scv.read("d3_gex_possorted_bam_GL9C5.loom")
adata3.obs.index = adata3.obs.index.str[:-1]
adata3.obs.index = adata3.obs.index.str[24:]
adata3.obs.index = adata3.obs.index + "-1"
adata3.var_names_make_unique()

adata = adata1.concatenate([adata2, adata3, adata4])
adata

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

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
# scv.pl.velocity_embedding(adata_day1 , basis='umap', alpha=.25, arrow_size =1.5, arrow_length=2)
scv.pl.velocity_embedding_stream(adM , basis='umap', color = 'type3', title = '', legend_loc = 'right')

# Dynamic model
####
scv.pp.moments(adM, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adM, var_names='velocity_genes')
scv.tl.velocity(adM, mode='dynamical')
scv.tl.velocity_graph(adM)

scv.tl.velocity_embedding(adM, basis='umap')
# scv.pl.velocity_embedding(adM, basis='umap', alpha=.25, arrow_size =2.5, arrow_length=2)
scv.pl.velocity_embedding_stream(adM, basis='umap', title = '', legend_loc = 'right', color = 'type3')


