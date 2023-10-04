# RNA velocity
# Stochastic model
import scvelo as scv
import scanpy as sc
adata_bone = adata_all[adata_all.obs['tissue'] == 'bone']
adata_marrow = adata_all[adata_all.obs['tissue'] == 'marrow']

# Bone fraction
scv.tl.velocity(adata_bone, mode='stochastic')
scv.tl.velocity_graph(adata_bone)
scv.tl.velocity_embedding(adata_bone, basis='umap')

scv.pl.velocity_embedding_stream(adata_bone, basis='umap', color = 'type', title = '', legend_loc = 'right')


# Marrow fraction
scv.tl.velocity(adata_marrow, mode='stochastic')
scv.tl.velocity_graph(adata_marrow)
scv.tl.velocity_embedding(adata_marrow, basis='umap')

scv.pl.velocity_embedding_stream(adata_marrow, basis='umap', color = 'type', title = '', legend_loc = 'right')


