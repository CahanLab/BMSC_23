python3
import scvelo as scv
import numpy as np
import pandas as pd
import scanpy as sc
import loompy

# Import Data
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/possorted_genome_bam_3QATU.loom", "possorted_genome_bam_3QATU.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/possorted_genome_bam_FCR1O.loom", "possorted_genome_bam_FCR1O.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/possorted_genome_bam_B18BA.loom", "possorted_genome_bam_B18BA.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/possorted_genome_bam_FPPRE.loom", "possorted_genome_bam_FPPRE.loom")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/regev_lab_cell_cycle_genes.txt", "regev_lab_cell_cycle_genes.txt")

MSC1 = sc.read("possorted_genome_bam_3QATU.loom", cache=True)
MSC2 = sc.read("possorted_genome_bam_FCR1O.loom", cache=True)
MSC3 = sc.read("possorted_genome_bam_B18BA.loom", cache=True)
MSC4 = sc.read("possorted_genome_bam_FPPRE.loom", cache=True)

MSC1.var_names_make_unique()
MSC2.var_names_make_unique()
MSC3.var_names_make_unique()
MSC4.var_names_make_unique()

adata= MSC1.concatenate([MSC2, MSC3, MSC4])
adata.var_names_make_unique()
adata

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1
thresh = np.percentile(adata.obs['n_counts'],95)
adata= adata[adata.obs['n_counts'] < thresh, :]

#Mito_genes filter
# For mice dataset
mito_genes = [name for name in adata.var_names if name.startswith('mt-')]
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X,axis=1)
ribo_genes = [name for name in adata.var_names if name.startswith('Rps') or name.startswith('Rpl')]
adata.obs['percent_ribo'] = np.sum(adata[:,ribo_genes].X,axis=1) / np.sum(adata.X,axis=1)
adata.obs['n_counts'] = adata.X.sum(axis=1)

sc.pl.violin(adata , ['n_genes', 'n_counts', 'percent_ribo'], jitter=0.4, multi_panel=True)


adata= adata[adata.obs['n_genes'] < 3000, :]
adata= adata[adata.obs['percent_ribo'] < 0.4, :]
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

adata_raw = adata.copy()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)


adata= adata[:, adata.var['highly_variable']]

def cycle_score(adata, species):
	cell_cycle_genes = [x.strip() for x in open("regev_lab_cell_cycle_genes.txt")] 
	if species == "mouse":
		cell_cycle_genes = [x.lower() for x in cell_cycle_genes]
		cell_cycle_genes = [x.title() for x in cell_cycle_genes]
		s_genes = cell_cycle_genes[:43]
		g2m_genes = cell_cycle_genes[43:]
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		if cell_cycle_genes == []:
			print("There is no cell cylce genes in the dataset")
		s_genes = [x for x in s_genes if x in adata.var_names]
		g2m_genes = [x for x in g2m_genes if x in adata.var_names]
		sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
		return(adata)
	elif species == "human":
		s_genes = cell_cycle_genes[:43]
		g2m_genes = cell_cycle_genes[43:]
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		if cell_cycle_genes == []:
			print("There is no cell cylce genes in the dataset")
		s_genes = [x for x in s_genes if x in adata.var_names]
		g2m_genes = [x for x in g2m_genes if x in adata.var_names]
		sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
		return(adata)
	else:
		print("please specify the species to be mouse or human")

def oricycle_plot(adata, species):
	cell_cycle_genes = [x.strip() for x in open("regev_lab_cell_cycle_genes.txt")] 
	if species == "mouse":
		cell_cycle_genes = [x.lower() for x in cell_cycle_genes]
		cell_cycle_genes = [x.title() for x in cell_cycle_genes]
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		adata_cc_genes = adata[:, cell_cycle_genes]
		sc.tl.pca(adata_cc_genes,  svd_solver='arpack')
		sc.pl.pca_scatter(adata_cc_genes, color='phase')
	elif species == "human":
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		adata_cc_genes = adata[:, cell_cycle_genes]
		sc.tl.pca(adata_cc_genes,  svd_solver='arpack')
		sc.pl.pca_scatter(adata_cc_genes, color='phase')
	else:
		print("please specify the species to be mouse or human")

def cycle_regress(adata):
	sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
	sc.pp.scale(adata)
	return(adata)

adata = cycle_score(adata, "mouse")
oricycle_plot(adata, "mouse")


adata = cycle_regress(adata)
oricycle_plot(adata, "mouse")


PCA
sc.tl.pca(adata,use_highly_variable=True, n_comps=60, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, 60, log=True)


adata_pre = adata.copy()

adata= adata_pre.copy()
sc.pp.neighbors(adata, n_neighbors=10,n_pcs=50)
sc.tl.louvain(adata, resolution = 0.5)
sc.tl.umap(adata)
adata_all= adata_raw.copy()
adata_all.obsm['X_pca'] = adata.obsm['X_pca']
adata_all.obsm['X_umap'] = adata.obsm['X_umap']
adata_all.obs['louvain'] = adata.obs['louvain']
sc.pl.umap(adata_all, color=['louvain'])



#Finding Marker Genes that consider all the genes (I keep the embedding the same)
sc.tl.rank_genes_groups(adata_all, 'louvain')
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(25)


# Subset non endothelial cells

a= adata_all[adata_all.obs['louvain'].isin(['1', '2','3', '4', '5', '6', '11', '14', '15', '16', '17', '18', '19', '20'])].obs.index
adata= MSC1.concatenate([MSC2, MSC3, MSC4])
adata.var_names_make_unique()
adata

adata = adata[adata.obs.index.isin(a)]
adata

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1
thresh = np.percentile(adata.obs['n_counts'],95)
adata= adata[adata.obs['n_counts'] < thresh, :]

#Mito_genes filter
# For mice dataset
mito_genes = [name for name in adata.var_names if name.startswith('mt-')]
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X,axis=1)
ribo_genes = [name for name in adata.var_names if name.startswith('Rps') or name.startswith('Rpl')]
adata.obs['percent_ribo'] = np.sum(adata[:,ribo_genes].X,axis=1) / np.sum(adata.X,axis=1)
adata.obs['n_counts'] = adata.X.sum(axis=1)

adata= adata[adata.obs['n_genes'] < 3000, :]
adata= adata[adata.obs['percent_ribo'] < 0.4, :]
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

adata_raw = adata.copy()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)


adata= adata[:, adata.var['highly_variable']]


def cycle_score(adata, species):
	cell_cycle_genes = [x.strip() for x in open("regev_lab_cell_cycle_genes.txt")] 
	if species == "mouse":
		cell_cycle_genes = [x.lower() for x in cell_cycle_genes]
		cell_cycle_genes = [x.title() for x in cell_cycle_genes]
		s_genes = cell_cycle_genes[:43]
		g2m_genes = cell_cycle_genes[43:]
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		if cell_cycle_genes == []:
			print("There is no cell cylce genes in the dataset")
		s_genes = [x for x in s_genes if x in adata.var_names]
		g2m_genes = [x for x in g2m_genes if x in adata.var_names]
		sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
		return(adata)
	elif species == "human":
		s_genes = cell_cycle_genes[:43]
		g2m_genes = cell_cycle_genes[43:]
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		if cell_cycle_genes == []:
			print("There is no cell cylce genes in the dataset")
		s_genes = [x for x in s_genes if x in adata.var_names]
		g2m_genes = [x for x in g2m_genes if x in adata.var_names]
		sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
		return(adata)
	else:
		print("please specify the species to be mouse or human")

def oricycle_plot(adata, species):
	cell_cycle_genes = [x.strip() for x in open("regev_lab_cell_cycle_genes.txt")] 
	if species == "mouse":
		cell_cycle_genes = [x.lower() for x in cell_cycle_genes]
		cell_cycle_genes = [x.title() for x in cell_cycle_genes]
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		adata_cc_genes = adata[:, cell_cycle_genes]
		sc.tl.pca(adata_cc_genes,  svd_solver='arpack')
		sc.pl.pca_scatter(adata_cc_genes, color='phase')
	elif species == "human":
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		adata_cc_genes = adata[:, cell_cycle_genes]
		sc.tl.pca(adata_cc_genes,  svd_solver='arpack')
		sc.pl.pca_scatter(adata_cc_genes, color='phase')
	else:
		print("please specify the species to be mouse or human")

def cycle_regress(adata):
	sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
	sc.pp.scale(adata)
	return(adata)

adata = cycle_score(adata, "mouse")
oricycle_plot(adata, "mouse")

adata = cycle_regress(adata)
oricycle_plot(adata, "mouse")

PCA
sc.tl.pca(adata,use_highly_variable=True, n_comps=60, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, 60, log=True)


adata_pre = adata.copy()

adata= adata_pre.copy()
sc.pp.neighbors(adata, n_neighbors=10,n_pcs=50)
sc.tl.louvain(adata, resolution = 0.2)
sc.tl.umap(adata)
adata_all= adata_raw.copy()
adata_all.obsm['X_pca'] = adata.obsm['X_pca']
adata_all.obsm['X_umap'] = adata.obsm['X_umap']
adata_all.obs['louvain'] = adata.obs['louvain']
sc.pl.umap(adata_all, color=['louvain'])

# Subset non endothelial cells

a= adata_all[adata_all.obs['louvain'].isin(['0', '1', '2','3', '5', '6', '8'])].obs.index
adata= MSC1.concatenate([MSC2, MSC3, MSC4])
adata.var_names_make_unique()
adata

adata = adata[adata.obs.index.isin(a)]
adata

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1
thresh = np.percentile(adata.obs['n_counts'],95)
adata= adata[adata.obs['n_counts'] < thresh, :]

#Mito_genes filter
# For mice dataset
mito_genes = [name for name in adata.var_names if name.startswith('mt-')]
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X,axis=1)
ribo_genes = [name for name in adata.var_names if name.startswith('Rps') or name.startswith('Rpl')]
adata.obs['percent_ribo'] = np.sum(adata[:,ribo_genes].X,axis=1) / np.sum(adata.X,axis=1)
adata.obs['n_counts'] = adata.X.sum(axis=1)

adata= adata[adata.obs['n_genes'] < 3000, :]
adata= adata[adata.obs['percent_ribo'] < 0.4, :]
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

adata_raw = adata.copy()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)


adata= adata[:, adata.var['highly_variable']]


def cycle_score(adata, species):
	cell_cycle_genes = [x.strip() for x in open("regev_lab_cell_cycle_genes.txt")] 
	if species == "mouse":
		cell_cycle_genes = [x.lower() for x in cell_cycle_genes]
		cell_cycle_genes = [x.title() for x in cell_cycle_genes]
		s_genes = cell_cycle_genes[:43]
		g2m_genes = cell_cycle_genes[43:]
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		if cell_cycle_genes == []:
			print("There is no cell cylce genes in the dataset")
		s_genes = [x for x in s_genes if x in adata.var_names]
		g2m_genes = [x for x in g2m_genes if x in adata.var_names]
		sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
		return(adata)
	elif species == "human":
		s_genes = cell_cycle_genes[:43]
		g2m_genes = cell_cycle_genes[43:]
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		if cell_cycle_genes == []:
			print("There is no cell cylce genes in the dataset")
		s_genes = [x for x in s_genes if x in adata.var_names]
		g2m_genes = [x for x in g2m_genes if x in adata.var_names]
		sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
		return(adata)
	else:
		print("please specify the species to be mouse or human")

def oricycle_plot(adata, species):
	cell_cycle_genes = [x.strip() for x in open("regev_lab_cell_cycle_genes.txt")] 
	if species == "mouse":
		cell_cycle_genes = [x.lower() for x in cell_cycle_genes]
		cell_cycle_genes = [x.title() for x in cell_cycle_genes]
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		adata_cc_genes = adata[:, cell_cycle_genes]
		sc.tl.pca(adata_cc_genes,  svd_solver='arpack')
		sc.pl.pca_scatter(adata_cc_genes, color='phase')
	elif species == "human":
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		adata_cc_genes = adata[:, cell_cycle_genes]
		sc.tl.pca(adata_cc_genes,  svd_solver='arpack')
		sc.pl.pca_scatter(adata_cc_genes, color='phase')
	else:
		print("please specify the species to be mouse or human")

def cycle_regress(adata):
	sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
	sc.pp.scale(adata)
	return(adata)

adata = cycle_score(adata, "mouse")
oricycle_plot(adata, "mouse")

adata = cycle_regress(adata)
oricycle_plot(adata, "mouse")


PCA
sc.tl.pca(adata,use_highly_variable=True, n_comps=60, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, 60, log=True)

adata_pre = adata.copy()

adata= adata_pre.copy()
sc.pp.neighbors(adata, n_neighbors=10,n_pcs=50)
sc.tl.louvain(adata, resolution = 0.2)
sc.tl.umap(adata)
adata_all= adata_raw.copy()
adata_all.obsm['X_pca'] = adata.obsm['X_pca']
adata_all.obsm['X_umap'] = adata.obsm['X_umap']
adata_all.obs['louvain'] = adata.obs['louvain']
sc.pl.umap(adata_all, color=['louvain'])


#Finding Marker Genes that consider all the genes (I keep the embedding the same)
sc.tl.rank_genes_groups(adata_all, 'louvain')
pd.DataFrame(adata_all.uns['rank_genes_groups']['names']).head(25)
