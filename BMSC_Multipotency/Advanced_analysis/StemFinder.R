# Please download All_day_Signac_all_20230330.RDS from the link below


library(devtools)
#devtools::install_github("pcahan1/stemfinder/stemFinder", ref="main", auth_token = "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")
library(stemFinder)
library(Seurat)
library(loomR)
library(gprofiler2)

adata <- readRDS("All_day_Signac_all_20230330.RDS")

DimPlot(adata, group.by = 'type')

head(adata,2)

DefaultAssay(adata) <- 'RNA'
#adata <- FindVariableFeatures(adata, nfeatures = 6000) 
adata <- FindVariableFeatures(adata, nfeatures = nrow(adata)) 
adata <- ScaleData(adata)
adata <- RunPCA(adata)

ElbowPlot(adata, ndims=50)

#K nearest neighbors
pcs = 20
k = round(sqrt(ncol(adata)))
adata = FindNeighbors(adata, dims = 1:pcs, k.param = k)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Recommended method to prepare murine cc genes: experienced timeout error
convertHumanGeneList <- function(x){
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
humanx <- unique(genesV2[, 2])
# Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}
s_genes_mouse <- convertHumanGeneList(s.genes)
g2m_genes_mouse <- convertHumanGeneList(g2m.genes)

#The alternative way
s_genes_mouse = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
g2m_genes_mouse = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

knn = adata@graphs$RNA_nn
expDat = as.matrix(adata@assays$RNA@scale.data)
cell_cycle_genes = c(s_genes_mouse, g2m_genes_mouse)[c(s_genes_mouse, g2m_genes_mouse) %in% rownames(expDat)]

#Compute single-cell potency
adata = run_stemFinder(adata, k = k, nn = knn, thresh = 0, markers = cell_cycle_genes)
head(adata,2) 

write.csv(adata@meta.data, "stemFinder_All_BMSC_DAY_20230330.csv")

FeaturePlot(adata, features = c('stemFinder_invert'), cols = c('blue','red'))


