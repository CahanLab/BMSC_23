# Download the file before start
"s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/All_BMSC_marrow_clustered_bc_20230329.loom"

# Prepare environment
library(devtools)
library(stemFinder)
library(Seurat)
library(loomR)
library(gprofiler2)
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SeuratDisk)
library(loomR)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(GenomicRanges)
library(future)
library(biovizBase)
library(Pando)
library(grr)
library(chromVARmotifs)
library(Pando)
library(magrittr)

# Preapre data in R
lfile <- connect(filename = "All_BMSC_marrow_clustered_bc_20230329.loom", skip.validate = TRUE)
geneNames<-lfile[["row_attrs"]][["var_names"]][]
cellNames<-lfile[["col_attrs"]][["obs_names"]][]
expPark<- t(lfile[["matrix"]][1:length(cellNames),])
louvain<-lfile[["col_attrs"]][["leiden0"]][]
cell_type<-lfile[["col_attrs"]][["type1"]][]
rownames(expPark)<-geneNames
colnames(expPark)<-cellNames
stPark<- data.frame(row.names = cellNames, cell_name=cellNames, cluster=louvain, type = cell_type)
genesPark<-rownames(expPark)

#row.names(ad$obs) <- ad$obs$new_name
#row.names(ad$var) <- ad$var$var_names
adata_Seurat <- CreateSeuratObject(counts = expPark, meta.data = stPark)

# Remove numeric genes: For GRN
test <- rownames(adata_Seurat)
test <- subset(test, !startsWith(test, '0'))
test <- subset(test, !startsWith(test, '1'))
test <- subset(test, !startsWith(test, '2'))
test <- subset(test, !startsWith(test, '3'))
test <- subset(test, !startsWith(test, '4'))
test <- subset(test, !startsWith(test, '5'))
test <- subset(test, !startsWith(test, '6'))
test <- subset(test, !startsWith(test, '7'))
test <- subset(test, !startsWith(test, '8'))
test <- subset(test, !startsWith(test, '9'))

adata <- subset(adata_Seurat, features = test)
adata<- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(adata)
adata<- ScaleData(adata, features = all.genes)

adata<- RunPCA(adata, features = VariableFeatures(object = adata))
adata<- RunUMAP(adata, dims = 1:10)


# Assign coordinates:
UMAP <- lfile[["col_attrs"]][["X_umap"]][,1:length(cellNames)]
adata[["umap"]]@cell.embeddings <- t(UMAP)
colnames(adata[["umap"]]@cell.embeddings) <- c('UMAP_1', 'UMAP_2')
rownames(adata[["umap"]]@cell.embeddings) <- colnames(adata)

DimPlot(adata, group.by = 'type')



head(adata,2)

DefaultAssay(adata) <- 'RNA'
#adata <- FindVariableFeatures(adata, nfeatures = 6000) 
adata <- FindVariableFeatures(adata, nfeatures = nrow(adata)) 
adata <- ScaleData(adata)
adata <- RunPCA(adata)

ElbowPlot(adata, ndims=50)



#K nearest neighbors
pcs = 30
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

write.csv(adata@meta.data, "stemFinder_All_BMSC_20230330.csv")





