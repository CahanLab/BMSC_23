
# Merge multiple ATAC with common feature sets
R
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SeuratDisk)
library(loomR)
library(patchwork)
library(GenomicRanges)
library(future)
library(biovizBase)
# For GRN
library(Pando)
library(grr)
#library(BSgenome.Hsapiens.UCSC.hg38) # human
# For motif analysis
library(chromVARmotifs)
library(Pando)
library(magrittr)
library(JASPAR2022)
library(TFBSTools)
library(motifmatchr)
library(chromVARmotifs)
library(magrittr)
# For peak linking
library(qlcMatrix)
library(ggforce)
# For plotting
library(ggplot2)
library(ggrepel)
library(plotly)

#plan("multiprocess", workers = 4)
#options(future.globals.maxSize = 16000 * 1024^2)

# Prepare bed files
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day0_0824_ARC/atac_peaks.bed", "Day0_atac_peaks.bed")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day0_1007_ARC/atac_peaks.bed", "Day0_re_atac_peaks.bed")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day1_0919_ARC/atac_peaks.bed", "Day1_atac_peaks.bed")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day3_0522_ARC/atac_peaks.bed", "Day3_atac_peaks.bed")

#Day0
peaks.500 <- read.table(
  file = "Day0_atac_peaks.bed", col.names = c("chr", "start", "end")
)
#Day0re
peaks.50 <- read.table(
  file = "Day0_re_atac_peaks.bed", col.names = c("chr", "start", "end")
)
#Day1
peaks.1k <- read.table(
  file = "Day1_atac_peaks.bed", col.names = c("chr", "start", "end")
)
#Day3
peaks.5k <- read.table(
  file = "Day3_atac_peaks.bed", col.names = c("chr", "start", "end")
)
# convert to genomic ranges
gr.500 <- makeGRangesFromDataFrame(peaks.500)
gr.50 <- makeGRangesFromDataFrame(peaks.50)
gr.1k <- makeGRangesFromDataFrame(peaks.1k)
gr.5k <- makeGRangesFromDataFrame(peaks.5k)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.500, gr.50, gr.1k, gr.5k))

# Filter out bad peaks based on length 
peakwidths <- width(combined.peaks) # the range of my data is 202-3772
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# load metadata
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day0_0824_ARC/Day0_metatable_20230308.csv", "Day0_metatable_20230308.csv")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day0_1007_ARC/Day0re_metatable_20230308.csv", "Day0re_metatable_20230308.csv")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day1_0919_ARC/Day1_metatable_20230308.csv", "Day1_metatable_20230308.csv")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day3_0522_ARC/Day3_metatable_20230308.csv", "Day3_metatable_20230308.csv")

md.500 <- read.table(
  file = "Day0_metatable_20230308.csv",
  stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1) # remove the first row

md.50 <- read.table(
  file = "Day0re_metatable_20230308.csv",
  stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1) # remove the first row


md.1k <- read.table(
  file = "Day1_metatable_20230308.csv",
  stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)

md.5k <- read.table(
  file = "Day3_metatable_20230308.csv",
  stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)

# create fragment objects
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day0_0824_ARC/atac_fragments.tsv.gz", "Day0_atac_fragments.tsv.gz")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day0_1007_ARC/atac_fragments.tsv.gz", "Day0re_atac_fragments.tsv.gz")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day1_0919_ARC/atac_fragments.tsv.gz", "Day1_atac_fragments.tsv.gz")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Multipotency/Day3_0522_ARC/atac_fragments.tsv.gz", "Day3_atac_fragments.tsv.gz")

frags.500 <- CreateFragmentObject(
  path = "Day0_atac_fragments.tsv.gz",
  cells = rownames(md.500)
)

frags.50 <- CreateFragmentObject(
  path = "Day0re_atac_fragments.tsv.gz",
  cells = rownames(md.50)
)

frags.1k <- CreateFragmentObject(
  path = "Day1_atac_fragments.tsv.gz",
  cells = rownames(md.1k)
)

frags.5k <- CreateFragmentObject(
  path = "Day3_atac_fragments.tsv.gz",
  cells = rownames(md.5k)
)

#Quantify peaks in each dataset
pbmc500.counts <- FeatureMatrix(
  fragments = frags.500,
  features = combined.peaks,
  cells = rownames(md.500)
)

pbmc50.counts <- FeatureMatrix(
  fragments = frags.50,
  features = combined.peaks,
  cells = rownames(md.50)
)

pbmc1k.counts <- FeatureMatrix(
  fragments = frags.1k,
  features = combined.peaks,
  cells = rownames(md.1k)
)

pbmc5k.counts <- FeatureMatrix(
  fragments = frags.5k,
  features = combined.peaks,
  cells = rownames(md.5k)
)

# create objects
md.500[1] <- NULL
pbmc500_assay <- CreateChromatinAssay(pbmc500.counts, sep = c(":", "-"), genome = "mm10", fragments = frags.500)
pbmc500 <- CreateSeuratObject(pbmc500_assay, assay = "ATAC", meta.data=md.500)

md.50[1] <- NULL
pbmc50_assay <- CreateChromatinAssay(pbmc50.counts, sep = c(":", "-"), genome = "mm10", fragments = frags.50)
pbmc50 <- CreateSeuratObject(pbmc50_assay, assay = "ATAC", meta.data=md.50)

md.1k[1] <- NULL
pbmc1k_assay <- CreateChromatinAssay(pbmc1k.counts, sep = c(":", "-"), genome = "mm10", fragments = frags.1k)
pbmc1k <- CreateSeuratObject(pbmc1k_assay, assay = "ATAC", meta.data=md.1k)

md.5k[1] <- NULL
pbmc5k_assay <- CreateChromatinAssay(pbmc5k.counts, sep = c(":", "-"), genome = "mm10", fragments = frags.5k)
pbmc5k <- CreateSeuratObject(pbmc5k_assay, assay = "ATAC", meta.data=md.5k)

# add information to identify dataset of origin
pbmc500$dataset <- 'Day0'
pbmc50$dataset <- 'Day0re'
pbmc1k$dataset <- 'Day1'
pbmc5k$dataset <- 'Day3'

# merge all datasets, adding a cell ID to make sure cell names are unique
adata <- merge(
  x = pbmc500,
  y = list(pbmc1k, pbmc5k, pbmc50),
  add.cell.ids = c("Day0", "Day1", "Day3", "Day0re"))
adata <- adata[rowSums(adata)!=0,]

# Compute embedding and identify ATAC clusters
adata <- RunTFIDF(adata)
adata<- FindTopFeatures(adata, min.cutoff = 20)
adata<- RunSVD(adata)
adata<- RunUMAP(adata, dims = 2:50, reduction = 'lsi')
adata<- FindNeighbors(object = adata, reduction = 'lsi', dims = 2:30)
adata<- FindClusters(object = adata, verbose = FALSE, algorithm = 1, resolution = 0.3)
#DimPlot(adata, label = TRUE, group.by = 'seurat_clusters') + NoLegend()
#DimPlot(adata, label = TRUE, group.by = 'dataset') + NoLegend()

# Remove scaffold debris https://github.com/stuart-lab/signac/issues/486
main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- which(as.character(seqnames(granges(adata))) %in% main.chroms)
adata<- adata[keep.peaks, ]
adata[["ATAC"]] <- subset(adata[["ATAC"]], features = rownames(adata[["ATAC"]])[keep.peaks])

# Add RNA data and additional information
# add annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
Annotation(adata) <- annotations

# add RNA to ATAC
lfile <- connect(filename = "/Users/raycheng/Dropbox (CahanLab)/BMSC_identity/Data/Sequencing/BMSC_All_Days/BMSC_All_day_clustered_bc_041923.loom", skip.validate = TRUE)
geneNames<-lfile[["row_attrs"]][["var_names"]][]
cellNames<-colnames(adata@assays$ATAC)
expPark<- t(lfile[["matrix"]][1:length(cellNames),])
louvain<-lfile[["col_attrs"]][["leiden0"]][]
ATAC<-lfile[["col_attrs"]][["ATAC"]][]
cell_type<-lfile[["col_attrs"]][["type"]][]
time<-lfile[["col_attrs"]][["time"]][]
rownames(expPark)<-geneNames
colnames(expPark)<- cellNames
stPark<- data.frame(row.names = cellNames, cell_name=cellNames, cluster=louvain, type = cell_type, time = time)
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

adata_Seurat <- subset(adata_Seurat, features = test)

# Assign coordinates:
UMAP <- lfile[["col_attrs"]][["X_umap"]][,1:length(cellNames)]
adata[["umap"]]@cell.embeddings <- t(UMAP)
colnames(adata[["umap"]]@cell.embeddings) <- c('UMAP_1', 'UMAP_2')
rownames(adata[["umap"]]@cell.embeddings) <- colnames(adata)

# Add RNA data
adata@assays$RNA <- adata_Seurat@assays$RNA
adata@meta.data <- cbind(adata@meta.data, adata_Seurat@meta.data)
# adata@meta.data<-adata@meta.data[,-16]

## Check embedding
DimPlot(object = adata, group.by = 'type', label = TRUE, repel = TRUE) + NoLegend()
#cols = c("#4688BD", "#F29839", "#59AA4A", "#CF4B3E", '#A080C4', '#976C60') #SCANPY
#cols = c("#D62728", "#BCBD22", "#2CA02C", "#17BECF", '#1F77B4', '#E377C2') #SCANPY


## create gene activity based on ATAC
gene.activities <- GeneActivity(adata)
adata[['ATAC_RNA']] <- CreateAssayObject(counts = gene.activities)
adata<- NormalizeData(object = adata, assay = 'ATAC_RNA', normalization.method = 'LogNormalize', scale.factor = median(adata$nCount_RNA))

#Motif analysis
# Database 2: CIS-BP - a list of complete mouse pfm (better)
pfm_CISBP_raw <- readRDS("/Users/raycheng/Dropbox (CahanLab)/BMSC_identity/Data/Sequencing/BMSC_All_Days/GRN/cisBP_mouse_pfms_2021.rds")

x <- character()
for(i in 1:length(pfm_CISBP_raw@listData)){
x[i] <- pfm_CISBP_raw@listData[[i]]@name
}

motif2tf <- data.frame(motif = names(pfm_CISBP_raw@listData), tf = x, origin = "CIS-BP", gene_id = gsub("_[[:alnum:][:punct:]]*", "", names(pfm_CISBP_raw@listData)), family = NA, name = NA, symbol = NA, motif_tf = NA) %>% subset(gene_id != "XP" & gene_id != "NP")
pfm_CISBP <- subset(pfm_CISBP_raw, names(pfm_CISBP_raw@listData) %in% motif2tf$motif)

adata <- AddMotifs(object = adata, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm_CISBP)

# Computing motif activities
adata <- RunChromVAR(object = adata, genome = BSgenome.Mmusculus.UCSC.mm10)

# Check Motif sequence
DefaultAssay(adata) <- 'ATAC'
FeaturePlot(
  object = adata,
  features = head(rownames(da_peaks1)),
  pt.size = 0.1
)

MotifPlot( object = adata, motifs = head(rownames(enriched.motifs)))
install.packages('ggseqlogo')

MotifPlot(
  object = adata,
  motifs = c("Tead1"),
  assay = 'ATAC')

# Save RDS
saveRDS(adata, file = "All_day_Signac_all_20230413.RDS") 

# Additional analysis

## Error:
Error: package or namespace load failed for 'BSgenome.Mmusculus.UCSC.mm10':
 .onLoad failed in loadNamespace() for 'BSgenome.Mmusculus.UCSC.mm10', details:
  call: validObject(.Object)
  error: invalid class "TwoBitFile" object: undefined class for slot "resource" ("characterORconnection")

Solution:  it seemed whenever id try to downlood mm10, it would stop too soon. Set the timeout higher!
BiocManager::install(version = "3.12")
options(timeout = 300) 
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)

#Clustering analysis (Not useful)
DefaultAssay(adata) <- 'ATAC'
adata <- FindTopFeatures(adata, min.cutoff = 'q0')
adata<- FindNeighbors(object = adata, reduction = 'lsi', dims = 1:30)
adata<- FindClusters(object = adata, verbose = FALSE, algorithm = 1, resolution = 0.5)
DimPlot(object = adata, label = TRUE) + NoLegend()
adata@meta.data$ATAC_clusters <- adata@meta.data$ATAC_snn_res.0.2
Idents(adata) <- "type"

DefaultAssay(adata) <- "ATAC_RNA"
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
adata <- ScaleData(adata)
adata <- RunPCA(adata, features = VariableFeatures(object = adata))
ElbowPlot(adata)

adata <- FindNeighbors(adata, dims = 1:10)
adata <- FindClusters(adata, resolution = 0.5)
adata<- RunUMAP(adata, dims = 2:50, reduction = 'lsi')
# Assign coordinates:
UMAP <- lfile[["col_attrs"]][["X_umap"]][,1:length(cellNames)]
adata[["umap"]]@cell.embeddings <- t(UMAP)
colnames(adata[["umap"]]@cell.embeddings) <- c('UMAP_1', 'UMAP_2')
rownames(adata[["umap"]]@cell.embeddings) <- colnames(adata)

DimPlot(object = adata, label = TRUE, reduction = 'umap') + NoLegend()
adata@meta.data$GA_clusters <- adata@meta.data$ATAC_RNA_snn_res.0.5

# Attribution plot
library(devtools)
library(singleCellNet)
library(dplyr)
library(loomR)
library(pheatmap)

stTM <- adata@meta.data
p1 = ggplot(stTM, aes(x= stTM$time, fill = stTM$ATAC_clusters)) + geom_bar(position = 'fill', width = 0.6) + scale_y_continuous(labels = scales::percent) 
p1

p1 = ggplot(stTM, aes(x= stTM$time, fill = stTM$GA_clusters)) + geom_bar(position = 'fill', width = 0.6) + scale_y_continuous(labels = scales::percent) 
p1

# Integrate with Scanpy

write.csv(adata@meta.data, "All_days_clusters.csv")
clusters = pd.read_csv("/Users/raycheng/Desktop/All_days_clusters.csv")

clusters.index = adata_all.obs.index.copy()

adata_all.obs['ATAC_clusters'] = clusters['ATAC_clusters'].copy()
adata_all.obs['GA_clusters'] = clusters['GA_clusters'].copy()

adata_all.obs['ATAC_clusters'] = adata_all.obs['ATAC_clusters'].astype("category")
adata_all.obs['GA_clusters'] = adata_all.obs['GA_clusters'].astype("category")
