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
#library(BSgenome.Hsapiens.UCSC.hg38) # human
library(chromVARmotifs)
library(Pando)
library(magrittr)

#plan("multiprocess", workers = 4)
#options(future.globals.maxSize = 16000 * 1024^2)
 
# Prepare bed files
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/AD/atac_peaks.bed", "AD_atac_peaks.bed")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/OS/atac_peaks.bed", "OS_atac_peaks.bed")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/CH/atac_peaks.bed", "CH_atac_peaks.bed")

#AD
peaks.500 <- read.table(file = "AD_atac_peaks.bed", col.names = c("chr", "start", "end"))
#OS
peaks.1k <- read.table(file = "OS_atac_peaks.bed", col.names = c("chr", "start", "end"))
#CH
peaks.5k <- read.table(file = 'CH_atac_peaks.bed', col.names = c("chr", "start", "end"))

# convert to genomic ranges
gr.500 <- makeGRangesFromDataFrame(peaks.500)
gr.1k <- makeGRangesFromDataFrame(peaks.1k)
gr.5k <- makeGRangesFromDataFrame(peaks.5k)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.500, gr.1k))

# Filter out bad peaks based on length 
peakwidths <- width(combined.peaks) # the range of my data is 107-2509
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# load metadata
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/AD/AD_metatable.csv", "AD_metatable.csv")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/OS/OS_metatable_20230328.csv", "OS_metatable.csv")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/CH/CH_metatable_20240425.csv", "CH_metatable.csv")

md.500 <- read.table(
  file = "AD_metatable.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
) # remove the first row

md.1k <- read.table(
  file = "OS_metatable.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)

md.5k <- read.table(file = "CH_metatable.csv", stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)


# create fragment objects
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/AD/atac_fragments.tsv.gz", "AD_atac_fragments.tsv.gz")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/OS/atac_fragments.tsv.gz", "OS_atac_fragments.tsv.gz")
download.file("s3://cahanlab/ray.cheng/BMSC_2023/BMSC_Differentiation/CH/atac_fragments.tsv.gz", "CH_atac_fragments.tsv.gz")

frags.500 <- CreateFragmentObject(
  path = "AD_atac_fragments.tsv.gz",
  cells = rownames(md.500)
)

frags.1k <- CreateFragmentObject(
  path = "OS_atac_fragments.tsv.gz",
  cells = rownames(md.1k)
)

frags.5k <- CreateFragmentObject(path = "CH_atac_fragments.tsv.gz", cells = rownames(md.5k))

#Quantify peaks in each dataset
pbmc500.counts <- FeatureMatrix(
  fragments = frags.500,
  features = combined.peaks,
  cells = rownames(md.500)
)

pbmc1k.counts <- FeatureMatrix(
  fragments = frags.1k,
  features = combined.peaks,
  cells = rownames(md.1k)
)

pbmc5k.counts <- FeatureMatrix(fragments = frags.5k, features = combined.peaks, cells = rownames(md.5k))

# create objects
md.500[1] <- NULL
pbmc500_assay <- CreateChromatinAssay(pbmc500.counts, sep = c(":", "-"), genome = "mm10", fragments = frags.500)
pbmc500 <- CreateSeuratObject(pbmc500_assay, assay = "ATAC", meta.data=md.500)

md.1k[1] <- NULL
pbmc1k_assay <- CreateChromatinAssay(pbmc1k.counts, sep = c(":", "-"), genome = "mm10", fragments = frags.1k)
pbmc1k <- CreateSeuratObject(pbmc1k_assay, assay = "ATAC", meta.data=md.1k)

md.5k[1] <- NULL
pbmc5k_assay <- CreateChromatinAssay(pbmc5k.counts, sep = c(":", "-"), genome = "mm10", fragments = frags.5k)
pbmc5k <- CreateSeuratObject(pbmc5k_assay, assay = "ATAC", meta.data=md.5k)

# add information to identify dataset of origin
pbmc500$dataset <- 'AD'
pbmc1k$dataset <- 'OS'
pbmc5k$dataset <- 'CH'

# merge all datasets, adding a cell ID to make sure cell names are unique
adata <- merge(
  x = pbmc500,
  y = list(pbmc1k, pbmc5k),
  add.cell.ids = c("AD", "OS", "CH"))
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

# Remove scaffold debris
 https://github.com/stuart-lab/signac/issues/486
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
lfile <- connect(filename = "/Users/raycheng/Dropbox (CahanLab)/BMSC_identity/Data/Sequencing/BMSC_All_Diff/BMSC_AD_OS_CH_clustered_bc_20230425.loom", skip.validate = TRUE)
geneNames<-lfile[["row_attrs"]][["var_names"]][]
cellNames<-colnames(adata@assays$ATAC)
expPark<- t(lfile[["matrix"]][1:length(cellNames),])
cell_type<-lfile[["col_attrs"]][["type"]][]
#ATAC<-lfile[["col_attrs"]][["ATAC"]][]
time<-lfile[["col_attrs"]][["time"]][]
lineage <-lfile[["col_attrs"]][["lineage"]][]
rownames(expPark)<-geneNames
colnames(expPark)<-cellNames
stPark<- data.frame(row.names = cellNames, cell_name=cellNames, type=cell_type,  time=time, lineage = lineage)
genesPark<-rownames(expPark)

#row.names(ad$obs) <- ad$obs$new_name
#row.names(ad$var) <- ad$var$var_names
adata_Seurat <- CreateSeuratObject(counts = expPark, meta.data = stPark)


# Remove Rik genes: For GRN
test <- rownames(adata_Seurat)

# Remove numeric genes: For GRN
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
# adata@meta.data<-adata@meta.data[,-18]

# Add fragments data
# Run merge ATAC section (Merge ATAC - fragments method)
# adata@assays$ATAC@fragments <- combined@assays$ATAC@fragments
# DefaultAssay(adata) <- "ATAC"

## Check embedding
#DefaultAssay(adata) <- "RNA"
#adata@meta.data$cluster <- as.factor(adata@meta.data$cluster)
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
  motifs = c("Pparg"),
  assay = 'ATAC')

# Save RDS
saveRDS(adata, file = "All_diff_Signac_all_20230425.RDS") 
