# Download the file before start
"s3://cahanlab/ray.cheng/BMSC_2023/Meta-analysis/All_BMSC_marrow_clustered_bc_20230329.loom"

# Prepare environment
library(Seurat)
library(tricycle)
library(SingleCellExperiment)
library(ggplot2)
library(scattermore)
library(scater)
library(cowplot)
library(loomR)
library(org.Mm.eg.db)

# Prepare Seurat object
lfile <- connect(filename = "All_BMSC_marrow_clustered_bc_20230329.loom", skip.validate = TRUE)
geneNames<-lfile[["row_attrs"]][["var_names"]][]
cellNames<-lfile[["col_attrs"]][["obs_names"]][]
expPark<- t(lfile[["matrix"]][1:length(cellNames),])
louvain<-lfile[["col_attrs"]][["type1"]][]
rownames(expPark)<-geneNames
colnames(expPark)<-cellNames
stPark<- data.frame(row.names = cellNames, cell_name=cellNames, cluster=louvain)
genesPark<-rownames(expPark)

adata <- CreateSeuratObject(counts = expPark, meta.data = stPark)

# Convert gene names into ENSEMBLE
master_gene_table <- mapIds(org.Mm.eg.db, keys = rownames(adata), keytype = "SYMBOL", column="ENSEMBL")
master_gene_table <- as.data.frame(master_gene_table)

adata[["RNA"]] <- AddMetaData(
  object = adata[["RNA"]],
  metadata = master_gene_table,
  col.name = 'ensembleID'
)

# Preapre SCE object
adata2 <- as.SingleCellExperiment(adata)
adata2 <- project_cycle_space(adata2, gname.type = c("SYMBOL"), species = c("mouse"))
adata2  <- estimate_cycle_position(adata2)

# Check if it works
top2a.idx <- which(rownames(adata2) == 'Top2a')
fit.l <- fit_periodic_loess(adata2$tricyclePosition,
                            assay(adata2, 'logcounts')[top2a.idx,],
                            plot = TRUE,
                       x_lab = "Cell cycle position \u03b8", y_lab = "log2(Top2a)",
                       fig.title = paste0("Expression of Top2a along \u03b8 (n=",
                                          ncol(adata2), ")"))
names(fit.l)
fit.l$fig + theme_bw(base_size = 14)


# Infer cell cycle stages
adata2 <- estimate_Schwabe_stage(adata2, gname.type = 'SYMBOL', species = 'mouse')
scater::plotReducedDim(adata2, dimred = "tricycleEmbedding",
                       colour_by = "CCStage") +
  labs(x = "Projected PC1", y = "Projected PC2",
       title = paste0("Projected cell cycle space (n=", ncol(adata2), ")")) +
  theme_bw(base_size = 14)

# Plot
p <- plot_emb_circle_scale(adata2, dimred = 1,
                           point.size = 2, point.alpha = 0.9) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
legend <- circle_scale_legend(text.size = 5, alpha = 0.9)
plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4))

# Assign coordinates:
UMAP <- lfile[["col_attrs"]][["X_umap"]][,1:length(cellNames)]
UMAP <- t(UMAP)
colnames(UMAP) <- c('UMAP_1', 'UMAP_2')
rownames(UMAP) <- colnames(adata)

# Alter the embedding
ori.embed <- adata2@int_colData@listData[["reducedDims"]]@listData[["tricycleEmbedding"]]
adata2@int_colData@listData[["reducedDims"]]@listData[["tricycleEmbedding"]] <- UMAP

scater::plotReducedDim(adata2, dimred = "tricycleEmbedding",
                       colour_by = "CCStage") +
  labs(x = "UMAP1", y = "UMAP2",
       title = paste0("Projected cell cycle space (n=", ncol(adata2), ")")) +
  theme_bw(base_size = 14)
