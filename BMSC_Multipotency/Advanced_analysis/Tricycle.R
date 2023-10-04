library(Seurat)
library(tricycle)
library(SingleCellExperiment)
library(ggplot2)
library(scattermore)
library(scater)
library(cowplot)

# Convert gene names into ENSEMBLE

adata <- readRDS("All_day_Signac_all_20230330.RDS")

library(org.Mm.eg.db)
DefaultAssay(adata) <- 'RNA'
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

# Alter the embedding
ori.embed <- adata2@int_colData@listData[["reducedDims"]]@listData[["tricycleEmbedding"]]
adata2@int_colData@listData[["reducedDims"]]@listData[["tricycleEmbedding"]] <- adata[["umap"]]@cell.embeddings

scater::plotReducedDim(adata2, dimred = "tricycleEmbedding",
                       colour_by = "CCStage") +
  labs(x = "UMAP1", y = "UMAP2",
       title = paste0("Projected cell cycle space (n=", ncol(adata2), ")")) +
  theme_bw(base_size = 14)

# Plot Kernel density
plot_ccposition_den(adata2$tricyclePosition,
                    adata2$type, 'type',
                    bw = 10, fig.title = "Kernel density of \u03b8", palette.v = c("#4688BD", "#F29839", "#59AA4A", "#CF4B3E", '#A080C4', '#976C60')) +
  theme_bw(base_size = 14)

plot_ccposition_den(adata2$tricyclePosition,
                    adata2$type, 'type', type = "circular",
                    bw = 10,  fig.title = "Kernel density of \u03b8", palette.v = c("#4688BD", "#F29839", "#59AA4A", "#CF4B3E", '#A080C4', '#976C60')) +
  theme_bw(base_size = 14)
