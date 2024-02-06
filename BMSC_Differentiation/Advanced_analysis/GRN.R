# Without Pseudotime
R
library(igraph)
library(qgraph)
library(loomR)
library(gam)
library(devtools)
library(singleCellNet)
library(RColorBrewer)
library(pheatmap)
library(epoch)
library(minet)

mydate<-utils_myDate()

loadLoomExpDiffMap<-function# load a loom object containing expression data
(path,
  cellNameCol='obs_names',
  xname='louvain'
  ){
  lfile <- connect(filename = path, skip.validate = TRUE)
  geneNames<-lfile[["row_attrs"]][["var_names"]][]
  cellNames<-lfile[["col_attrs"]][["obs_names"]][]
  expMat<- t(lfile[["matrix"]][1:length(cellNames),])
  rownames(expMat)<-geneNames
  colnames(expMat)<-cellNames
  cluster <- lfile[['col_attrs']][[xname]][]
  sampTab <- data.frame(cell_name=cellNames, cluster = cluster)
  lfile$close_all()
  rownames(sampTab) = as.vector(sampTab$cell_name)
  list(expDat = expMat, sampTab = sampTab)
}

#Load data
mmTFs<-utils_loadObject("/Users/raycheng/Dropbox (CahanLab)/BMSC_Meta/Data/Reference/mmTFs_123019.rda")
list12<-loadLoomExpDiffMap("BMSC_AD_OS_CH_clustered_bc_20231204.loom", xname='type')

expDat<-list12[['expDat']]
sampTab<-list12[['sampTab']] 

# Expanded
sampTab<-sampTab[which(sampTab$cluster %in% c("BMSC1")),]

# Adipogenic
sampTab<-sampTab[which(sampTab$cluster %in% c("BMSC1", 'BMSC2', 'BMSC3', 'BMSC4')),]
# Chondrogenic
sampTab<-sampTab[which(sampTab$cluster %in% c("BMSC1", 'BMSC2', 'BMSC5')),]
# Osteogenic
sampTab<-sampTab[which(sampTab$cluster %in% c("BMSC1", 'BMSC6', 'BMSC7', 'BMSC8')),]

# Remove the redundant clusteres, Remove expression level = 0
expDat <- expDat[,rownames(sampTab)]
expDat= expDat[rowSums(expDat)>0, ] 

# Reconstruct and perform optional crossweighting
system.time(grnDF <- reconstructGRN(expDat, mmTFs, method="pearson", zThresh=3))

# Save
write.csv(grnDF, "All_diff_stem_20230501.csv", sep= "\t", row.names = TRUE, col.names = TRUE)

compute_pagerank2<-function(dynnet,weight_column="zscore",directed_graph=FALSE){
    df<-dynnet
    df<-df[,c("TF","TG",weight_column)]
    colnames(df)<-c("TF","TG","weight")
    tfnet<-graph_from_data_frame(df,directed=directed_graph)
    pagerank<-data.frame(page_rank(tfnet,directed=directed_graph)$vector)
    colnames(pagerank)<-c("page_rank")
    pagerank$gene<-rownames(pagerank)
    pagerank<-pagerank[,c("gene","page_rank")]
    pagerank<-pagerank[order(pagerank$page_rank,decreasing=TRUE),]
    pagerank$is_regulator<-FALSE
    pagerank$is_regulator[pagerank$gene %in% unique(df$TF)]<-TRUE
    ranks<-pagerank
ranks
}

gene_rank<-compute_pagerank2(grnDF,weight_column="zscore")


# Adipogenic
grnDF <- read.csv("All_diff_adipo_20230501.csv", row.name = 1)

# Select features
# TF list
tf_list <- c("Cebpb", 'Ppara', 'Stat5b', 'Pparg', 'Nfia', 'Ebf1', 'Nr2c2', 'Atf4', 'Tef', 'Ppard', 'Ctcf', 'Nfix', 'Cebpa', 'Zeb2', 'Ebf1', 'Sp7', 'Runx1', 'Runx3', 'Runx2', 'Nfic', 'Ar', 'Vdr', 'Prrx1')

tglist1 <- c('Ibsp', 'Mmp13', 'Gas6', 'Plxna2', 'Igfbp5', 'Col3a1', 'Col6a1', 'Igfbp4', 'Alpl')
tglist2 <- c('Cxcl12', 'Kitl', 'Tgfbr3', 'Negr1', 'Vcam1', 'Plpp3', 'Nrp1', 'Lepr', 'Apoe', 'Ghr')
tglist3 <- c('Fabp4', 'Adipoq', 'Lpl', 'Cd36', 'Acsl1', 'Pde3b', 'Fabp5', 'Adipor2')

tg_list <- c(tglist1, tglist2, tglist3)

gene_list <- c(tf_list, tg_list)

grnDF_adipo <- grnDF[which(grnDF$TF %in% gene_list & grnDF$TG %in% gene_list),]
#grnDF_adipo <- grnDF_adipo[which(grnDF_adipo$corr > 0.1 | grnDF_adipo$corr < -0.1),]

final_DF <- grnDF_adipo
final_DF$name <- paste(pmax(final_DF$TG, final_DF$TF), pmin(final_DF$TG, final_DF$TF), sep = "_")
final_DF <- final_DF[!duplicated(final_DF$name), ]

# Painting
TG <- as.data.frame(final_DF$TG)
TF <- as.data.frame(final_DF$TF)
colnames(TG) <- "node"
colnames(TF) <- "node"
nodes <- rbind(TG, TF)
nodes  <- unique(nodes$node)
nodes <- as.data.frame(nodes)

# Assign node types
nodes$type <- 0
nodes[which(nodes$nodes %in% tglist1 ), ]$type <- 2
nodes[which(nodes$nodes %in% tglist2 ), ]$type <- 3
nodes[which(nodes$nodes %in% tglist3 ), ]$type <- 4
nodes[which(nodes$nodes %in% mmTFs), ]$type <- 1

# Assign edge types
links <- final_DF
links$edge <- 1
links[which(links$corr <0),]$edge <- 2 # represent inhibitory

# Build network
net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)

# Plot node text color
colrs <- c("#000000", "#AA7942", "#AA7942", "#AA7942") 
V(net)$color <- colrs[V(net)$type]
V(net)$label.cex = 0.5
#V(net)$label.cex = (V(net)$size+1)/25

# Node color based on groups
# c('#E97E73', '#CA9734', '#9EAF3C', '#7DC874', '#73C9B0', '#64BEE1', '#7FA5F9', '#CE78F4', '#F292D0')
# V(net)$group <- 1
colrs <- c("#F6D8D8", "#CA9734", "#9EAF3C", "#7DC874") 
V(net)$color2 <- colrs[V(net)$type]

# Plot edge color
colrs2 <- c("#FC2B2D", "#106BFF")
E(net)$color <- colrs2[E(net)$edge]

# pick layout
#graph_attr(net, "layout") <- layout_with_gem
#graph_attr(net, "layout") <- layout_with_lgl
graph_attr(net, "layout") <- layout_with_kk
#graph_attr(net, "layout") <- layout_with_fr
#graph_attr(net, "layout") <- layout_with_graphopt
#graph_attr(net, "layout") <- layout_with_dh

png(file = '/Users/raycheng/Desktop/plot1.png', width = 5, height = 5, units = 'in', res = 300)
par(bg=NA)
plot.igraph(net,
edge.curved=curve_multiple(net), 
vertex.label.color = adjustcolor(V(net)$color, alpha.f = 1), 
vertex.label.font=1, 
vertex.label.degree= -pi/2, 
vertex.frame.color = "transparent", 
vertex.color = adjustcolor(V(net)$color2, alpha.f = 0.5), 
edge.color = adjustcolor(E(net)$color, alpha.f = .25), 
edge.width=1)
dev.off()

# Chondroenic
grnDF <- read.csv("All_diff_chondro_20230501.csv", row.name = 1)

# Select features
# TF list
tf_list <- c('Sp7', 'Runx1', 'Runx3', 'Runx2', 'Nfic', 'Ar', 'Vdr', 'Prrx1', 'Cbfb', 'Tead1', 'Prrx1', 'Ybx1', 'Twist1', 'Ebf1', 'Ppara', 'Pparg', 'Ebf3', 'Ebf2')

tglist1 <- c('Ibsp', 'Ogn', 'Mmp13', 'Gpx3', 'Ptx3', 'Gas6', 'Plxna2', 'Igfbp5', 'Col3a1', 'Col6a1', 'Igfbp4', 'Alpl', 'Omd', 'Aspn', 'Chrdl1', 'Col8a1', 'Col8a2')
tglist2 <- c('Cxcl12', 'Malat1', 'Zeb2', 'Runx1', 'Epas1', 'Tgfbr3', 'Mark1', 'Cald1', 'Tln2', 'Kitl', 'Hp', 'Tnc', 'Fgf7')

tg_list <- c(tglist1, tglist2)

gene_list <- c(tf_list, tg_list)

grnDF_chondro <- grnDF[which(grnDF$TF %in% gene_list & grnDF$TG %in% gene_list),]
#grnDF_adipo <- grnDF_adipo[which(grnDF_adipo$corr > 0.1 | grnDF_adipo$corr < -0.1),]

final_DF <- grnDF_chondro
final_DF$name <- paste(pmax(final_DF$TG, final_DF$TF), pmin(final_DF$TG, final_DF$TF), sep = "_")
final_DF <- final_DF[!duplicated(final_DF$name), ]

# Painting
TG <- as.data.frame(final_DF$TG)
TF <- as.data.frame(final_DF$TF)
colnames(TG) <- "node"
colnames(TF) <- "node"
nodes <- rbind(TG, TF)
nodes  <- unique(nodes$node)
nodes <- as.data.frame(nodes)

# Assign node types
nodes$type <- 0
nodes[which(nodes$nodes %in% tglist1 ), ]$type <- 2
nodes[which(nodes$nodes %in% tglist2 ), ]$type <- 3
nodes[which(nodes$nodes %in% mmTFs), ]$type <- 1

# Assign edge types
links <- final_DF
links$edge <- 1
links[which(links$corr <0),]$edge <- 2 # represent inhibitory

# Build network
net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)

# Plot node text color
colrs <- c("#000000", "#AA7942", "#AA7942") 
V(net)$color <- colrs[V(net)$type]
V(net)$label.cex = 0.5
#V(net)$label.cex = (V(net)$size+1)/25

# Node color based on groups
# c('#E97E73', '#CA9734', '#9EAF3C', '#7DC874', '#73C9B0', '#64BEE1', '#7FA5F9', '#CE78F4', '#F292D0')
# V(net)$group <- 1
colrs <- c("#F6D8D8", "#CA9734", "#73C9B0") 
V(net)$color2 <- colrs[V(net)$type]

# Plot edge color
colrs2 <- c("#FC2B2D", "#106BFF")
E(net)$color <- colrs2[E(net)$edge]

# pick layout
#graph_attr(net, "layout") <- layout_with_gem
#graph_attr(net, "layout") <- layout_with_lgl
graph_attr(net, "layout") <- layout_with_kk
#graph_attr(net, "layout") <- layout_with_fr
#graph_attr(net, "layout") <- layout_with_graphopt
#graph_attr(net, "layout") <- layout_with_dh

png(file = '/Users/raycheng/Desktop/plot2.png', width = 5, height = 5, units = 'in', res = 300)
par(bg=NA)
plot.igraph(net,
edge.curved=curve_multiple(net), 
vertex.label.color = adjustcolor(V(net)$color, alpha.f = 1), 
vertex.label.font=1, 
vertex.label.degree= -pi/2, 
vertex.frame.color = "transparent", 
vertex.color = adjustcolor(V(net)$color2, alpha.f = 0.5), 
edge.color = adjustcolor(E(net)$color, alpha.f = .25), 
edge.width=1)
dev.off()


# Osteogenic
grnDF <- read.csv("All_diff_osteo_20230501.csv", row.name = 1)

# Select features
# TF list
tf_list <- c("Fos", 'Fosb', 'Fosl2', 'Ebf1', 'Glis3', 'Jund', 'Bach2', 'Junb', 'Bach1', 'Nfkb2', 'Nfkb1', 'Rela', 'Relb', 'Runx1', 'Runx2', 'Mef2a', 'Glis3', 'Zeb1', 'Pbx1', 'Ebf3', 'Cebpb', 'Pparg', 'Mef2a')

tglist1 <- c('Gzmd', 'Gzme', 'Ccn2', 'Tagln', 'Timp1', 'S100a6', 'S100a11', 'Lgals1', 'Tpm2')
tglist2 <- c('Col11a1', 'Col11a2', 'Fgfr1', 'Fgfr2', 'Edil3', 'Plcb1', 'Exoc4', 'Fndc3b', 'Tnc')
tglist3 <- c('Cxcl12', 'Lifr', 'Zeb2', 'Ghr', 'Pdzrn4', 'Kitl', 'Negr1', 'Limch1', 'Shroom3', 'Ibsp')

tg_list <- c(tglist1, tglist2, tglist3)

gene_list <- c(tf_list, tg_list)

grnDF_osteo<- grnDF[which(grnDF$TF %in% gene_list & grnDF$TG %in% gene_list),]
#grnDF_adipo <- grnDF_adipo[which(grnDF_adipo$corr > 0.1 | grnDF_adipo$corr < -0.1),]

final_DF <- grnDF_osteo
final_DF$name <- paste(pmax(final_DF$TG, final_DF$TF), pmin(final_DF$TG, final_DF$TF), sep = "_")
final_DF <- final_DF[!duplicated(final_DF$name), ]
#final_DF <- final_DF[which(final_DF$zscore > 6),]

# Painting
TG <- as.data.frame(final_DF$TG)
TF <- as.data.frame(final_DF$TF)
colnames(TG) <- "node"
colnames(TF) <- "node"
nodes <- rbind(TG, TF)
nodes  <- unique(nodes$node)
nodes <- as.data.frame(nodes)

# Assign node types
nodes$type <- 0
nodes[which(nodes$nodes %in% tglist1 ), ]$type <- 2
nodes[which(nodes$nodes %in% tglist2 ), ]$type <- 3
nodes[which(nodes$nodes %in% tglist3 ), ]$type <- 4
nodes[which(nodes$nodes %in% mmTFs), ]$type <- 1

# Assign edge types
links <- final_DF
links$edge <- 1
links[which(links$corr <0),]$edge <- 2 # represent inhibitory

# Build network
net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)

# Plot node text color
colrs <- c("#000000", "#AA7942", "#AA7942", "#AA7942") 
V(net)$color <- colrs[V(net)$type]
V(net)$label.cex = 0.5
#V(net)$label.cex = (V(net)$size+1)/25

# Node color based on groups
# c('#E97E73', '#CA9734', '#9EAF3C', '#7DC874', '#73C9B0', '#64BEE1', '#7FA5F9', '#CE78F4', '#F292D0')
# V(net)$group <- 1
colrs <- c("#F6D8D8", "#64BEE1", "#7FA5F9", "#CE78F4") 
V(net)$color2 <- colrs[V(net)$type]

# Plot edge color
colrs2 <- c("#FC2B2D", "#106BFF")
E(net)$color <- colrs2[E(net)$edge]

# pick layout
#graph_attr(net, "layout") <- layout_with_gem
#graph_attr(net, "layout") <- layout_with_lgl
graph_attr(net, "layout") <- layout_with_kk
#graph_attr(net, "layout") <- layout_with_fr
#graph_attr(net, "layout") <- layout_with_graphopt
#graph_attr(net, "layout") <- layout_with_dh

png(file = '/Users/raycheng/Desktop/plot3.png', width = 5, height = 5, units = 'in', res = 300)
par(bg=NA)
plot.igraph(net,
edge.curved=curve_multiple(net), 
vertex.label.color = adjustcolor(V(net)$color, alpha.f = 1), 
vertex.label.font=1, 
vertex.label.degree= -pi/2, 
vertex.frame.color = "transparent", 
vertex.color = adjustcolor(V(net)$color2, alpha.f = 0.5), 
edge.color = adjustcolor(E(net)$color, alpha.f = .25), 
edge.width=1)
dev.off()
