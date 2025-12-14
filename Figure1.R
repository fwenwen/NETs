options(stringsAsFactors = F)
library(cowplot)
library(ggpubr)
library(dplyr)
library(readr)
library(tidyr)
library(ggforce)
library(pals)
library(pheatmap)
library(scales)
library(ggthemes)
library(Seurat)
library(tidyverse)
library(ggh4x)
library(ggsci)
library(magrittr)
library(grid)
setwd("~/NETs/Figure1/")
getwd()

seuratObj <- readRDS("~/NETs/seuratObj.rds")

cols <- c(`Fibroblasts` ="#7570B3" ,
          `Endothelial cells` ="#7B4173",
          `Unciliated epithelial cells` ="#74C476",
          `Ciliated epithelial cells` = "#E7B800",
          `Macrophages`="#b10318",
          `Neutrophils`="#FD8D3C",
          `Monocytes`="#FB9A99",
          `Mast cells` = "#CE6DBD",
          `Dendritic cells`="#FED9A6",
          `NK cells` ="#377EB8",
          `T cells` = "#12a2a8",
          `B cells` = "#9c755F",
          `Plasma cells`="#767676"
)

table(seuratObj$celltype)
seuratObj$celltype <- factor(seuratObj$celltype,levels = names(cols))
DimPlot(seuratObj,group.by = "celltype",cols = cols,raster = F)+NoAxes()+ggtitle("")

ridiculous_strips <- strip_themed(
  text_x = elem_list_text(colour=c("#26A9E0","#B07AA1","#F6921E"),
                          face=c("bold","bold","bold"),size=c(10,10,10)),
  background_x = elem_list_rect(fill = c("#E8F2FC","#EBEBFB","#FFF2E7")))


pltd <- data.frame(seuratObj@reductions$umap@cell.embeddings[,1:2], 
                   cluster=as.character(seuratObj$celltype[rownames(seuratObj@reductions$umap@cell.embeddings)]))
colnames(pltd) <- c('x','y','cluster')

ggplot(data=pltd, aes(x,y)) + geom_point(aes(colour=cluster), alpha=.05, size=.2)+ 
  scale_colour_manual(values=cols) + theme_pubr() + NoLegend() +
  geom_density_2d(contour_var = "ndensity", aes(colour=cluster))

v3 <- ggplot(data=cell_eb, aes(x,y)) + geom_point(aes(colour = factor(cluster)), alpha=0.05, size=0.25)+
  geom_density_2d(contour_var = "ndensity", aes(colour=cluster))+
  facet_wrap2(vars(group),scales = "free_y",axes = "y",remove_labels = "x",strip = ridiculous_strips)+
  guides(color=guide_legend(title="location",nrow=1))+
  scale_color_manual(values = cols)+
  # scale_fill_jco()+
  coord_cartesian(clip="off")+
  theme_bw()+
  theme(panel.spacing.x = unit(0.1,"cm"),
        panel.spacing.y = unit(0.15,"cm"),
        axis.title = element_blank(),
        strip.background = element_rect(fill="grey80"),
        strip.text.x = element_text(size=9,color="black"),
        axis.text = element_text(color="black"),
        plot.margin=unit(c(0.2,0.5,0.2,0.2),units=,"cm"),
        legend.text=element_text(color="black",size=9,face="bold"),
        legend.key=element_blank(),  
        legend.title = element_blank(),
        legend.position = "top",
        legend.spacing.y=unit(0.1,'cm'),
        legend.key.width=unit(0.5,'cm'),
        legend.key.height=unit(0.5,'cm'), 
        legend.background=element_blank())+NoLegend()
v3

pdf("./Figure/Figure1_dim_facet_density.pdf",height = 3.27,width = 10)
print(v3)
dev.off()

png("./Figure/Figure1_dim_facet_density.png",height = 327,width = 1000)
print(v3)
dev.off()

# Dot plot
clusterExp <- AverageExpression(seuratObj, assays = 'RNA', slot = 'data')$RNA
clusterExp <- t(apply(clusterExp, 1, function(x){(x - min(x)) / max(x)}))

Idents(seuratObj) <- seuratObj$celltype
celltypes <- levels(seuratObj$celltype)
genes <- marker_ref$marker

percent_df <- lapply(celltypes, function(ct){
  cells <- WhichCells(seuratObj, idents = ct)
  m <- seuratObj@assays$RNA@data[genes, cells, drop=FALSE]
  
  percent <- Matrix::rowMeans(m > 0) * 100
  
  data.frame(
    gene = genes,
    celltype = ct,
    percent = percent
  )
}) %>% bind_rows()

exp_df <- clusterExp_norm[genes, celltypes, drop=FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "celltype", values_to = "Exp")


per_df <- percent_df %>%
  left_join(exp_df, by = c("gene","celltype"))


per_df <- per_df %>%
  mutate(
    gene = factor(gene, levels = genes),
    gene1 = gene,
    celltype = factor(celltype, levels = celltypes)
  )

group_labels <- c(
  rep("Fibroblasts", 5),
  rep("Endothelial cells", 3),
  rep("Unciliated epithelial cells", 4),
  rep("Ciliated epithelial cells", 3),
  rep("Macrophages", 3),
  rep("Neutrophils", 4),
  "Monocytes",
  rep("Mast cells", 4),
  rep("Dendritic cells", 3),
  rep("NK cells", 3),
  rep("T cells", 2),
  rep("B cells", 3),
  rep("Plasma cells", 4)
)

per_df$cellGroup <- factor(group_labels[match(per_df$gene, genes)],
                           levels = unique(group_labels))


pdot <-
  ggplot(per_df,aes(x = gene1,y = celltype)) +
  geom_point(aes(fill = Exp,size = percent),
             color = 'black',
             shape = 21) +
  theme_bw(base_size = 14) +
  xlab('') + ylab('') +
  scale_fill_gradient2(low = 'white',mid = '#EB1D36',high = '#990000',
                       midpoint = 0.5,
                       name = 'Mean expression') +
  scale_size(range = c(1,13)) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black'),
    aspect.ratio = 0.5,
    plot.margin = margin(1,1,1,1,"cm"),
    axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = 'italic')
  ) +
  coord_cartesian(clip = 'off')

### add segment
library(jjAnno)
pdot <- annoRect(object = pdot, annoPos = 'top',
                 aesGroup = TRUE, aesGroName = 'cellGroup',
                 yPosition = c(15.8,16.3),
                 addText = TRUE, textRot = 45, hjust = 0,
                 pCol = as.character(cols), pFill = as.character(cols),
                 rectWidth = 0.8)

p2 <- annoRect(object = pdot, annoPos = 'botomn',
               aesGroup = TRUE, aesGroName = 'cellGroup',
               yPosition = c(-2,0.22),
               pCol = as.character(cols), pFill = as.character(cols),
               rectWidth = 0.8, alpha = 0.6)

pdf("./Figure/Figure1_dot_cluster.pdf",height = 12.27,width = 21.87)
print(p2)
dev.off()

# qc
UMI_plot <- data.frame(Cell=seuratObj$sampleID, nUMI=log10(seuratObj$nCount_RNA))

p1 <- ggplot(UMI_plot, aes(x=Cell, y=nUMI)) +
  stat_boxplot(geom='errorbar', color='darkgrey', width=0.4) +
  geom_boxplot(fill=sampleCols, size=0.5, width=0.8, outlier.shape=NA) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1))

nGene_plot <- data.frame(Cell=seuratObj$sampleID, nGene=log10(seuratObj$nFeature_RNA))

p2 <- ggplot(nGene_plot, aes(x=Cell, y=nGene)) +
  stat_boxplot(geom='errorbar', color='darkgrey', width=0.4) +
  geom_boxplot(fill=sampleCols, size=0.5, width=0.8, outlier.shape=NA) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1))

p1 + p2

# feature
immu_marks <- list(
  Macrophages = c("C1QB","APOE","C1QA"),
  `Mast cells` = c("TPSAB1","CPA3","TPSB2"),
  Monocytes = c("FCN1"),
  `Dendritic cells` = c("IL3RA","IRF7"),
  Neutrophils = c("FCGR3B","CSF3R"),
  `T cells` = c("CD3D","CD3E"),
  `NK cells` = c("NKG7","NCAM1","CCL5","GZMA","CST7"),
  `B cells` = c("MS4A1","TNFRSF13C","LY9","BANK1"),
  `Plasma cells` = c("IGHA1","SSR4","MZB1","FKBP11")
)

for(x in names(immu_marks)){
  immuObj <- AddModuleScore(
    object = immuObj,
    features = list(immu_marks[[x]]),
    name = sprintf("%s_score", x)
  )
}
mapal <- colorRampPalette(c(rev(brewer.blues(8))[-(6:8)], brewer.ylorrd(5)))(256)

plist <- list()
for(m in c("Macrophages_score1","Monocytes_score1","Neutrophils_score1",
           "Mast.cells_score1","Dendritic.cells_score2","NK.cells_score1",
           "T.cells_score1","B.cells_score1","Plasma.cells_score1")){
  fp <- FeaturePlot(immuObj, reduction = "umap", features = m,
                    order = TRUE, cols = c("lightgrey", brewer.orrd(10)),
                    raster = FALSE, min.cutoff = "q10",
                    pt.size = 0.5, combine = TRUE) +
    NoLegend() + NoAxes()
  plist <- c(plist, list(fp))
}

ggarrange(plotlist = plist, legend = "none", ncol = 3, nrow = 3)

# percent
library(ggradar)

cellstat <- table(immuObj$celltype, immuObj$group)
cellstat <- t(cellstat) * 100 / colSums(cellstat)
pltd <- data.frame(group = rownames(cellstat),
                   do.call("cbind", lapply(colnames(cellstat), function(i) cellstat[, i])))
colnames(pltd)[-1] <- colnames(cellstat)

cellstat <- table(immuObj$celltype, immuObj$sample)
cellstat <- t(cellstat) * 100 / colSums(cellstat)

samples <- unique(data.frame(sample = immuObj$sample, group = immuObj$group))
rownames(samples) <- samples$sample

pltd <- cellstat %>%
  as.data.frame() %>%
  mutate(cluster = Var2, value = Freq) %>%
  mutate(group = factor(samples[Var1, "group"], levels = names(gcols)))

ggplot(pltd, aes(x = cluster, y = value, color = cluster)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "crossbar", width = 0.5) +
  facet_grid(group ~ .) +
  scale_color_manual(values = immuCols) +
  theme_pubr() + NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

data_summary <- function(data, varname, groupnames){
  sem <- function(x) sd(x, na.rm = TRUE) / sqrt(length(x))
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm = TRUE),
      se = sem(x[[col]]))
  }
  data_sum <- plyr::ddply(data, groupnames, .fun = summary_func, varname)
  plyr::rename(data_sum, c("mean" = varname))
}

pltd1 <- data_summary(pltd, varname = "value", groupnames = c("group", "cluster"))

pltd1 %>%
  ggplot(aes(x = cluster, y = value, fill = cluster)) +
  geom_bar(stat = "identity", color = NA, position = position_dodge()) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = .2,
                position = position_dodge(.9)) +
  facet_grid(group ~ .) +
  scale_fill_manual(values = immuCols) +
  theme_pubr() + NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pval2 <- setNames(lapply(levels(pltd$cluster), function(i){
  summary(aov(value ~ group, data = pltd %>% filter(cluster == i)))[[1]]$`Pr(>F)`[1]
}) %>% unlist(), levels(pltd$cluster))


seuratObj$majorCluster <- factor(seuratObj$majorCluster, levels = names(majorColor))

cellstat3 <- table(seuratObj$sampleID, seuratObj$majorCluster)
cellstat3 <- as.data.frame(cellstat3 * 100 / rowSums(cellstat3))
colnames(cellstat3) <- c("patient", "cellType", "percentage")

patient_meta <- unique(data.frame(patient = seuratObj$sampleID, group = seuratObj$group))
rownames(patient_meta) <- patient_meta$patient
cellstat3$group <- patient_meta[as.character(cellstat3$patient), "group"]


ggbarplot(
  cellstat3, x="patient", y="percentage", fill="cellType",
  color="white", add="mean_se", palette=majorColor
) + theme(axis.text.x=element_text(angle=90, hjust=1))

# 
library(tidyverse)
library(RColorBrewer)
library(scales)
library(igraph)

statdf <- data.frame()
for(i in levels(seuratObj$sampleID)){
  pvalues <- read.delim(sprintf("./out/%s/pvalues.txt", i), check.names = F)
  pvalues <- pvalues[, 12:dim(pvalues)[2]]
  statdf0 <- as.data.frame(colSums(pvalues < 0.05))
  colnames(statdf0) <- c("number")
  statdf0$indexb <- str_replace(rownames(statdf0), "^.*\\|", "")
  statdf0$indexa <- str_replace(rownames(statdf0), "\\|.*$", "")
  rankname <- sort(unique(statdf0$indexa))
  statdf0 <- statdf0[statdf0$number > 0, ]
  statdf0$sample <- i
  write.csv(
    statdf0,
    sprintf(
      "./tables/%s_number.csv",
      i
    ),
    col.names = T
  )
  statdf <- rbind(statdf, statdf0)
}
write.csv(statdf, "./tables/pair_number.csv", col.names = T)

statdf <- read.csv("./tables/pair_number.csv")

for(j in levels(seuratObj$celltype)[c(5:13)]){
  mast.list <- list()
  for(i in levels(seuratObj$celltype)[c(1:4)]){
    statdf_select <- statdf[statdf$indexb %in% i & statdf$indexa %in% j, ]
    library(pals)
    statdf_select1 <- statdf_select %>%
      group_by(sample) %>%
      summarise(number1 = sum(number))
    statdf_select1$group <- factor(
      sampleGroup[as.character(statdf_select1$sample), "group"],
      levels = levels(seuratObj$group)
    )
    statdf_select1$log_number <- log(statdf_select1$number1)
    my_comparisons <- list(c("EN","Normal"), c("EH","EN"), c("EH","Normal"))
    mast.list[[i]] <- ggboxplot(
      statdf_select1,
      x = "group",
      y = "log_number",
      fill = "group",
      color = "black",
      palette = gcols
    ) +
      ggtitle(paste(j, i, sep = "_")) +
      labs(x = "group", y = "Interaction count (log)") +
      NoLegend() +
      stat_compare_means(
        comparisons = my_comparisons,
        method = "wilcox.test",
        method.args = list(alternative = "greater")
      )
  }
  pdf(sprintf("./out/%s_box.pdf", j), height = 4.27, width = 14.27)
  print(ggarrange(plotlist = mast.list, ncol = 4, nrow = 1))
  dev.off()
}

statdf$pairs <- sub("\\d+$", "", statdf$X)
statdf$group <- sub("\\d+$", "", statdf$sample)

result2 <- statdf %>%
  group_by(pairs, group) %>%
  summarize(mean_X_group = mean(number)) %>%
  ungroup()

matrix_data <- spread(result2, key = pairs, value = mean_X_group) %>% as.matrix()
rownames(matrix_data) <- matrix_data[, 1]
matrix_data <- matrix_data[, -1]

saveRDS(matrix_data, "./tables/matrix_data.rds")

pdf("out/Neu.CCI.chat.pdf", width = 6, height = 2)
op <- par(mfrow = c(1, 3), mar = rep(0.5, 4), xpd = TRUE)
c64 <- levels(seuratObj$celltype)[c(1:4, 6)]
for (i in levels(seuratObj$group)) {
  mat <- matrix(0, nrow = length(c64), ncol = length(c64), dimnames = list(c64, c64))
  mat[, "Neutrophils"] <- matrix_data[i, sprintf("Neutrophils|%s", c64)] %>% as.numeric()
  mat["Neutrophils", "Neutrophils"] <- 0
  CellChat::netVisual_circle(
    mat,
    color.use = cols[c64],
    weight.scale = T,
    edge.weight.max = max(10),
    title.name = i
  )
}
par(op)
dev.off()

