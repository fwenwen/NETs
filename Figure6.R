library(dplyr)
library(Seurat)
library(patchwork)
library(pals)
library(pheatmap)
library(ggplot2)
library(ggthemes)

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

raincloudPlot <- function(pdata, tt, cols=NULL){
  if(nrow(pdata) < 5){
    p <- NULL
  }else{
    p <- ggplot(pdata, aes(x = group, y = value, fill = group, color = group)) +
      ggtitle(tt) +
      ylab(NULL) + xlab(NULL) +
      theme_cowplot() +
      scale_shape_identity() +
      theme_pubr() +
      # theme(legend.position = "none",
      #       plot.title = element_text(size = 20),
      #       axis.title = element_text(size = 15),
      #       axis.text = element_text(size = 15),
      #       axis.text.x = element_text(angle = 0, 
      #                                  hjust = 0,
      #                                  vjust = 0)) +
      scale_color_manual(values=cols) +
      scale_fill_manual(values=cols) +
      geom_point(position = position_jitter(0.15), 
                 size = 2, 
                 alpha = 1, 
                 aes(shape = 16)) +
      geom_flat_violin(position = position_nudge(x = 0.2, y = 0),
                       adjust = 2,
                       alpha = 0.6, 
                       trim = TRUE, 
                       scale = "width") +
      geom_boxplot(aes(x = as.numeric(group) + 0.2, y = value), 
                   notch = FALSE, 
                   width = 0.1, 
                   varwidth = FALSE, 
                   outlier.shape = NA, 
                   alpha = 0.3, 
                   colour = "black", 
                   show.legend = FALSE) #+
      # stat_compare_means()
  }
  return(p)
}

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}


geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
            
          },
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)) #利用transform函数为数据框mydata增加数据
            newdata <- rbind(plyr::arrange(transform(data, x = xmaxv), -y),plyr::arrange(transform(data, x = xminv), y))
            newdata_Polygon <- rbind(newdata, newdata[1,])
            newdata_Polygon$colour<-NA
            newdata_Path <- plyr::arrange(transform(data, x = xmaxv), -y)
            ggplot2:::ggname("geom_flat_violin", grobTree(
              GeomPolygon$draw_panel(newdata_Polygon, panel_scales, coord),
              GeomPath$draw_panel(newdata_Path, panel_scales, coord))
            )
          },
          draw_key = draw_key_polygon,
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          required_aes = c("x", "y")
  )

seuratObj <- readRDS("~/mouse/monkeyObj.rds")

Idents(seuratObj) <- seuratObj$group
avgExp0 <- AverageExpression(seuratObj, assays = 'SCT')$SCT 
head(avgExp0)
colSums(avgExp0)
avgExp <- avgExp0[rowSums(avgExp0) > 0,]
sampleCors <- cor(avgExp, method = 'p')
head(avgExp)
pheatmap::pheatmap(sampleCors, cluster_rows = T, cluster_cols = T, 
                   cellwidth = 20, cellheight = 20)
# 
library(foreach)
library(enrichplot)
library(clusterProfiler)
library(org.Mfascicularis.eg.db)
minpct <- .1
gseOut <- foreach(cls=levels(as.factor(top.markers$cluster))) %do% {
  submarkers <- top.markers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
  geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$gene), decreasing = T)
  ego <- gseGO(geneList     = geneList,
               OrgDb        = org.Mfascicularis.eg.db,
               ont          = "BP",
               keyType      = "GID",
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 1,
               verbose      = FALSE)
  ego@result <- cbind(ego@result, cluster=cls)
  ego
}
gseResult <- do.call("rbind", lapply(gseOut, function(x){x@result}))

table(gseResult$cluster)
topPath <- gseResult %>% dplyr::mutate(logFDR=-log10(pvalue)) %>% 
  dplyr::filter(logFDR > 1.5 & enrichmentScore>0)

SLC12A2_topPath <- topPath[grep("ENSMFAG00000001384",topPath$core_enrichment),]

dim(SLC12A2_topPath)
gseMat <- -log10(reshape2::acast(gseResult[gseResult$enrichmentScore>0, ], 
                                 ID~cluster, value.var='p.adjust'))
dim(gseMat)
gseMat[1:3,1:3]


gseMat[is.na(gseMat)] <- 0

gseMat0 <- gseMat[SLC12A2_topPath$ID,names(gcols)]
gseMat0[1:3,1:3]
dim(gseMat0)
library(GOfuncR)
labels_row <- sprintf("%s %s",get_names(rownames(gseMat0))$go_id, get_names(rownames(gseMat0))$go_name)
# labels_row[apply(gseMat0, 1, function(x){all(x<2.6) | sum(x>2.6) > 2})] <- ""
labels_row[nchar(labels_row) > 40] <- ""
gseMat0[1:3,1:3]
mapal <- colorRampPalette(brewer.ylorrd(5))(256)
annotation_col <- data.frame(group=colnames(gseMat0))
rownames(annotation_col) <- annotation_col$group
ann_colors <- list(group = gcols)

pheatmap::pheatmap(gseMat0, cluster_cols = F, cluster_rows = T,
                   color=mapal, labels_row = labels_row,
                   annotation_col = annotation_col, annotation_colors = ann_colors,
                   fontsize_row = 3, border_color = NA)

EH_gse <- gseOut[[3]]
gos_EH <- SLC12A2_topPath$ID
EH_gse@result[EH_gse@result$ID %in% gos_EH,]$Description
# 
y <- EH_gse@result[EH_gse@result$ID %in% gos_EH,]$Description
SLC12A2_topPath.gene <- unique(unlist(strsplit(SLC12A2_topPath[SLC12A2_topPath$Description %in% y,]$core_enrichment,"/")))
# 
library(enrichplot)
cnetplot(EH_gse, node_label="all",foldChange = foldchange,showCategory = y,colorEdge = TRUE)
# 
Idents(epiObj) <- epiObj$group
clusterExp <- AverageExpression(epiObj, assays = 'SCT', slot = 'data')$SCT
table(rownames(clusterExp) %in% SLC12A2_topPath.gene)
markerExp <- t(apply(clusterExp[as.character(SLC12A2_topPath.gene), ], 1, scale))
colnames(markerExp) <- colnames(clusterExp)

labels_row <- top.markers_EH$name
annotation_col <- data.frame(group=colnames(markerExp))
rownames(annotation_col) <- annotation_col$group
ann_colors <- list(group = gcols)

mapal <- colorRampPalette(c(rev(brewer.blues(8))[-(6:8)],brewer.ylorrd(5)))(256)
library(ComplexHeatmap)
ha = rowAnnotation(foo = anno_mark(at = match(labels_row,rownames(markerExp)),
                                   labels = labels_row))

ComplexHeatmap::pheatmap(markerExp, cluster_cols = F, cluster_rows = T, 
                               color=mapal, labels_row = rep("",nrow(markerExp)),
                               right_annotation = ha,angle_col = "45",
                               annotation_col = annotation_col, annotation_colors = ann_colors,
                               fontsize_row = 6, border_color = NA)

