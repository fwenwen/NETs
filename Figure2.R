Idents(seuratObj) <- seuratObj$celltype
table(seuratObj$celltype)
NeuObj <- subset(seuratObj,idents="Neutrophils")
NeuObj@meta.data <- droplevels(NeuObj@meta.data)
DimPlot(NeuObj)
NeuObj <- NormalizeData(NeuObj) %>% FindVariableFeatures() %>% 
  ScaleData() %>% RunPCA(verbose=FALSE)
NeuObj <- RunUMAP(NeuObj, reduction = "harmony", dims = 1:30, min.dist=0.5)
NeuObj <- FindNeighbors(NeuObj, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution=0.3)
NeuObj <- FindClusters(NeuObj, resolution = 0.5,graph.name = "RNA_snn")
NeuObj$cluster <- factor(paste0("Neus",as.numeric(NeuObj$RNA_snn_res.0.5)),
                            levels = paste0("Neus",1:5))
#
pltd <- data.frame(NeuObj@reductions$umap@cell.embeddings[,1:2], 
                   cluster=as.character(NeuObj$cluster[rownames(NeuObj@reductions$umap@cell.embeddings)]))
colnames(pltd) <- c('x','y','cluster')

ggplot(data=pltd, aes(x,y)) + geom_point(aes(colour = factor(cluster)), alpha=0.05, size=0.1)+
  geom_density_2d(contour_var = "ndensity", aes(colour=cluster))+
  # facet_wrap2(vars(group),scales = "free_y",axes = "y",remove_labels = "x",strip = ridiculous_strips)+
  guides(color=guide_legend(title="location",nrow=1))+
  scale_color_manual(values = ncols)+
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

# percent
ncellstat <- table(NeuObj$group, NeuObj$cluster)
ncellstat <- as.data.frame(ncellstat * 100 / rowSums(ncellstat))
colnames(ncellstat) <- c("group", "cellType", "percentage")
ggplot(ncellstat, aes(fill=cellType, x=group, y=percentage)) + 
  geom_bar(position="fill", stat="identity",color="white") + theme_pubr() +
  scale_fill_manual(values = ncols)+coord_flip()

ggplot(data = ncellstat,
       aes(x=group, y=percentage, fill = cellType ,
           stratum = cellType , alluvium = cellType)) +
  geom_stratum(width = 0.5, color=NA) + 
  geom_alluvium(width = 0.5, curve_type = "line", 
                color="black", alpha=0.5, linewidth=0.25) +
  scale_fill_manual(values = ncols,
                    name = "Cell Type") +
  scale_x_discrete(expand = expansion(add = c(0.1, 0.1))) +
  # scale_y_continuous(name = "percentage",
  #                    limits = c(-0.01, 1.01), 
  #                    expand = c(0, 0),
  #                    breaks = seq(0, 1, 0.25),
  #                    labels = str_c(seq(0, 100, 25), "%")) +
  # coord_flip() +
  # labs(title = "Line stacked percentage histogram",
  #      subtitle = "plot charts with ggalluvial") +
  hrbrthemes::theme_ipsum(base_family = "serif") +
  theme(panel.grid = element_line(linetype = "dashed"),
        plot.background = element_rect(fill = "white", color = "white"),
        plot.margin = margin(10, 10, 10, 10))

# expression
top.markers <- read.delim("./top.markers.txt",header = T)
top20 <-  top.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=50, wt=avg_log2FC)
Idents(NeuObj) <- NeuObj$cluster
clusterExp <- AverageExpression(NeuObj, assays = 'SCT', slot = 'data')$SCT
markerExp <- t(apply(clusterExp[as.character(unique(c(top20$gene))), ], 1, scale))
colnames(markerExp) <- colnames(clusterExp)
Neu_top5 <- top20 %>% dplyr::group_by(cluster) %>% dplyr::top_n(5,wt = avg_log2FC)
Neu_top5
labels_row <- unique(as.character(Neu_top5$gene))
table(labels_row %in% rownames(markerExp))

annotation_col <- data.frame(cluster=colnames(markerExp))
rownames(annotation_col) <- annotation_col$cluster
ann_colors <- list(cluster = ncols)
mapal <- colorRampPalette(c(rev(brewer.blues(8))[-(6:8)],brewer.ylorrd(5)))(256)

library(ComplexHeatmap)
ha = rowAnnotation(foo = anno_mark(at = match(labels_row,rownames(markerExp)),
                                   labels = labels_row))

ComplexHeatmap::pheatmap(markerExp, cluster_cols = F, cluster_rows = F, 
                         color=mapal, labels_row = rep("",287),
                         right_annotation = ha,angle_col = "45",
                         annotation_col = annotation_col, annotation_colors = ann_colors,
                         fontsize_row = 6, border_color = NA)

