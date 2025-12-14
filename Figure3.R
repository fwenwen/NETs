library(ggthemes)
library(dplyr)
library(Seurat)
library(ggforce)
library(monocle)

mycds <- readRDS("./cds.rds")
top.markers <- read.delim("./Ncluster.markers.txt")
sig_diff.genes <- subset(top.markers,p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)$gene
sig_diff.genes <- unique(sig_diff.genes)

diff_test <- differentialGeneTest(mycds[sig_diff.genes,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 0.01))
length(sig_gene_names)
plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=5,
                             show_rownames=T, return_heatmap=T)

# BEAM
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
mycds_sub <- mycds[disp.genes,]
plot_cell_trajectory(mycds_sub, color_by = "State")
beam_res <- BEAM(mycds_sub, branch_point = 1, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
# mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 1e-2)),]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
beam_plot <- plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, 
                                         branch_colors = c("#979797", "#F05662", "#7990C8"),
                                         num_clusters = 3, show_rownames = T)
saveRDS(mycds_sub_beam,"mycds_sub_beam.rds")
saveRDS(beam_plot,"beam_plot.rds")

mapal <- colorRampPalette(c(rev(brewer.blues(8))[-(6:8)],brewer.ylorrd(5)))(200)

pdf("beam_plot.pdf",height = 8.27,width = 8.27)
heatmap_byPsedotime <- plot_genes_branched_heatmap(mycds_sub_beam,
                                                   branch_point = 1,
                                                   num_clusters = 4, 
                                                   cores = 1,
                                                   branch_labels = c("Cell fate 1", "Cell fate 2"),
                                                   # hmcols = colorRampPalette(rev(c("#CF384D","#ED6345","#FA9A58","#F3FAAD","#D1EC9C","#96D5A4","#5BB6A9","#3682BA")))(62),# colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                                   # hmcols = jet(100),
                                                   hmcols = colorRampPalette(c("#211D1E","#174263","#0E66A4","#5D9FCF","#C8DFED","#FCE6CE","#F6B66C","#ED8522","#D75728","#C42630"))(62),
                                                   branch_colors = c("lightgrey", '#B07AA1', '#76B7B2'), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                                   use_gene_short_name = TRUE,
                                                   show_rownames = TRUE,
                                                   return_heatmap = T)
dev.off()

saveRDS(heatmap_byPsedotime,"heatmap_byPsedotime.rds")
meta <- data.frame(Cluster=heatmap_byPsedotime$annotation_row$Cluster,
                   marker=heatmap_byPsedotime$ph$tree_row$labels)
head(meta)
geneCluster <- setNames(as.numeric(as.character(meta$Cluster)),meta$marker)
gstat <- geneCluster[heatmap_byPsedotime$ph$tree_row$labels[heatmap_byPsedotime$ph$tree_row$order]]
pst.genes <- droplevels(top.markers %>% filter(cluster %in% pData(mycds)$cluster)) %>%
  arrange(desc(avg_log2FC)) %>% filter(!duplicated(gene))
dim(pst.genes)

k <- 4
plist <- list()
for(i in 1:k){
  pd <- as.data.frame(table(pst.genes$cluster[pst.genes$gene %in% names(gstat)[gstat==i]]))
  plist[[i]] <- ggpie(pd, "Freq", label = "Var1",
                      fill = "Var1", color = "white",
                      palette = ncols) + ggtitle(i) + NoLegend()
}
n18 <- ggarrange(plotlist = plist, ncol = 4, nrow = 1)
pdf("gene.pie.beam.pdf",height = 3,width = 10)
print(n18)
dev.off()
# 
phtd <- data.frame(x=as.numeric(gstat), y=as.numeric(gstat), row.names = names(gstat))
pst.genes <- as.data.frame(pst.genes) 
rownames(pst.genes) <- as.character(pst.genes$gene)
pst.top <- pst.genes %>% dplyr::group_by(cluster) %>%
  filter(avg_log2FC > 0.5) %>% 
  dplyr::top_n(n=20, wt=avg_log2FC)
pst.top <- as.data.frame(sort(na.omit(gstat[pst.top$gene])))
head(pst.top)
pst.top$name <- rownames(pst.top)
for(i in 1:k){
  print(paste(">>>>", i))
  print(paste(pst.top[pst.top[,1]==i, 'name'], collapse = ", "))
}
annotation_row <- data.frame(cluster=ifelse(rownames(phtd) %in% rownames(pst.genes), 
                                            as.character(pst.genes[rownames(phtd),'cluster']), "Neus0"), 
                             row.names=rownames(phtd))
head(annotation_row)
ann_colors <- list(cluster=c("Neus0"=NA, ncols))
n22 <- pheatmap::pheatmap(phtd, cutree_rows = k, cluster_rows = F,
                          labels_row =pst.top$name,
                          annotation_row = annotation_row, 
                          annotation_colors = ann_colors)

pdf("./k_heatmap.pdf", height=8.27, width=8.27)
print(n22)
dev.off()
# 
enrichOut <- lapply(1:k, function(i){
  geneSet <- names(geneCluster)[which(geneCluster==i)]
  enrichGO(OrgDb=org.Hs.eg.db,
           gene = geneSet,
           pvalueCutoff = 0.05,
           qvalueCutoff = 1,
           keyType = 'SYMBOL',
           pAdjustMethod = 'fdr',
           ont = "BP")
})
names(enrichOut) <- sprintf("K%s",1:k)
enrichOut <- lapply(1:k, function(i){
  data.frame(Cluster=sprintf("K%s",i), enrichOut[[i]])
})
enrichRes <- do.call("rbind", enrichOut) 
write.table(enrichRes, file="neu.pseudotime.GO.beam.txt", quote = F, row.names = F, sep="\t")
# 

pData(mycds)$Pseudotime2 <- rank(pData(mycds)$Pseudotime, ties.method = "first")
stateOrder <- pData(mycds) %>% arrange(Pseudotime) %>% 
  dplyr::select(State) %>% unique() %>% unlist() %>% as.character()
stateOrder <- c("2", "3",  "5", "4","1")
stateOrder <- setNames(sprintf("S%01d",1:length(stateOrder)), stateOrder)
pData(mycds)$State2 <- factor(stateOrder[as.character(pData(mycds)$State)])
stateColor <- setNames(kovesi.rainbow_bgyr_35_85_c73(length(stateOrder)), stateOrder)

cs1 <- table(pData(mycds)$cluster, pData(mycds)$State2)
cs2 <- table(pData(mycds)$group, pData(mycds)$State2)
cs1 <- as.data.frame(cs1 * 100 / rowSums(cs1)) 
cs2 <- as.data.frame(cs2 * 100 / rowSums(cs2))

n14 <- ggdensity(pData(mycds), x = "Pseudotime2", alpha=0.1, scales='free_y',
                 add = "mean", rug = F, facet.by="cluster", 
                 color = "group", fill = "group", palette = gcols)
n15 <- ggdensity(pData(mycds), x = "Pseudotime2", alpha=0.1, scales='free_y',
                 add = "mean", rug = F, facet.by="group", 
                 color = "cluster", fill = "cluster", 
                 palette = ncols)
n16 <- ggdensity(pData(mycds), x = "Pseudotime2", alpha=0.1, rug = F, 
                 color = "group", fill = "group", palette = gcols)
n17 <- ggdensity(pData(mycds), x = "Pseudotime2", alpha=0.1, rug = F, 
                 color = "cluster", fill = "cluster", palette = ncols)

n71 <- cs1 %>% mutate(State=factor(Var2, levels=rev(stateOrder))) %>% 
  ggplot(aes(x=Var1, y=Freq, fill=State)) + 
  geom_bar(stat="identity") + rotate() + xlab("") +
  scale_fill_manual(values=stateColor) + theme_pubr()
n72 <- cs1 %>% mutate(State=Var2) %>% 
  ggplot(aes(x=Var1, y=Freq, fill=State)) + 
  geom_bar(stat="identity") + rotate() + xlab("") +
  scale_fill_manual(values=stateColor) + theme_pubr()
n81 <- cs2 %>% mutate(State=factor(Var2, levels=rev(stateOrder))) %>% 
  ggplot(aes(x=Var1, y=Freq, fill=State)) + 
  geom_bar(stat="identity") + rotate() + xlab("") +
  scale_fill_manual(values=stateColor) + theme_pubr()
n82 <- cs2 %>% mutate(State=Var2) %>% 
  ggplot(aes(x=Var1, y=Freq, fill=State)) + 
  geom_bar(stat="identity") + rotate() + xlab("") +
  scale_fill_manual(values=stateColor) + theme_pubr()

n9 <- monocle::plot_cell_trajectory(mycds, cell_size=1, color_by="State2", show_branch_points=F) + 
  scale_color_manual(values=stateColor)

n10 <- monocle::plot_cell_trajectory(mycds, cell_size=1, color_by="Pseudotime", show_branch_points=F) +
  gradient_color(kovesi.rainbow(100))
n11 <- monocle::plot_cell_trajectory(mycds, cell_size=1, color_by="cluster", show_branch_points=F) + 
  scale_color_manual(values=ncols)
n12 <- monocle::plot_cell_trajectory(mycds, cell_size=1, color_by="group", show_branch_points=F) + 
  scale_color_manual(values=gcols)

phtd <- data.frame(x=as.numeric(gstat), y=as.numeric(gstat), row.names = names(gstat))
pst.genes <- as.data.frame(pst.genes) 
rownames(pst.genes) <- as.character(pst.genes$gene)
pst.top <- pst.genes %>% dplyr::group_by(cluster) %>%
  filter(avg_log2FC > 0.5) %>% 
  dplyr::top_n(n=10, wt=avg_log2FC)
pst.top <- as.data.frame(sort(na.omit(gstat[pst.top$gene])))
head(pst.top)
pst.top$name <- rownames(pst.top)
for(i in 1:k){
  print(paste(">>>>", i))
  print(paste(pst.top[pst.top[,1]==i, 'name'], collapse = ", "))
}
annotation_row <- data.frame(cluster=ifelse(rownames(phtd) %in% rownames(pst.genes), 
                                            as.character(pst.genes[rownames(phtd),'cluster']), "Neus0"), 
                             row.names=rownames(phtd))
head(annotation_row)
ann_colors <- list(cluster=c("Neus0"=NA, ncols))
n22 <- pheatmap::pheatmap(phtd, cutree_rows = k, labels_row =pst.top$name,
                          annotation_row = annotation_row, 
                          annotation_colors = ann_colors)

pdf("./k_heatmap.pdf", height=8.27, width=8.27)
print(n22)
dev.off()

pdf("./Neu.pseudotime.pdf", height=8.27, width=8.27)
pal.bands(kovesi.rainbow_bgyr_35_85_c73, main="Kovesi") 
nl <- ggarrange(as_ggplot(get_legend(n72)), 
                as_ggplot(get_legend(n9)), as_ggplot(get_legend(n10)),
                as_ggplot(get_legend(n11)), as_ggplot(get_legend(n12)),
                nrow = 3, ncol = 2)
print(nl)
print(n14)
print(n15)
print(ggarrange(n71 + NoLegend(), n81 + NoLegend(), nrow=2, ncol=1))
print(ggarrange(n9 + NoAxes() + NoLegend(), 
                n10 + NoAxes() + NoLegend(), 
                n12 + NoAxes() + NoLegend(), 
                n11 + NoAxes() + NoLegend(), 
                nrow=2, ncol=2))

print(ggarrange(n10 + NoAxes(), 
  n11 + NoAxes() , 
  nrow=1, ncol=2))

plot.new()
print(pseudotime_out)
print(n18)
dev.off()

# 
plot_expression_trajectory <- function(gene){
  pData(mycds)$GeneScore <- GetAssayData(NeuObj, assay = 'RNA', slot="data")[gene,rownames(pData(mycds))]
  ## pData(mycds)$GeneScore <- ifelse(pData(mycds)$GeneScore>0, pData(mycds)$GeneScore, NA)
  no <- monocle::plot_cell_trajectory(mycds, cell_size=1, color_by="GeneScore", show_branch_points=F) +
    gradient_color(mapal)
  no
}
n13 <- plot_expression_trajectory('S100A12')


n23 <- plot_expression_trajectory('CD14')
n24 <- plot_expression_trajectory('TGM2')+NoAxes()+NoLegend()#+ggtitle("TGM2")
n25 <- plot_expression_trajectory('LST1')
n26 <- plot_expression_trajectory('FCGR3B')
n27 <- plot_expression_trajectory('ESR1')+NoAxes()+NoLegend()#+ggtitle("ESR1")

n24+n27

print(ggarrange(n11 + NoAxes() + NoLegend(), 
                n23 + NoAxes() + NoLegend(), 
                n24 + NoAxes() + NoLegend(), 
                n25 + NoAxes() + NoLegend(), 
                n26 + NoAxes() + NoLegend(), 
                n27 + NoAxes() + NoLegend(), 
                nrow=2, ncol=3))

print(ggarrange(n13 + NoAxes() + NoLegend(), 
                n23 + NoAxes() + NoLegend(), 
                n24 + NoAxes() + NoLegend(), 
                n25 + NoAxes() + NoLegend(), 
                n26 + NoAxes() + NoLegend(), 
                n27 + NoAxes() + NoLegend(), 
                nrow=2, ncol=3))
