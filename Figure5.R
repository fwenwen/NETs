library(dplyr)
library(Seurat)
library(patchwork)
library(pals)
library(pheatmap)
library(ggplot2)
library(ggthemes)

seuratObj <- readRDS("~/mouse/mouseObj.rds")
DimPlot(seuratObj,group.by = "orig.ident")
Cols <- c(`Unciliated epithelial cells` = "#A3BFDB", 
          `Ciliated epithelial cells` = "#AED8A8",
          `Fibroblasts` = "#CEA9D2" ,
          `Smooth muscle cells` = "#FFBF8E" ,
          `Endothelial cells` = "#C99C86",
          `Cycling cells`="#8DD3C7",
          `Dendritic cells` = "#F2E2C9",Macrophages = "#F781B8" , 
          Neutrophils = "#CDCDCD" , 
          `NK cells` = "#FA9193",`T cells` =  "#FEC1DF", 
          `B cells` = "#FFED6F" ) 

DimPlot(seuratObj,group.by = "celltype",cols = Cols)

p1 <- DimPlot(seuratObj,group.by = "celltype",pt.size = 0.2,cols = Cols)
p2 <- DimPlot(seuratObj,group.by = "group",pt.size = 0.2,cols = gcols)
p1+p2
ggarrange(p1+NoAxes()+NoLegend()+ggtitle(""),
          p2+NoAxes()+NoLegend()+ggtitle(""),
          ncol = 2)

pdf("figures/meta.pdf", height=4, width=8)
ggarrange(p1+NoAxes()+NoLegend()+ggtitle(""),
          p2+NoAxes()+NoLegend()+ggtitle(""),
          ncol = 2)
dev.off()


pdf("figures/meta_legend.pdf", height=3, width=5)
ggarrange(as_ggplot(get_legend(p1)),
          as_ggplot(get_legend(p2)),
          ncol = 2)
dev.off()

# 
cellstat1 <- table(seuratObj$group, seuratObj$celltype)
cellstat1 <- (cellstat1 * 100 / rowSums(cellstat1))
cellstat1 <- as.data.frame(cellstat1)
colnames(cellstat1) <- c("group", "cellType", "percentage")
cellstat1$group <- factor(cellstat1$group,
                          levels = c("Normal","EN","EH"))

ggplot(data=cellstat1, aes(x=group, y=percentage, fill=group)) +
  geom_point(size=2, shape=21, color="black") +  
  scale_fill_manual(values=gcols) +
  theme_classic() +
  facet_grid(.~cellType,scales = "free")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

print(ggplot(data=cellstat1[cellstat1$cellType%in% levels(seuratObj$celltype)[7:12],], aes(x=group, y=percentage, color=group, group=group)) +  
        geom_point(size=3) +  
        geom_line(aes(group=cellType), size=0.8, alpha=1) +  
        scale_color_manual(values=gcols) +  
        theme_classic() +  
        facet_wrap(~ cellType, scales = "fixed", nrow = 1) +  
        theme(axis.title.x = element_blank(), axis.text.x = element_blank()))
