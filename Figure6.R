library(Seurat) 
library(SeuratObject)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(hrbrthemes)
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
library(cowplot)
library(scales)
library(VennDiagram)
library(AUCell)
library(ggalt)
library(factoextra)
library(FactoMineR)
library(SCopeLoomR) #devtools::install_github("aertslab/SCopeLoomR")
library(AUCell)
library(SCENIC) #devtools::install_github("aertslab/SCENIC")
library(ggalluvial)

species_col <- c("steelblue3","palegreen3","orangered1")
my_color <- colorRampPalette(c("#0da9ce", "white", "#e74a32"))(100)

load(file="./human.seurat.rda")  #human seurat
load(file="./mouse.seurat.rda")  #mouse seurat
load(file="./fly.seurat.rda")

#SCENIC
hm_SCENIC <- open_loom("./hm_SCENIC.loom")

hm_regulons <- get_regulons(hm_SCENIC, column.attr.name="Regulons")
hm_regulons <- regulonsToGeneLists(hm_regulons)
class(hm_regulons)

hm_regulonAUC <- get_regulons_AUC(hm_SCENIC, column.attr.name='RegulonsAUC')

cellinfo <- human_germ@meta.data[,c('Cluster',"nFeature_RNA","nCount_RNA")]
colnames(cellinfo)=c('Cluster','nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'Cluster'))

selectedResolution <- "Cluster"
hm_sub_regulonAUC <- hm_regulonAUC

hm_rss <- calcRSS(AUC=getAUC(hm_sub_regulonAUC),
                  cellAnnotation=cellTypes[colnames(hm_sub_regulonAUC),
                                           selectedResolution])

hm_rss=na.omit(hm_rss)

#mouse
ms_SCENIC <- open_loom("./ms_SCENIC.loom")
ms_regulons <- get_regulons(ms_SCENIC, column.attr.name="Regulons")
ms_regulons <- regulonsToGeneLists(ms_regulons)
class(ms_regulons)

ms_regulonAUC <- get_regulons_AUC(ms_SCENIC, column.attr.name='RegulonsAUC')
ms_regulonAucThresholds <- get_regulon_thresholds(ms_SCENIC)
cellinfo <- mouse_germ@meta.data[,c('Cluster',"nFeature_RNA","nCount_RNA")]
colnames(cellinfo)=c('Cluster','nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'Cluster'))

selectedResolution <- "Cluster"
ms_sub_regulonAUC <- ms_regulonAUC

ms_rss <- calcRSS(AUC=getAUC(ms_sub_regulonAUC),
                  cellAnnotation=cellTypes[colnames(ms_sub_regulonAUC),
                                           selectedResolution])

ms_rss=na.omit(ms_rss)

#fly
dr_SCENIC <- open_loom("./dr_SCENIC.loom")
dr_regulons <- get_regulons(dr_SCENIC, column.attr.name="Regulons")
dr_regulons <- regulonsToGeneLists(dr_regulons)
class(dr_regulons)

dr_regulonAUC <- get_regulons_AUC(dr_SCENIC, column.attr.name='RegulonsAUC')
dr_regulonAucThresholds <- get_regulon_thresholds(dr_SCENIC)

cellinfo <- fly_germ@meta.data[,c('Cluster',"nFeature_RNA","nCount_RNA")]
colnames(cellinfo)=c('Cluster','nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'Cluster'))

selectedResolution <- "Cluster"
dr_sub_regulonAUC <- dr_regulonAUC

dr_rss <- calcRSS(AUC=getAUC(dr_sub_regulonAUC),
                  cellAnnotation=cellTypes[colnames(dr_sub_regulonAUC),
                                           selectedResolution])

dr_rss=na.omit(dr_rss)

#sankey
sankey.df <- read.csv("./Sankey_18TF.csv",header = T)
ggplot(sankey.df, aes(axis1 = Human, axis2 = Mouse, axis3 = Fly)) +
  geom_flow(aes(fill=Human), width = 0.4, curve_type = "cubic", alpha = 0.5) +
  geom_stratum(width = 0.4) +
  geom_text(stat = "stratum", infer.label = TRUE, size = 5)+
  #scale_fill_brewer(palette ="Set3")+
  theme(axis.text.y = element_blank(),legend.position = "none",panel.background = element_blank(),line = element_blank())+
  scale_x_discrete(NULL)

#regulon
hm_rss.hmdf <- data.frame(hm_rss)[,c(1,3,5,2,4)]
ms_rss.hmdf <- data.frame(ms_rss)[,c(5,2,3,1,4)]
dr_rss.hmdf <- data.frame(dr_rss)[,c(2,3,1,4,5)]

hm_hmap <- pheatmap(subset(hm_rss.hmdf, rownames(hm_rss.hmdf) %in% hm.rss.df$Topic),
                    fontsize=12, 
                    fontsize_row = 10, 
                    main = "Human", 
                    scale = "row",
                    angle_col = 45,
                    color = my_color,
                    treeheight_col = 0,  
                    border_color = NA,
                    cluster_cols = F,
                    cluster_rows = F)

ms_hmap <- pheatmap(subset(ms_rss.hmdf, rownames(ms_rss.hmdf) %in% ms.rss.df$Topic),
                    fontsize=12, 
                    fontsize_row = 10, 
                    main = "Mouse",
                    scale = "row", 
                    angle_col = 45,
                    color = my_color,
                    treeheight_col = 0,  
                    border_color = NA,
                    cluster_cols = F,
                    cluster_rows = F)

dr_hmap <- pheatmap(subset(dr_rss.hmdf, rownames(dr_rss.hmdf) %in% dr.rss.df$Topic),
                    fontsize=12, 
                    fontsize_row = 10, 
                    main = "Fly",
                    scale = "row",
                    angle_col = 45,
                    color = my_color,
                    treeheight_col = 0,  
                    border_color = NA,
                    cluster_cols = F,
                    cluster_rows = F)

#YY1
gene.toplot <- as.vector(c("YY1","YY1..."))
plist<-list()
for (i in 1:length(gene.toplot)) {
  p = FeaturePlot(human_germ, features = gene.toplot[i], pt.size = 0.1, reduction = "umap", cols=c("gray85","#6A51A3","#3F007D")) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.title = element_text(size = 25,face = "bold.italic")) + NoLegend()
  plist[[i]] = p
}
p1 <- plot_grid(plist[[1]],plist[[2]], ncol=1)

gene.toplot <- as.vector(c("Yy1","YY1..."))
plist<-list()
for (i in 1:length(gene.toplot)) {
  p = FeaturePlot(mouse_germ, features = gene.toplot[i], pt.size = 0.1, reduction = "umap", cols=c("gray85","#6A51A3","#3F007D")) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.title = element_text(size = 25,face = "bold.italic")) + NoLegend()
  plist[[i]] = p
}
p2 <- plot_grid(plist[[1]],plist[[2]], ncol=1)

gene.toplot <- as.vector(c("pho","YY1..."))
plist<-list()
for (i in 1:length(gene.toplot)) {
  p = FeaturePlot(fly_germ, features = gene.toplot[i], pt.size = 0.1, reduction = "umap", cols=c("gray85","#6A51A3","#3F007D")) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.title = element_text(size = 25,face = "bold.italic")) + NoLegend()
  plist[[i]] = p
}
p3 <- plot_grid(plist[[1]],plist[[2]], ncol=1)


