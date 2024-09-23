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

species_col <- c("steelblue3","palegreen3","orangered1")

#load markers
load(file="./human_AllMarkers.rda")
load(file="./mouse2_AllMarkers.rda")
load(file="./fly_AllMarkers.rda")


human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
droso <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse_to_human_spermatogenesis.list<-getLDS(values = mouse_spermatogenesis.markers,mart = mouse,attributes = "mgi_symbol",filters = "mgi_symbol", martL=human,attributesL = "hgnc_symbol") 
fly_to_human_spermatogenesis.list<-getLDS(values = elife_fly_spermatogenesis.markers,mart = droso,attributes = "external_gene_name",filters = "external_gene_name", martL=human,attributesL = "hgnc_symbol") 

venn.diagram(
  x = list(human_spermatogenesis.markers, mouse_spermatogenesis.list$HGNC.symbol, fly_to_human_spermatogenesis.list$HGNC.symbol),
  category.names = c("Human" , "Mouse" , "Fly"),
  filename = NULL,
  fill = species_col,
  cex = 2.5, # number size
  cat.cex = 2.5, # set size
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.dist = c(0.06, 0.06, 0.05),
  print.mode = c("raw","percent")
)


#load germ cell seurat object
load(file="./hu_seurat_annot.rda")  
load(file="./ms_seurat_annot.rda")
load(file="./dm_seurat_annot.rda")

#germ cell marker
#human
gene.toplot <- c("UTF1","HORMAD1","ACRV1","PRM3")
gene.toplot <- as.vector(gene.toplot)
plist<-list()

for (i in 1:length(gene.toplot)) {
  p = FeaturePlot(human_germ, features = gene.toplot[i], pt.size = 0.1, reduction = "umap", cols=c("gray85","#6A51A3","#3F007D")) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),plot.title = element_text(face = "bold.italic")) + NoLegend()
  plist[[i]] = p
}
plot_grid(plist[[1]],plist[[2]],plist[[3]], plist[[4]],ncol=1) 

#mouse
gene.toplot <- c("Uchl1","Sycp3","Acrv1","Prm3")
gene.toplot <- as.vector(gene.toplot)
plist<-list()

for (i in 1:length(gene.toplot)) {
  p = FeaturePlot(mouse_germ, features = gene.toplot[i], pt.size = 0.1, reduction = "umap", cols=c("gray85","#6A51A3","#3F007D")) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),plot.title = element_text(face = "bold.italic")) + NoLegend()
  plist[[i]] = p
}
plot_grid(plist[[1]],plist[[2]],plist[[3]], plist[[4]],ncol=1)

#fly
gene.toplot <- c("bam","fzo","twe","p-cup")
gene.toplot <- as.vector(gene.toplot)
plist<-list()

for (i in 1:length(gene.toplot)) {
  p = FeaturePlot(fly_germ, features = gene.toplot[i], pt.size = 0.1, reduction = "umap", cols=c("gray85","#6A51A3","#3F007D")) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),plot.title = element_text(face = "bold.italic")) + NoLegend()
  plist[[i]] = p
}
plot_grid(plist[[1]],plist[[2]],plist[[3]], plist[[4]],ncol=1)

#monocle
#load cds
load(file = "./human.monocle_2.rda")
load(file = "./mouse.monocle.rda")
load(file = "./fly.monocle.rda")

#prepare dataframe
human.pseudo <- data.frame(FetchData(object = human_germ, vars = c("Cluster")))
mouse.pseudo <- data.frame(FetchData(object = mouse_germ, vars = c( "Cluster")))
fly.pseudo <- data.frame(FetchData(object = fly_germ, vars = c( "Cluster")))

human.pseudo$Pseudotime <-ceiling( rank(hm_cds$Pseudotime)/(max(rank(hm_cds$Pseudotime)))*50)
mouse.pseudo$Pseudotime <-ceiling( rank(ms_cds$Pseudotime)/(max(rank(ms_cds$Pseudotime)))*50)
fly.pseudo$Pseudotime <-ceiling( rank(dr_cds$Pseudotime)/(max(rank(dr_cds$Pseudotime)))*50)
hm_cds$Germ <- Idents(human_germ)
ms_cds$Germ <- Idents(mouse_germ)
dr_cds$Germ <- Idents(fly_germ)

plot_cell_trajectory(hm_cds, color_by = "Germ", cell_size = 1, show_branch_points = F) + scale_color_manual(values = cell_col[1:3]) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
plot_cell_trajectory(ms_cds, color_by = "Germ", cell_size = 1, show_branch_points = F) + scale_color_manual(values = cell_col[1:3]) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()
plot_cell_trajectory(dr_cds, color_by = "Germ", cell_size = 1, show_branch_points = F) + scale_color_manual(values = cell_col[1:3]) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) + NoLegend()

#ridge plot
ggplot(human.pseudo)+
  geom_density(aes(x=Pseudotime, fill = Cluster),
               alpha = 0.9, position = "fill")+
  theme_test()+
  scale_fill_manual(values = cell_col[1:3]) + NoLegend()
ggplot(mouse.pseudo)+
  geom_density(aes(x=Pseudotime, fill = Cluster),
               alpha = 0.9, position = "fill")+
  theme_test()+
  scale_fill_manual(values = cell_col[1:3]) + NoLegend()
ggplot(fly.pseudo)+
  geom_density(aes(x=Pseudotime, fill = Cluster),
               alpha = 0.9, position = "fill")+
  theme_test()+
  scale_fill_manual(values = cell_col[1:3]) + NoLegend()











