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

um_cols <- c("#E01516","#4552A0","darkgoldenrod2","darkgrey")
species_col <- c("steelblue3","palegreen3","orangered1")

#load seurat object
load(file="./human_all_seurat.rda")  #all testis cell seurat
load(file="./mouse_all_seurat.rda")
load(file="./fly_all_seurat.rda")

Idents(human_all_testis) <- factor(Idents(human_all_testis),levels = c("Spermatogonia","Spermatocyte","Spermatid","Somatic cells"))
Idents(mouse_all_testis) <- factor(Idents(mouse_all_testis),levels = c("Spermatogonia","Spermatocyte","Spermatid","Somatic cells"))
Idents(fly_all_testis) <- factor(Idents(fly_all_testis),levels = c("Spermatogonia","Spermatocyte","Spermatid","Somatic cells"))

#umap
DimPlot(human_all_testis, reduction = "umap", label=F, cols = um_cols, pt.size = 0.2)+theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +NoLegend()
DimPlot(mouse_all_testis, reduction = "umap", label=F, cols = um_cols, pt.size = 0.2)+theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +NoLegend()
DimPlot(fly_all_testis, reduction = "umap", label=F, cols = um_cols, pt.size = 0.2)+theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +NoLegend()

#marker featureplot
#human
gene.toplot <- c("UCHL1","SYCP3","ACRV1","PRM3","PDGFRA","PECAM1","CD14")
gene.toplot <- as.vector(gene.toplot)
plist<-list()

for (i in 1:length(gene.toplot)) {
  p = FeaturePlot(human_all_testis, features = gene.toplot[i], pt.size = 0.1, reduction = "umap", cols=c("gray85","#6A51A3","#3F007D")) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.title = element_text(size = 25,face = "bold.italic")) + NoLegend()
  plist[[i]] = p
}
plot_grid(plist[[1]],plist[[2]],plist[[3]], plist[[4]],plist[[5]],plist[[6]],plist[[7]],ncol=7) #dim=1000:1080

#mouse
gene.toplot <- c("Uchl1","Sycp3","Acrv1","Prm3","Sox9","Insl3","Tcf21")
gene.toplot <- as.vector(gene.toplot)
plist<-list()

for (i in 1:length(gene.toplot)) {
  p = FeaturePlot(mouse_all_testis, features = gene.toplot[i], pt.size = 0.1, reduction = "umap", cols=c("gray85","#6A51A3","#3F007D")) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.title = element_text(size = 25,face = "bold.italic")) + NoLegend()
  plist[[i]] = p
}
plot_grid(plist[[1]],plist[[2]],plist[[3]], plist[[4]],plist[[5]],plist[[6]],plist[[7]],ncol=7) #dim=1000:1080

#fly
gene.toplot <- c("aub","fzo","twe","p-cup","MtnA","Hsp23","dlg1")
gene.toplot <- as.vector(gene.toplot)
plist<-list()

for (i in 1:length(gene.toplot)) {
  p = FeaturePlot(fly_all_testis, features = gene.toplot[i], pt.size = 0.1, reduction = "umap", cols=c("gray85","#6A51A3","#3F007D")) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.title = element_text(size = 25,face = "bold.italic")) + NoLegend()
  plist[[i]] = p
}
plot_grid(plist[[1]],plist[[2]],plist[[3]], plist[[4]],plist[[5]],plist[[6]],plist[[7]],ncol=7) #dim=1000:1080

#PCA
#transfer to ortholog
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
droso <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

mouse_to_human.list<-getLDS(values = rownames(mouse_all_testis),mart = mouse, attributes = c("mgi_symbol"),filters = "mgi_symbol", martL=human, attributesL = c("hgnc_symbol")) 
mouse_to_human.type<-getLDS(values = rownames(mouse_all_testis),mart = mouse, attributes = c("mgi_symbol"),filters = "mgi_symbol", martL=human, attributesL = c("mmusculus_homolog_orthology_type")) 
fly_to_human.list<-getLDS(values = rownames(fly_all_testis@assays$RNA),mart = droso,attributes = "external_gene_name",filters = "external_gene_name", martL=human, attributesL = c("hgnc_symbol")) 
fly_to_human.type<-getLDS(values = rownames(fly_all_testis@assays$RNA),mart = droso,attributes = "external_gene_name",filters = "external_gene_name", martL=human, attributesL = c("dmelanogaster_homolog_orthology_type")) 

mouse_to_human.df <- subset(merge(mouse_to_human.list,mouse_to_human.type), Mouse.homology.type == "ortholog_one2one")
fly_to_human.df <- subset(merge(fly_to_human.list,fly_to_human.type), Drosophila.melanogaster.homology.type == "ortholog_one2one")

cons.one2one <- intersect(rownames(human_all_testis),intersect(mouse_to_human.df$HGNC.symbol,fly_to_human.df$HGNC.symbol))

one2one.df<-getBM(values = cons.one2one, 
                  mart = human,
                  filters = "external_gene_name", 
                  attributes = c("external_gene_name", 
                                 "mmusculus_homolog_associated_gene_name",
                                 "dmelanogaster_homolog_associated_gene_name",
                                 "mmusculus_homolog_orthology_type",
                                 "dmelanogaster_homolog_orthology_type")) 

one2one.df <- subset(one2one.df, mmusculus_homolog_orthology_type == "ortholog_one2one")
one2one.df <- subset(one2one.df, dmelanogaster_homolog_associated_gene_name %in% rownames(fly_all_testis@assays$integrated)) 

human_all_testis$samplecelltype <- paste("Human",Idents(human_all_testis),sep = "-")
mouse_all_testis$samplecelltype <- paste("Mouse",Idents(mouse_all_testis),sep = "-")
fly_all_testis$samplecelltype <- paste("Fly",Idents(fly_all_testis),sep = "-")

#load PCA dataframe
load(file="./pca.input.rda")
pca <- PCA(t(PCA.input))
fviz_pca_var(pca,geom = c("point"))
data=pca$ind$coord
data1=data.frame(data[,c(1,2)])
data1$sample=rownames(data1)
data1$Species=stringr::str_split_fixed(data1$sample,"\\.",n=3)[,2]
data1$Cell=stringr::str_split_fixed(data1$sample,"\\.",n=3)[,3]

data1$Cell <- factor(data1$Cell,levels = c("Spermatogonia", "Spermatocyte", "Spermatid","Somatic cells"))
data1$Species <- factor(data1$Species,levels = c("Human", "Mouse", "Fly"))

ggplot(data1) +
  geom_point(aes(Dim.1,Dim.2,
                 color= Species,
                 shape= Cell), size=7, stroke = 2.5) +
  scale_shape_manual(values=c(1,4,7,11)) +
  scale_color_manual(values = species_col) +
  stat_ellipse(level = 0.8, size = 1, show.legend = F,
               aes(Dim.1,Dim.2, color=Species), 
               geom ="path") +
  theme_test() +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14), legend.position = "right",
        panel.border = element_rect(color = "black",size = 1)) +
  guides(color = guide_legend(override.aes = list(size = 3)),shape = guide_legend(override.aes = list(size = 4, stroke = 1)))






















