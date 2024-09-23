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


load(file = "./sperm_fail.gene.rda")
load(file = "./human.monocle_2.rda")
pt.matrix <- exprs(hm_cds)[match(sperm_fail.gene,rownames(hm_cds@featureData)),order(hm_cds$Pseudotime)] 
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- sperm_fail.gene

htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11), rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 5, fontface = 'italic'),
  km = 4,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

#fly
load(file = "./candidates_list.rda")
load(file = "./fly.monocle.rda")
pt.matrix <- exprs(dr_cds)[match(candidates_list,rownames(dr_cds@featureData)),order(dr_cds$Pseudotime)] 
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- candidates_list

#draw gene expression according pseudotime
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11), rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 5, fontface = 'italic'),
  km = 3,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
htkm








