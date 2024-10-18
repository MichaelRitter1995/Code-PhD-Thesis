library(SPATA2)
library(stringr)
library(ggplot2)
library(dplyr)
library(Seurat)

setwd("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample")

sample_name <- "GB5"
x <- "GB5"
pathSeu <- ("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/visium_Analysis/GBM/per_sample/")
seurat_obj <- readRDS(paste0(pathSeu, x, "/",x,"_seurat_obj.rds"))

IHC <- c("#FFFFFF", "#DDD4CF","#E1C07E", "#9A401C", "#6F2C36", "#3B1112")

DefaultAssay(seurat_obj) <- "SCT"
SpatialFeaturePlot(seurat_obj, features = "MKI67")+scale_fill_gradientn(colours = IHC)+scale_alpha(name = waiver(), range=c(0,1))

SpatialFeaturePlot(seurat_obj, features = "nCount_Spatial")+scale_fill_gradientn(colours = viridis::magma(10), limits=c(0,50000), na.value = "red")


sample_name <- "GB1"
x <- "GB1"
pathSeu <- ("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/visium_Analysis/GBM/per_sample/")
seurat_obj <- readRDS(paste0(pathSeu, x, "/",x,"_seurat_obj.rds"))

IHC <- c("#FFFFFF", "#DDD4CF","#E1C07E", "#9A401C", "#6F2C36", "#3B1112")

DefaultAssay(seurat_obj) <- "SCT"
SpatialFeaturePlot(seurat_obj, features = "GFAP")+scale_fill_gradientn(colours = IHC)+scale_alpha(name = waiver(), range=c(0,1))

sample_name <- "GB10"
x <- "GB10"
pathSeu <- ("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/visium_Analysis/GBM/per_sample/")
seurat_obj <- readRDS(paste0(pathSeu, x, "/",x,"_seurat_obj.rds"))

IHC <- c("#FFFFFF", "#DDD4CF","#E1C07E", "#9A401C", "#6F2C36", "#3B1112")

DefaultAssay(seurat_obj) <- "SCT"
SpatialFeaturePlot(seurat_obj, features = "GFAP")+scale_fill_gradientn(colours = IHC)+scale_alpha(name = waiver(), range=c(0,1))
