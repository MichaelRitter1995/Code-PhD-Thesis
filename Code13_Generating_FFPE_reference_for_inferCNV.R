setwd("/omics/groups/OE0146/internal/Micha")
library(Seurat)
library(harmony)
library(ggplot2)
library(SCpubr)
library(ggplotify)

#Load the relaxed and strict merged object, which is including all spotrs aof all samples used as reference
relaxed <- readRDS("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/yiheng/p16/preproc/seurat_relaxed.rds")
strict <- readRDS("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/yiheng/p16/preproc/seurat_strict.rds")

p1 <- DimPlot(strict, reduction = "umap", group.by = "orig.ident")
p1
p2 <- DimPlot(relaxed, reduction = "umap", group.by = "orig.ident")
p2


#Select only healthy spots
strict <- subset(x = strict, subset = reference == TRUE)
relaxed <- subset(x = relaxed, subset = reference == TRUE)


#merge both objects

merged <- merge(strict, y = relaxed, add.cell.ids = c("strict", "relaxed"), project = "Reference_set")
DefaultAssay(merged) <- "SCT"

#harmony

Features <- SelectIntegrationFeatures(list(relaxed, strict))
VariableFeatures(merged) <- Features
#merged <- FindVariableFeatures(object = merged)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
harm <- RunHarmony(merged, group.by.vars = "orig.ident")
harm <- RunUMAP(harm, reduction = "harmony", dims = 1:15)
harm <- FindNeighbors(harm, reduction = "harmony", dims = 1:15)
harm <- FindClusters(harm, resolution = 0.5)

p2 <- DimPlot(harm, group.by = "orig.ident")
p2


###Add Cell Type Markers###

Cell_Types<- list(
  c("MBP","OLIG2","OLIG3","CLDN11","MOG","SOX10","GJB1","MAG","CNP","PLP1"),
  c("ENO2","NEFL","NEFM","NEFH","INA","SLC17A7","RBFOX3","SNAP25","SYT1","NEUROD2"),
  c("NCAN", "NTSR2", "ADGRV1", "SLC1A2", "SLCO1C1", "ACSBG1", "GFAP", "MLC1", "PIRT", "SLC7A10"),
  c("GABRA1", "NEFH", "INA", "NEFM", "STMN2", "VSNL1", "CALY", "SV2C", "RGS16", "UCHL1", "THY1", "NGS1", "RAB3C")
)
names(Cell_Types) <- c("Oligodendrocytes","Neurons", "Astrocytes", "strange_cluster")
harm <- AddModuleScore(object =harm, features= Cell_Types, name = c("Oligodendrocytes", "Neurons", "Astrocytes"), assay ='Spatial')

###Plot Results###
SCpubr::do_DimPlot(harm)
SCpubr::do_DimPlot(harm, group.by = "orig.ident")
SCpubr::do_FeaturePlot(harm, features = "nCount_Spatial", max.cutoff = 10000)+ggplot2::scale_color_viridis_c(option = "turbo")
SCpubr::do_FeaturePlot(harm, features = "Oligodendrocytes1", max.cutoff = 40)+ggplot2::scale_color_viridis_c(option = "turbo")
SCpubr::do_FeaturePlot(harm, features = "Neurons2", min.cutoff = -10, max.cutoff = 10)+ggplot2::scale_color_viridis_c(option = "turbo")

SCpubr::do_DimPlot(harm, reduction = "pca", group.by = "orig.ident")
SCpubr::do_FeaturePlot(harm, reduction = "pca", features = "Oligodendrocytes1", max.cutoff = 40)+ggplot2::scale_color_viridis_c(option = "turbo")
SCpubr::do_FeaturePlot(harm,  reduction = "pca", features = "Neurons2", min.cutoff = -10, max.cutoff = 10)+ggplot2::scale_color_viridis_c(option = "turbo")

pbmc.markers <- FindAllMarkers(harm, only.pos = TRUE, recorrect_umi=FALSE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster4 <- pbmc.markers[which(pbmc.markers$cluster==4),]
cluster4 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster3 <- pbmc.markers[which(pbmc.markers$cluster==3),]
cluster3 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
