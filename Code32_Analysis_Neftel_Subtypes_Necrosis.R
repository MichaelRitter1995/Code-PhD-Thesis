library(SPATA2)
library(stringr)
library(ggplot2)
library(Seurat)
library(stringr)
library(viridis)
library(SCpubr)

setwd("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")
pathVis <- ("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")

#Change sample name to the right sample
sample_name <- "GB12"
x <-  "GB12"

#sample_name <- "GB1"
#x <- "GB1"

#sample_name <- "MNG1"
#x <- "MNG1"
#spata_obj <- updateSpataObject(readRDS(file = paste0("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/bhuvic/Necrotic_Meningioma/SPATA/",sample_name,"_spata.rds")))


#Load spata object
spata_obj <- updateSpataObject(readRDS(paste0(pathVis, x, "/SPATA_obj_", x, ".RDS")))

#Plot Annotation
anno <- read.csv(paste0(x,"/Annotated.csv"), fill=TRUE)
names(anno) <- c("barcodes", "Histology")
anno[anno$Histology=="", ]$Histology <- "Unknown"
spata_obj <- SPATA2::addFeatures(spata_obj, anno, overwrite = T)

nr_of_hist <- anno$Histology %>% unique() %>% length()
colors_histology <- colorRampPalette(RColorBrewer::brewer.pal(13, "Set3"))(nr_of_hist)
names(colors_histology) <- anno$Histology %>% unique()

p1 <- 
  plotSurface(spata_obj, color_by="Histology", pt_alpha=0.5)+
  scale_color_manual(values=colors_histology)+
  ggpLayerThemeCoords() +
  ggpLayerAxesSI(spata_obj, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
  ggtitle("Annotated Histology")
p1
ggsave(filename = paste0(x, "/",i,"_",x,"_Pathologist_Annotation.pdf"),plot=p1, width = 6, height = 5)

#Read matrix containing results of deconvolution of Subtypes
Subtypes <- t(read.csv(paste0(x,"/",x,"_predictions_Neftel.csv"), fill=TRUE, header=TRUE, stringsAsFactors = FALSE))
colnames(Subtypes) <- Subtypes[1,]
Subtypes <- Subtypes[2:nrow(Subtypes),]

#Add Deconvolution scores to data.frame
for(i in 1:ncol(Subtypes)){
    type <- data.frame(cbind(rownames(Subtypes), Subtypes[,i]))
    colnames(type) <- c("barcodes", paste0(colnames(Subtypes)[[i]], "_D_Chr"))
    type[,2] <- as.numeric(type[,2])
    type[,1] <- as.character(type[,1])
    type$barcodes <- gsub("\\.", "-", type$barcodes)
    spata_obj <- SPATA2::addFeatures(spata_obj, type)
    #Plot and save results
    plot <- plotSurface(spata_obj, paste0(colnames(Subtypes)[[i]], "_D_Chr"), smooth = T, pt_alpha=0.75)+
      ggpLayerThemeCoords() +
      ggpLayerAxesSI(spata_obj, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
      ggtitle(paste0(colnames(Subtypes)[[i]], "_D_Chr"))
    ggsave(filename = paste0(x, "/",i,"_",x,"_", colnames(Subtypes)[[i]], "_D_Chr.pdf"),plot=plot, width = 6, height = 5)
    }

#Read List of marker
Marker <- read.csv("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/Reference_sets/GB_Subtypes_Neftel.csv")

for(i in 1:length(unique(Marker$List))){
  Marker_type <- Marker[which(Marker$List ==unique(Marker$List)[[i]]),]
  #Expr <- spata_obj@data[["s_GB12"]][["scaled"]][which(rownames(spata_obj@data[["s_GB12"]][["scaled"]]) %in% Marker_type$Name),]
  Expr <- spata_obj@data[["MNG1"]][["scaled"]][which(rownames(spata_obj@data[["MNG1"]][["scaled"]]) %in% Marker_type$Name),]
  Expr_mean <- colMeans(Expr)
  Expr_mean <- data.frame(cbind(colnames(Expr), Expr_mean))
  Expr_mean[,2] <- as.numeric(Expr_mean[,2])
  colnames(Expr_mean) <- c("barcodes", paste0(unique(Marker$List)[[i]],"_mean"))
  spata_obj <- SPATA2::addFeatures(spata_obj, Expr_mean, overwrite = T)
  plot <- plotSurface(spata_obj, paste0(unique(Marker$List)[[i]],"_mean"), smooth = T, pt_alpha=0.75)+
    ggpLayerThemeCoords() +
    ggpLayerAxesSI(spata_obj, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
    ggtitle(paste0(unique(Marker$List)[[i]],"_mean"))
  #ggsave(filename = paste0(x, "/",i,"_",x,"_", unique(Marker$List)[[i]], "_mean.pdf"),plot=plot, width = 6, height = 5)
  #for MNG
  ggsave(filename = paste0("/omics/groups/OE0146/internal/Micha/Analysis/MNG/",i,"_",x,"_", unique(Marker$List)[[i]], "_mean.pdf"),plot=plot, width = 6, height = 5)
}

#Plot necrotic Markers from MNG
Marker_MNG <- read.csv("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/Reference_sets/Necrotic_Markers.csv")
Marker_MNG <- Marker_MNG[which(Marker_MNG$List =="Necrotic_Markers"),]
Expr <- spata_obj@data[["s_GB12"]][["scaled"]][which(rownames(spata_obj@data[["s_GB12"]][["scaled"]]) %in% Marker_MNG$Name),]
Expr_mean <- colMeans(Expr)
Expr_mean <- data.frame(cbind(colnames(Expr), Expr_mean))
Expr_mean[,2] <- as.numeric(Expr_mean[,2])
colnames(Expr_mean) <- c("barcodes", "Marker_Necr_MNG")
#Add Marker means
spata_obj <- SPATA2::addFeatures(spata_obj, Expr_mean, overwrite = T)

p3 <- plotSurface(spata_obj, "Marker_Necr_MNG", smooth = T, pt_alpha=0.75)+
  ggpLayerThemeCoords() +
  ggpLayerAxesSI(spata_obj, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
  ggtitle("Marker_Necr_MNG")
ggsave(filename = paste0(x, "/",x,"_Marker_Necr_MNG_mean.pdf"),plot=p3)

#Extract overlapping genes between both lists
Overlap <- Marker[which(Marker$Name %in% Marker_MNG$Name),]
Gene_list <- Overlap$Name

path_sn <- "/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/Combined_Analysis/sn/analysis_merged/"
harmon <- readRDS(paste0(path_sn, "harmon_strict_DR_Annot.rds"))
FeaturePlot(harmon, features = Gene_list)
harmon <- AddModuleScore(object =harmon, features= Gene_list, name = "Necr_Marker", assay ='SCT')
MES2 <- Marker[which(Marker$List=="MES2-type"),2]
MES1 <- Marker[which(Marker$List=="MES1-type"),2]
MES2 <- MES2[which(MES2 %in% harmon@assays[["SCT"]]@var.features)]
MES1 <-MES1[which(MES1 %in% harmon@assays[["SCT"]]@var.features)]
harmon <- AddModuleScore(object =harmon, features= Gene_list, name = "Necr_Marker", assay ='SCT')
harmon <- AddModuleScore(object =harmon, features= MES1, name = "MES1", assay ='SCT')
harmon <- AddModuleScore(object =harmon, features= MES2, name = "MES2", assay ='SCT')

Features <- c("Necr_Marker1", "MES11", "MES21")

for(i in 1:3){
  Feat <- Features[i]
  plot1 <- SCpubr::do_FeaturePlot(harmon, features = Feat, pt.size = 0.25, max.cutoff = 2)
  plot2 <- SCpubr::do_ViolinPlot(harmon, features= Feat, group.by="Subtypes")
  ggsave(plot = plot1,filename = paste0(path_sn,Feat,".pdf"), width = 8, height = 8)
  ggsave(plot = plot2,filename = paste0(path_sn,Feat,"_Vln_plot.pdf"), width = 8, height = 8)
}

p4 <- SCpubr::do_FeaturePlot(harmon, features = c("Necr_Marker1"), pt.size = 0.25, max.cutoff = 2)
p4
ggsave(plot = p4,filename = paste0(path_sn,"Necr_Marker_MNG.pdf"), width = 8, height = 8)

p5 <- SCpubr::do_FeaturePlot(harmon, features = c("MES11"), pt.size = 0.25, max.cutoff = 4)
p5
ggsave(plot = p5,filename = paste0(path_sn,"MES1_Marker.pdf"), width = 8, height = 8)

p6 <- SCpubr::do_FeaturePlot(harmon, features = c("MES21"), pt.size = 0.25, max.cutoff = 2)
p6
ggsave(plot = p6,filename = paste0(path_sn,"MES2_Marker.pdf"), width = 8, height = 8)

Xenium_Marker <- c("ADM", "CD44", "C1R", "VEGFA", "GBE1", "ENO2")

#Annotate MES Subtypes
harmon$Subtypes <- harmon$Annot_v1
harmon$Subtypes[which(harmon$seurat_clusters %in% c(17))] <- "MES2-like tumour"
harmon$Subtypes[which(harmon$seurat_clusters %in% c(20))] <- "MES1-like tumour"

SCpubr::do_DimPlot(harmon, group.by="Subtypes", pt.size = 0.25)
SCpubr::do_DimPlot(harmon, group.by="Annot_v2", pt.size = 0.25)

for(i in 1:length(Xenium_Marker)){
  gene <- Xenium_Marker[i]
  plot1 <- SCpubr::do_FeaturePlot(harmon, features = gene, pt.size = 0.25)
  plot2 <- SCpubr::do_ViolinPlot(harmon, features= gene, group.by="Subtypes")
  ggsave(plot = plot1,filename = paste0(path_sn,gene,".pdf"), width = 8, height = 8)
  ggsave(plot = plot2,filename = paste0(path_sn,gene,"_Vln_plot.pdf"), width = 8, height = 8)
}
