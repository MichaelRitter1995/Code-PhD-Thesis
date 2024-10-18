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
sample_name <- "195096"
x <-  "195096"

spata_obj <- updateSpataObject(readRDS(paste0(pathVis, x, "/SPATA_obj_", x, ".RDS")))

#Plot Annotation
anno <- read.csv(paste0(x,"/IFZ.csv"), fill=TRUE)
names(anno) <- c("barcodes", "Annotation")
spata_obj <- SPATA2::addFeatures(spata_obj, anno, overwrite = T)

nr_of_hist <- anno$Annotation %>% unique() %>% length()
colors_histology <- colorRampPalette(RColorBrewer::brewer.pal(13, "Set3"))(nr_of_hist)
names(colors_histology) <- anno$Annotation %>% unique()

p1 <- 
  plotSurface(spata_obj, color_by="Annotation", pt_alpha=0.5)+
  scale_color_manual(values=colors_histology)+
  ggpLayerThemeCoords() +
  ggpLayerAxesSI(spata_obj, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
  ggtitle("Annotation")
p1
ggsave(filename = paste0(x, "/",x,"_Histology_Annotation.pdf"),plot=p1, width = 6, height = 5)

#Create a List of Markers:
Oligo <- c("MBP","OLIG2","OLIG3","CLDN11","MOG","SOX10","GJB1","MAG","CNP","PLP1")
Neurons <-  c("ENO2","NEFL","NEFM","NEFH","INA","SLC17A7","RBFOX3","SNAP25","SYT1","NEUROD2", "MAP2", "RBFOX1", "SNAP25")

for (i in 1:2) {
  ifelse(i==1, Marker <- Oligo, Marker <- Neurons)
  ifelse(i==1, Name <- "Oligo", Name <- "Neurons")
  Expr <- spata_obj@data[["s_195096"]][["scaled"]][which(rownames(spata_obj@data[["s_195096"]][["scaled"]]) %in% Marker),]
  Expr_mean <- colMeans(Expr)
  Expr_mean <- data.frame(cbind(colnames(Expr), Expr_mean))
  Expr_mean[,2] <- as.numeric(Expr_mean[,2])
  colnames(Expr_mean) <- c("barcodes", paste0(Name,"_mean"))
  spata_obj <- SPATA2::addFeatures(spata_obj, Expr_mean, overwrite = T)
  plot <- plotSurface(spata_obj, paste0(Name,"_mean"), smooth = T, pt_alpha=0.75)+
    ggpLayerThemeCoords() +
    ggpLayerAxesSI(spata_obj, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
    ggtitle(paste0(unique(Name,"_mean")))
  plot
  ggsave(filename = paste0(x, "/",x, "_", Name, "_mean.pdf"),plot=plot, width = 6, height = 5)
}

spata_obj <-barcodesToImageAnnotation(spata_obj,
    barcodes = anno[which(anno$Annotation=="Tumour"),1],
    id = 'Tumour',
    tags = 'Tumour',
    overwrite = TRUE)

ias_layer_bins <- 
  ggpLayerEncirclingIAS(
    object = spata_obj,
    distance = "2.25mm", 
    n_bins_circle = 20, 
    id = "Tumour_1",
    line_size = 1
  )

bcsp_dist <- getCCD(spata_obj, unit = "mm")

p2 <- plotSurfaceIAS(
  object =spata_obj, 
  id = "Tumour_1", # the ID of the image annotation of interest
  distance = "3mm", # covers the whole sample
  binwidth = bcsp_dist, # two layers of spots per bine
  display_bins_angle = FALSE,
  n_bins_circle = 15, )
)
p2
ggsave(filename = paste0(x, "/",x,"_Expected_CNVs_Screening_Annotation.pdf"),plot=p2, width = 6, height = 5)
