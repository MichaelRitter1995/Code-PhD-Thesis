#Plots Thesis

library(SPATA2)
library(stringr)
library(ggplot2)

setwd("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")
pathVis <- ("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")

#Change sample name to the right sample
sample_name <- "195096"
x <-  "195096"

#load spata object
spata_obj <- updateSpataObject(readRDS(paste0(pathVis, x, "/SPATA_obj_", x, ".RDS")))


#Plot HE
p1 <- 
  plotSurface(spata_obj, pt_alpha=0)+
  ggpLayerThemeCoords() +
  ggpLayerAxesSI(spata_obj, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
  ggtitle("H&E")
p1

#Plot Annotation
anno <- read.csv(paste0(x,"/Annotated.csv"), fill=TRUE)
names(anno) <- c("barcodes", "Histology")
anno[anno$Histology=="", ]$Histology <- "Unknown"
spata_obj <- SPATA2::addFeatures(spata_obj, anno, overwrite = T)

nr_of_hist <- anno$Histology %>% unique() %>% length()
colors_histology <- colorRampPalette(RColorBrewer::brewer.pal(13, "Set3"))(nr_of_hist)
names(colors_histology) <- anno$Histology %>% unique()

p2 <- 
  plotSurface(spata_obj, color_by="Histology", pt_alpha=0.5)+
  scale_color_manual(values=colors_histology)+
  ggpLayerThemeCoords() +
  ggpLayerAxesSI(spata_obj, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
  ggtitle("Annotated Histology")
p2


#Plot 7/10 Score
combined_score <- data.frame(cbind(spata_obj@cnv[[paste0("s_",sample_name)]][["cnv_df"]][["barcodes"]], as.numeric(1+(spata_obj@cnv[[paste0("s_",sample_name)]][["cnv_df"]][["Chr7"]]-spata_obj@cnv[[paste0("s_",sample_name)]][["cnv_df"]][["Chr10"]])/2)))
names(combined_score) <- c("barcodes", "combined_score")
combined_score$combined_score <- as.numeric(combined_score$combined_score)
spata_obj <- SPATA2::addFeatures(spata_obj, combined_score, overwrite = T)

p3 <- 
  plotSurface(spata_obj, "combined_score", smooth = T, pt_alpha=0.75)+
  ggpLayerThemeCoords() +
  ggpLayerAxesSI(spata_obj, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
  ggtitle("Combined 7/10 score")+
  scale_colour_gradientn(colours = c(col[1], "white", col[9]), limits=c(0.99,1.05), name="Combined Chr 7 and 10 score")
p3

