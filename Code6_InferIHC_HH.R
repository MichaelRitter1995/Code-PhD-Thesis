library(ggplot2)
library(SPATA2)
library(dplyr)
library(stringr)
setwd("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")

#############Define missing functions################

#Define runSuperresolution function

runSuperresolution <- function(object, path, nrep=1000, burn.in=100, features){
  
  nr_clusters <- joinWithFeatures(object, features = features) %>% pull(!!sym(features)) %>% unique() %>% length()
  
  library(BayesSpace)
  sce <- BayesSpace::readVisium(path)
  sce <- sce[, colSums(SingleCellExperiment::counts(sce)) > 0]
  confuns::give_feedback(msg ="Preprocess Data")
  sce <- BayesSpace::spatialPreprocess(sce, platform="Visium")
  confuns::give_feedback(msg ="Run Spatial Clustering")
  sce <- BayesSpace::spatialCluster(sce, q=nr_clusters, nrep=nrep, burn.in=burn.in)
  confuns::give_feedback(msg ="Run SuperResolution")
  sce.enhanced <- BayesSpace::spatialEnhance(sce,q=nr_clusters, platform = c("Visium"),nrep=nrep, burn.in=burn.in)
  
  
  return(list(sce, sce.enhanced))

}


#Define plotInferedIHC function

plotInferedIHC <- function(features,
                           object, 
                           sce,
                           sce.enhanced,
                           pointsize=6,
                           limits=c(0.4, 1),
                           range=c(0, 0.3),
                           type=c("dots", "voronoi")[1]){
  
  
  
  #object <- object %>% SPATA2::flipAll(axis="h")
  x.r <- object %>% getCoordsDf() %>% pull(x) %>% range()
  y.r <-  object %>% getCoordsDf() %>% pull(y)%>% range()
  
  
  
  p <- plotSurface(object, pt_alpha = 0)
  sce.enhanced <- BayesSpace <- enhanceFeatures(sce.enhanced, sce , feature_names=features)
  data <- BayesSpace::featurePlot(sce.enhanced, feature = features)
  
  data$data$fill <- SPATA2::hlpr_normalize_vctr(data$data$fill)
  
  x <- data$data$y.vertex 
  y <- data$data$x.vertex 
  
  img_range <- SPATA2::getImageRange(object)$x
  x <- img_range[2] - x + img_range[1]
  
  data$data$y <- y %>% scales::rescale(.,y.r)
  data$data$x <- x %>% scales::rescale(.,x.r)
  
  if(type=="dots"){
    p <- 
      p+
      scattermore::geom_scattermore(data=data$data , 
                                    mapping=aes(y=y,# %>% scales::rescale(.,x.r), 
                                                x=x,# %>% scales::rescale(.,y.r), 
                                                color = fill,
                                                alpha=fill), pointsize = pointsize,
                                    pixels = c(2048, 2048))+
      scale_color_gradientn(colors = IHC(50), limits=limits, oob=scales::squish)+
      scale_alpha(limits=limits, oob=scales::squish, range=range)
    
  }else{
    p <- 
      ggplot()+
      ggforce::geom_voronoi_tile(data=data$data , 
                                 mapping=aes(x=x.vertex %>% scales::rescale(.,x.r), 
                                             y=y.vertex %>% scales::rescale(.,y.r),
                                             group= -1L,
                                             fill = fill,alpha=fill))+
      scale_fill_gradientn(colors = IHC(50), limits=limits, oob=scales::squish)+
      scale_alpha(limits=limits, oob=scales::squish, range=range)
  }
  
  
  p+
    SPATA2::ggpLayerAxesSI(object)+
    SPATA2::ggpLayerThemeCoords()+
    coord_fixed()
  
  
  
  
  
}

###Define colour palette

IHC <- colorRampPalette(c("#FFFFFF", "#DDD4CF","#E1C07E", "#9A401C", "#6F2C36", "#3B1112"))

path_outs <- "/omics/groups/OE0146/internal/Christina/Visium_seq/data_spaceranger/"


#########################Start Analysis###############################################


ID_list=data.frame(c("GB_STX2", "GB_STX7", "GB_STX1", "GB_STX6", "GB_STX4",
                     "195096", "GB2", "GB6",  "GB11", "GB13", "GB3",
                     "N20_2382", "GB1", "GB10", "GB12", "GB_STX3", "N21_2972",
                     "GB9", "GB5", "N20_1260_IB", "GB7"))



sample_name <- "195096"
  

for (i in 1:length(ID_list)) {
  sample_name = ID_list[[i]]
  print(sample_name)
  spata_obj <- updateSpataObject(readRDS(paste0(sample_name, "/SPATA_obj_", sample_name, ".RDS")))
  IHC_inf <- runSuperresolution(spata_obj, paste0(path_outs, sample_name, "/outs"), features = "bayes_space")
  p1 <- plotInferedIHC(spata_obj, IHC_inf[[1]], IHC_inf[[2]], features="GFAP", type="dots", limits=c(0.2, 0.8), pointsize=8)
  ggsave(filename = paste0("7_", sample_name, "_GFAP_inferredIHC.pdf"), plot = p1, path= paste0(sample_name, "/"), height = 6, width = 6)
  
  p2 <- plotInferedIHC(spata_obj, IHC_inf[[1]], IHC_inf[[2]], features="MKI67", type="dots", limits=c(0.2, 0.8), pointsize=8)
  ggsave(filename = paste0("8_", sample_name, "_Ki67_inferredIHC.pdf"), plot = p2, path= paste0(sample_name, "/"), height = 6, width = 6)
  
  p3 <- plotInferedIHC(spata_obj, IHC_inf[[1]], IHC_inf[[2]], features="RBFOX3", type="dots", limits=c(0.2, 0.8), pointsize=8)
  ggsave(filename = paste0("9_", sample_name, "_NeUN_inferredIHC.pdf"), plot = p3, path= paste0(sample_name, "/"), height = 6, width = 6)

  p4 <- plotInferedIHC(spata_obj, IHC_inf[[1]], IHC_inf[[2]], features="MBP", type="dots", limits=c(0.2, 0.8), pointsize=8)

  p5 <- plotSurface(spata_obj, color_by="Oligodendrocyte", pt_alpha=0.5)+
    ggpLayerThemeCoords() +
    ggpLayerAxesSI(spata_obj, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
    ggtitle("Oligodendrocytes Deconvolution")
 
  spata_obj@fdata[[paste0("s_", sample_name)]][["7/10_Score"]] <- 1+(spata_obj@fdata[[paste0("s_", sample_name)]][["Chr7"]]-spata_obj@fdata[[paste0("s_", sample_name)]][["Chr10"]])/2 
  
  p6 <- plotSurface(spata_obj, color_by="7/10_Score", pt_alpha=0.5, smooth=T)+
   ggpLayerThemeCoords() +
   ggpLayerAxesSI(spata_obj, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
   ggtitle("7/10 Score")
 
  plot_panel <- p5+p4+p6
 
  ggsave(filename = paste0("10_", sample_name, "_Oligodendrocytes_Deconv_IHC.pdf"), plot = plot_panel, path= paste0(sample_name, "/"), height = 6, width = 18)
 
}

sample_name <- "195096"
x <- "GB_STX6"

path_Chr <- "/omics/groups/OE0146/internal/Christina/Visium_seq/analysis/GBM/"
seurat_obj <- readRDS(paste0(path_Chr, sample_name, "/", sample_name, "_seurat_obj.rds"))

for (i in 1:nrow(ID_list)) {
  x = ID_list[i,1]
  ifelse(x=="GB_STX3" | x=="N21_2972", sample_name <- "GB_STX3_2972", sample_name <- x)
  print(sample_name)
  seurat_obj <- readRDS(paste0(path_Chr, sample_name, "/", sample_name, "_seurat_obj.rds"))
  seurat_obj@active.assay <- "SCT"
  p7 <- SpatialPlot(seurat_obj, feature="EGFR")+scale_fill_gradientn(colours = viridis::mako(10), limits = c(0,5))
  ggsave(filename = paste0("10_", x, "_EGFR_Expression.pdf"), plot = p7, path= paste0(x, "/"), height = 6, width = 6)
  
  }

p7
