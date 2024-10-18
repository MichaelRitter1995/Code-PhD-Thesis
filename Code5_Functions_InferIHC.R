#Load missing functions for inferIHC

###add_hull

add_hull <- function(object,plot,pt.b=6,pt.w=5){
  
  data <- plot$data
  dim_names <- names(data)
  
  p <- 
    ggplot(data)+
    scattermore::geom_scattermore(data=SPATA2::getCoordsDf(object), mapping = aes(x,y),pointsize = pt.b, color="black")+
    scattermore::geom_scattermore(data=SPATA2::getCoordsDf(object), mapping = aes(x,y),pointsize = pt.w, color="white")
  
  p <- 
    p+plot$layers[[1]]+
    SPATA2::ggpLayerAxesSI(object)
  
  return(p)
  
}



###runSuperresolution

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



###plitInferedIHC

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
