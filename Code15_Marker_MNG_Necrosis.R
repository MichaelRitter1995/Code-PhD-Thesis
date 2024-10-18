library(SPATA2)
library(rstatix)

pathVis <- ("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")

ID_list=read.table("ID_list.txt")
ID_list <- ID_list[2:nrow(ID_list),1]

Marker_MNG <- read.csv("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/Reference_sets/Necrotic_Markers.csv")
Marker_MNG <- Marker_MNG[which(Marker_MNG$List =="Necrotic_Markers"),]

list_of_dataframes=lapply(ID_list,function (x){ 
                            spata_obj <- updateSpataObject(readRDS(paste0(pathVis, x, "/SPATA_obj_", x, ".RDS")))
                            Expr <- spata_obj@data[[paste0("s_",x)]][["scaled"]][which(rownames(spata_obj@data[[paste0("s_",x)]][["scaled"]]) %in% Marker_MNG$Name),]
                            Expr_mean <- colMeans(Expr)
                            Expr_mean <- data.frame(cbind(colnames(Expr), Expr_mean))
                            Expr_mean[,2] <- as.numeric(Expr_mean[,2])
                            colnames(Expr_mean) <- c("barcodes", "Marker_Necr_MNG")
                            
                            Annotation <- read.csv(paste0(x,"/Annotated.csv"), fill=TRUE)
                            Annotation[,ncol(Annotation)+1]=rep(x, nrow(Annotation)) 
                            colnames(Annotation)=c("barcodes", "Annotation", "Sample ID")
                            fin <- merge(Annotation, Expr_mean, by="barcodes")
                          })

MNG_Necr_all=dplyr::bind_rows(list_of_dataframes)
MNG_Necr_all[,ncol(MNG_Necr_all)+1]<-ifelse(MNG_Necr_all$Annotation %in% c("IFZ", "Necrosis", "Main Tumour", "Healthy"), MNG_Necr_all$Annotation, "Other")

MNG_Necr_per_sample <- data.frame()

for(i in 1:length(unique(MNG_Necr_all$`Sample ID`))){
  samp <- unique(MNG_Necr_all$`Sample ID`)[[i]]
  subset <- MNG_Necr_all[which(MNG_Necr_all$`Sample ID`==samp),]
  for(j in 1:length(unique(subset$Annotation))){
    st <- unique(subset$Annotation)[[j]]
    Hist <- subset[which(subset$Annotation == st),]
    Hist <- na.omit(Hist)
    if(nrow(Hist)>1){
      Mean_hist <- mean(sapply(Hist[,4], as.numeric))
      Mean_fin <-cbind(Hist[1,2:3], t(Mean_hist), Hist[1,5])
    } else {
      Mean_fin <-Hist[1,c(2:5)]
      colnames(Mean_fin) <- colnames(MNG_Necr_per_sample)
    }
    MNG_Necr_per_sample <- rbind(MNG_Necr_per_sample, Mean_fin)
  }
}
colnames(MNG_Necr_per_sample) <- c("Annotation", "ID", "Mean_scaled_Expression", "Histology")
t1 <- ggplot(MNG_Necr_per_sample)+theme_bw()+labs(x=NA, y= "Mean_scaled_Expression", main="Necrotic Marker Expression")
t2 <- geom_boxplot(aes(x=Histology, y=Mean_Expression, color=Histology), outlier.shape = NA)
t3 <- geom_jitter(aes(x=Histology, y=Mean_Expression, color=Histology), position = position_jitterdodge())


t1+t3+t2+theme_classic()
for(i in 1:length(unique(MNG_Necr_per_sample$Histology))){
  test <- shapiro.test(MNG_Necr_per_sample[which(MNG_Necr_per_sample$Histology==unique(MNG_Necr_per_sample$Histology)[i]),3])
  print(test)
  }

MNG_Necr_per_sampleS <- MNG_Necr_per_sample[which(MNG_Necr_per_sample$Histology != "Other"),]

stat.test <- MNG_Necr_per_sampleS %>%
  t_test(Mean_scaled_Expression ~ Histology) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test

kruskal.test(Mean_scaled_Expression~Histology, data=MNG_Necr_per_sampleS)
dunn_test(Mean_scaled_Expression~Histology, data=MNG_Necr_per_sampleS, p.adjust.method = "bonferroni")


