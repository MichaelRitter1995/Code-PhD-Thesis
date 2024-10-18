###################Seurat##########################

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(purrr)
library(Seurat)
library(SPATA2)
library(hdf5r)
library(stringr)
library(RColorBrewer)
reticulate::use_condaenv('r-reticulate')
setwd("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/Combined_Analysis/visium")


path_Vis <- "/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium/aligned/Xenium/GB/"

#Seu_obj1 <- Load10X_Spatial(paste0(path_Vis,"Xenium_CytAssist_GB1/outs"))
#Seu_obj2 <- Load10X_Spatial(paste0(path_Vis,"Xenium_CytAssist_GB2/outs"))

#Seu_obj1 <- SCTransform(Seu_obj1, assay = "Spatial", verbose = FALSE)
#saveRDS(Seu_obj1, file = "Xenium_CytAssist_GB1/Seurat_obj_Xenium_CytAssist_GB1.rds")

#spata_obj1 <- 
#  SPATA2::asSPATA2(
#    object = Seu_obj1,
#    sample_name = "Xenium_Visium_GB1",
#    spatial_method = "Visium"
#  )

#saveRDS(spata_obj1, file = "Xenium_CytAssist_GB1/SPATA_obj_Xenium_CytAssist_GB1.rds")

#Seu_obj2 <- SCTransform(Seu_obj2, assay = "Spatial", verbose = FALSE)
#saveRDS(Seu_obj2, file = "Xenium_CytAssist_GB2/Seurat_obj_Xenium_CytAssist_GB2.rds")


#spata_obj2 <- SPATA2::asSPATA2(
#    object = Seu_obj2,
#    sample_name = "Xenium_Visium_GB2",
#    spatial_method = "Visium"
#  )

#saveRDS(spata_obj2, file = "Xenium_CytAssist_GB2/SPATA_obj_Xenium_CytAssist_GB2.rds")

spata_obj1 <- readRDS("Xenium_CytAssist_GB1/SPATA_obj_Xenium_CytAssist_GB1.rds")
spata_obj2 <- readRDS("Xenium_CytAssist_GB2/SPATA_obj_Xenium_CytAssist_GB2.rds")

seq1 <- read.csv("Xenium_CytAssist_GB1/Xenium_CytAssist_GB1_Mutations.csv")
seq2 <- read.csv("Xenium_CytAssist_GB2/Xenium_CytAssist_GB2_Mutations.csv")

spata_obj1@fdata[["Xenium_Visium_GB1"]][["Mutations"]] <- seq1$Mutations
spata_obj2@fdata[["Xenium_Visium_GB2"]][["Mutations"]] <- seq2$Mutations

p <- plotSurface(spata_obj1, color_by="Mutations")+
  scale_color_manual(values=brewer.pal(n = 8, name = "Set3")[c(4,5,9)]) +
  ggpLayerThemeCoords() +
  ggpLayerAxesSI(spata_obj1, unit = "mm", breaks = str_c(1:11, "mm"), add_labs = TRUE)
ggsave(filename = "Mutations/Xenium_Visium_GB1_Annot.pdf", plot = p, height = 6, width = 9)

p0 <- plotSurface(spata_obj2, color_by="Mutations") +
scale_color_manual(values=brewer.pal(n = 8, name = "Set3")[c(9,4,9,6)]) +
  ggpLayerThemeCoords() +
  ggpLayerAxesSI(spata_obj2, unit = "mm", breaks = str_c(1:11, "mm"), add_labs = TRUE)
ggsave(filename = "Mutations/Xenium_Visium_GB2_Annot.pdf", plot = p0, height = 6, width = 9)

Genes <- c("TP53", "PTEN", "KMT2C", "ARID1A")
Samples <- list(spata_obj1, spata_obj2)

#Plot the wildtype and mutation specific probes on the sample

for(j in 1:2){
  y <- Samples[[j]]
  for(i in 1:length(Genes)){
    x <- Genes[[i]]
    ifelse(x == "KMT2C",
           x1 <- "KMT2C ",
           x1 <- x)
    if(paste0("mut",x1) %in% rownames(y@data[[y@samples]][["counts"]])){
    p1 <- plotSurface(y, color_by = paste0("mut",x1), smooth=T, pt_size=0.75) +
      scale_fill_gradientn(colours=viridis::magma(10)) +
      ggpLayerThemeCoords() +
      ggpLayerAxesSI(y, unit = "mm", breaks = str_c(1:11, "mm"), add_labs = TRUE)
    ggsave(filename = paste0("Mutations/mut",x1,"_spatial_GB",j,".pdf"), plot = p1, height = 6, width = 9)
    }
    if(paste0("wt",x) %in% rownames(y@data[[y@samples]][["counts"]])){
    p2 <- plotSurface(y, color_by = paste0("wt",x), smooth=T, pt_size=0.75) +
      scale_fill_gradientn(colours=viridis::magma(10))+
      ggpLayerThemeCoords() +
      ggpLayerAxesSI(y, unit = "mm", breaks = str_c(1:11, "mm"), add_labs = TRUE)
    ggsave(filename = paste0("Mutations/wt",x1,"_spatial_GB",j,".pdf"), plot = p2, height = 6, width = 9)
    }
    p3 <- plotSurface(y, color_by = paste0(x), smooth=T, pt_size=0.75) +
      scale_fill_gradientn(colours=viridis::magma(10)) + 
      ggpLayerThemeCoords() +
      ggpLayerAxesSI(y, unit = "mm", breaks = str_c(1:11, "mm"), add_labs = TRUE)
    ggsave(filename = paste0("Mutations/",x1,"_spatial_GB",j,".pdf"), plot = p3, height = 6, width = 9)
  }
}

ggsave(filename = paste0("Mutations/KMT2C_spatial_GB2.pdf"), plot = p3, height = 6, width = 9)

#Quantify the probes

Seu_obj1 <- readRDS("Xenium_CytAssist_GB1/Seurat_obj_Xenium_CytAssist_GB1.rds")
Seu_obj2 <- readRDS("Xenium_CytAssist_GB2/Seurat_obj_Xenium_CytAssist_GB2.rds")

list_of_data.frames <- lapply(Genes, function(x){
  
  #Extract Reads
  ifelse(x == "KMT2C",
    x1 <- "KMT2C ",
  x1 <- x)
  
  total1 <- Seu_obj1@assays[["SCT"]]@counts[which(rownames(Seu_obj1@assays[["SCT"]]) == x),]
  total2 <- Seu_obj2@assays[["SCT"]]@counts[which(rownames(Seu_obj2@assays[["SCT"]]) == x),]
  
  ifelse(paste0("mut",x1) %in% rownames(Seu_obj1),
  mut1 <- Seu_obj1@assays[["SCT"]]@counts[which(rownames(Seu_obj1@assays[["SCT"]]) == paste0("mut",x1)),],
  mut1 <- rep(0, ncol(Seu_obj1))
  )
  
  ifelse(paste0("mut",x1) %in% rownames(Seu_obj2),
  mut2 <- Seu_obj2@assays[["SCT"]]@counts[which(rownames(Seu_obj2@assays[["SCT"]]) == paste0("mut",x1)),],
  mut2 <- rep(0, ncol(Seu_obj2))
  )
  
  ifelse(paste0("mut",x) %in% rownames(Seu_obj1),
  wt1 <- Seu_obj1@assays[["SCT"]]@counts[which(rownames(Seu_obj1@assays[["SCT"]]) == paste0("wt",x)),],
  wt1 <- rep(0, ncol(Seu_obj1))
  )
  
  ifelse(paste0("wt",x) %in% rownames(Seu_obj2),
  wt2 <- Seu_obj2@assays[["SCT"]]@counts[which(rownames(Seu_obj2@assays[["SCT"]]) == paste0("wt",x)),],
  wt2 <- rep(0, ncol(Seu_obj2))
  )
  
  #Merge raw data
  raw1 <- data.frame(cbind(names(total1), total1, wt1, mut1))
  raw2 <- data.frame(cbind(names(total2), total2, wt2, mut2))
  colnames(raw1) <- c("Barcode", "total_reads", "wt_reads", "mut_reads")
  colnames(raw2) <- c("Barcode", "total_reads", "wt_reads", "mut_reads")
  fin <- rbind(merge(raw1, seq1, by="Barcode"), merge(raw2, seq2, by="Barcode"))
  
  #Clean-up raw data.frame
  fin <- cbind(c(rep("Slide1", ncol(Seu_obj1)), rep("Slide2", ncol(Seu_obj2))), fin)
  colnames(fin)[1] <- "Slide"
  
  #account for 0 and transform for calclulation
  fin[which(fin$total_reads ==0),3] <- 1
  fin[,3:5] <- lapply(fin[,3:5], as.numeric)
  
  #Calculate abundance
  fin[,ncol(fin)+1] <- fin$wt_reads/fin$total_reads
  fin[,ncol(fin)+1] <- fin$mut_reads/fin$total_reads
  
  #Finalise data.frame
  fin <- fin[which(fin$Mutations !="Not_sequenced"),]
  fin[,ncol(fin)+1] <- rep(x, nrow(fin))
  colnames(fin)[7:9] <- c("norm_wt", "norm_mut", "Gene")
  fin
})

TP53 <- list_of_data.frames[[1]]
TP53 <- TP53[which(TP53$Mutations!=""),]
TP53_mut <- TP53[,c(6,8)]
TP53_wt <- TP53[,6:7]
names(TP53_mut) <- names(TP53_wt)
TP53_fin <- rbind(TP53_wt, TP53_mut)
TP53_fin[,3] <- c(rep("wildtype", nrow(TP53_fin)/2),  rep("mutation", nrow(TP53_fin)/2))
TP53_fin[,4] <- TP53_fin[,1]
TP53_fin$Mutations <- "wtTP53"
TP53_fin$Mutations[TP53_fin$V4 == "PTEN_TP53"] <- "mutTP53"
TP53_fin$norm_wt <- as.numeric(TP53_fin$norm_wt)

PTEN <- list_of_data.frames[[2]]
PTEN <- PTEN[which(PTEN$Mutations!=""),]
PTEN_mut <- PTEN[,c(6,8)]
PTEN_wt <- PTEN[,6:7]
names(PTEN_mut) <- names(PTEN_wt)
PTEN_fin <- rbind(PTEN_wt, PTEN_mut)
PTEN_fin[,3] <- c(rep("wildtype", nrow(PTEN_fin)/2),  rep("mutation", nrow(PTEN_fin)/2))
PTEN_fin[,4] <- PTEN_fin[,1]
PTEN_fin$Mutations <- "wtPTEN"
PTEN_fin$Mutations[PTEN_fin$V4 == "PTEN_TP53"] <- "mutPTEN"
PTEN_fin$norm_wt <- as.numeric(PTEN_fin$norm_wt)

KMT2C <- list_of_data.frames[[3]]
KMT2C <- KMT2C[which(KMT2C$Mutations!=""),]
KMT2C_mut <- KMT2C[,c(6,8)]
KMT2C_wt <- KMT2C[,6:7]
names(KMT2C_mut) <- names(KMT2C_wt)
KMT2C_fin <- rbind(KMT2C_wt, KMT2C_mut)
KMT2C_fin[,3] <- c(rep("wildtype", nrow(KMT2C_fin)/2),  rep("mutation", nrow(KMT2C_fin)/2))
KMT2C_fin[,4] <- KMT2C_fin[,1]
KMT2C_fin$Mutations <- "mutKMT2C"
KMT2C_fin$norm_wt <- as.numeric(KMT2C_fin$norm_wt)

ARID1A <- list_of_data.frames[[4]]
ARID1A <- ARID1A[which(ARID1A$Mutations!=""),]
ARID1A_mut <- ARID1A[,c(6,8)]
ARID1A_wt <- ARID1A[,6:7]
names(ARID1A_mut) <- names(ARID1A_wt)
ARID1A_fin <- rbind(ARID1A_wt, ARID1A_mut)
ARID1A_fin[,3] <- c(rep("wildtype", nrow(ARID1A_fin)/2),  rep("mutation", nrow(ARID1A_fin)/2))
ARID1A_fin[,4] <- ARID1A_fin[,1]
ARID1A_fin$Mutations <- "wtARID1A"
ARID1A_fin$Mutations[ARID1A_fin$V4 == "ARID1A"] <- "mutARID1A"
ARID1A_fin$norm_wt <- as.numeric(ARID1A_fin$norm_wt)

p1 <- ggplot(data=TP53_fin, aes(x=Mutations, y=norm_wt, colour=V3))+geom_jitter(position=position_jitterdodge(dodge.width = 1)) +theme_classic() + scale_color_manual(values=brewer.pal(n = 8, name = "Set3")[4:5])
p2 <- ggplot(data=PTEN_fin, aes(x=Mutations, y=norm_wt, colour=V3))+geom_jitter(position=position_jitterdodge(dodge.width = 1))+theme_classic() + scale_color_manual(values=brewer.pal(n = 8, name = "Set3")[4:5])
p3 <- ggplot(data=KMT2C_fin, aes(x=Mutations, y=norm_wt, colour=V3))+geom_jitter(position=position_jitterdodge(dodge.width = 1))+theme_classic() + scale_color_manual(values=brewer.pal(n = 8, name = "Set3")[4:5])
p4 <- ggplot(data=ARID1A_fin, aes(x=Mutations, y=norm_wt, colour=V3))+geom_jitter(position=position_jitterdodge(dodge.width = 1))+theme_classic() + scale_color_manual(values=brewer.pal(n = 8, name = "Set3")[4:5])
ggsave(filename = paste0("TP53_mutation.pdf"), plot = p1, height = 6, width = 9)
ggsave(filename = paste0("PTEN_mutation.pdf"), plot = p2, height = 6, width = 9)
ggsave(filename = paste0("KMT2C_mutation.pdf"), plot = p3, height = 6, width = 9)
ggsave(filename = paste0("ARID1A_mutation.pdf"), plot = p4, height = 6, width = 9)

ARID1A_fin[,5] <- "wt_pwt"
ARID1A_fin$V5[ARID1A_fin$Mutations == "wtARID1A" & ARID1A_fin$V3 == "mutation"] <- "wt_pmut"
ARID1A_fin$V5[ARID1A_fin$Mutations == "mutARID1A" & ARID1A_fin$V3 == "wildtype"] <- "mut_pwt"
ARID1A_fin$V5[ARID1A_fin$Mutations == "mutARID1A" & ARID1A_fin$V3 == "mutation"] <- "mut_pmut"
kruskal.test(norm_wt~V5, data=ARID1A_fin)
dunn_test(norm_wt~V5, data=ARID1A_fin, p.adjust.method = "bonferroni")

TP53_fin[,5] <- "wt_pwt"
TP53_fin$V5[TP53_fin$Mutations == "wtTP53" & TP53_fin$V3 == "mutation"] <- "wt_pmut"
TP53_fin$V5[TP53_fin$Mutations == "mutTP53" & TP53_fin$V3 == "wildtype"] <- "mut_pwt"
TP53_fin$V5[TP53_fin$Mutations == "mutTP53" & TP53_fin$V3 == "mutation"] <- "mut_pmut"
kruskal.test(norm_wt~V5, data=TP53_fin)
dunn_test(norm_wt~V5, data=TP53_fin, p.adjust.method = "bonferroni")

#Due to probe malfunctions, mutation present in all samples only two groups need to be compared

wilcox.test(norm_wt~Mutations, data=PTEN_fin)
wilcox.test(norm_wt~V3, data=KMT2C_fin)
