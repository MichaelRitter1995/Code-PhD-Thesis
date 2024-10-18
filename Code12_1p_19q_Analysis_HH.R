library(SPATA2)
library(ggplot2)

setwd("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")

pathVis <- ("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")

### Only STX included in analysis ###

#Create dataframe including sample name, ID and Entity
ID_list=data.frame(cbind(c("GB_STX7", "GB_STX1", "GB_STX6", "GB_STX4", "GB_STX3", "N21_2972", "GB_STX2",
                           "AC_NOS_STX1", "AC_STX6", "AC_STX3", "AC_STX4","AC_STX5","AC_STX7","AC_STX1", "AC_STX2",
                           "OG_STX4", "OG_STX1", "OG_STX2", "OG_STX3",
                           "N20_1060", "N19_2708", "N21_377", "PA_STX1"),
                         c("GB_STX7", "GB_STX1", "GB_STX6", "GB_STX4", rep("GB_STX3_2972",2), "GB_STX2",
                           "AC_NOS_STX1", "AC_STX6", rep("AC", 4), rep("AC_OG", 3), rep("OG", 3), rep("MG", 4)),
                         c(rep("Glioblastoma",7), rep("Astrocytoma",8), rep("Oligodendroglioma",4), rep("Midline Glioma", 3), "Astrocytoma")))
colnames(ID_list) <- c("ID","Slide", "Entity")

CNV_all <- data.frame()
CNV_med <- data.frame(matrix(ncol=4))

#Cycle through all samples

for(i in 1:nrow(ID_list)){
      #Define names
      sample_name <-ID_list[i,1]
      print(sample_name)
      x <- ID_list[i,2]
      x1 <- x
      
      #Adjust variables for multiple patients per sample
      file_name <- ifelse(ID_list[i,2] %in% list("AC", "OG", "AC_OG", "MG"), paste0(ID_list[i,2],"/",ID_list[i,1]), paste0(ID_list[i,2]))
      if(sample_name=="GB_STX3" | sample_name=="N21_2972"){
        file_name <- sample_name
        x1 <- sample_name
      }
      
      #Load Annotation and CNV objects
      spata_obj <- updateSpataObject(readRDS(paste0(pathVis, file_name, "/SPATA_obj_", x1, ".RDS")))
      bc_CNV <- read.csv(paste0(file_name,"/CNV.csv"))
      
      #Subset and get number of rows
      spata_sub <- subsetByBarcodes(spata_obj, barcodes = bc_CNV$Barcode)
      n <- nrow(bc_CNV)
      
      #Extract chromosomal alterations
      Chr7 <- data.frame(spata_sub@fdata[[paste0("s_", x)]][["Chr7"]], rep("Chr7", n), rep(ID_list[i,3], n), rep(sample_name, n))
      colnames(Chr7) <- c("Score", "Chr", "Entity", "Sample")
      Chr10 <-data.frame(spata_sub@fdata[[paste0("s_", x)]][["Chr10"]], rep("Chr10", n), rep(ID_list[i,3], n), rep(sample_name, n))
      colnames(Chr10) <- c("Score", "Chr", "Entity", "Sample")
      Chr1p <- data.frame(spata_sub@fdata[[paste0("s_", x)]][["Chr1p"]], rep("1p", n), rep(ID_list[i,3], n), rep(sample_name, n))
      colnames(Chr1p) <- c("Score", "Chr", "Entity", "Sample")
      Chr1q <- data.frame(spata_sub@fdata[[paste0("s_", x)]][["Chr1q"]], rep("1q", n), rep(ID_list[i,3], n), rep(sample_name, n))
      colnames(Chr1q) <- c("Score", "Chr", "Entity", "Sample")
      Chr19p <- data.frame(spata_sub@fdata[[paste0("s_", x)]][["Chr19p"]], rep("19p", n), rep(ID_list[i,3], n), rep(sample_name, n))
      colnames(Chr19p) <- c("Score", "Chr", "Entity", "Sample")
      Chr19q <- data.frame(spata_sub@fdata[[paste0("s_", x)]][["Chr19q"]], rep("19q", n), rep(ID_list[i,3], n), rep(sample_name, n))
      colnames(Chr19q) <- c("Score", "Chr", "Entity", "Sample")
      
      #merge into one dataframe
      ifelse(nrow(CNV_all)==0, 
             CNV_all <- rbind(Chr7, Chr10, Chr1p, Chr1q, Chr19p, Chr19q), 
             CNV_all <- rbind(CNV_all, Chr7, Chr10, Chr1p, Chr1q, Chr19p, Chr19q)
                              )
      #Clacluate median across all spots
      CNV_med[nrow(CNV_med)+1,] <- c(median(Chr7[,1]), "Chr7", ID_list[i,3], sample_name)
      CNV_med[nrow(CNV_med)+1,] <- c(median(Chr10[,1]), "Chr10", ID_list[i,3], sample_name)
      CNV_med[nrow(CNV_med)+1,] <- c(median(Chr1p[,1]), "1p", ID_list[i,3], sample_name)
      CNV_med[nrow(CNV_med)+1,] <- c(median(Chr1q[,1]), "1q", ID_list[i,3], sample_name)
      CNV_med[nrow(CNV_med)+1,] <- c(median(Chr19q[,1]), "19q", ID_list[i,3], sample_name)
      CNV_med[nrow(CNV_med)+1,] <- c(median(Chr19p[,1]), "19p", ID_list[i,3], sample_name)
}

###Load Healthy Control

Healthy_list=data.frame(cbind(c("GB_STX7", "GB_STX1", "GB_STX6", "GB_STX4", "N21_2972", "GB_STX2",
                           "AC_NOS_STX1", "AC_STX6"),
                         c("GB_STX7", "GB_STX1", "GB_STX6", "GB_STX4", "GB_STX3_2972", "GB_STX2",
                           "AC_NOS_STX1", "AC_STX6"),
                         c(rep("Healthy",8))))
colnames(Healthy_list) <- c("ID","Slide", "Entity")

#Process like samples above

for(i in 1:nrow(Healthy_list)){
  sample_name <- Healthy_list[i,1]
  print(sample_name)
  x <- Healthy_list[i,2]
  file_name <- ifelse(x %in% list("AC", "OG", "AC_OG", "MG"), paste0(x,"/",Healthy_list[i,1]), paste0(x))
  if(sample_name=="GB_STX3" | sample_name=="N21_2972"){
    file_name <- sample_name
    x <- "GB_STX3_2972"
  }
  spata_obj <- updateSpataObject(readRDS(paste0(pathVis, file_name, "/SPATA_obj_", sample_name, ".RDS")))
  
  bc_CNV <- read.csv(paste0(file_name,"/Reference.csv"))
  spata_sub <- subsetByBarcodes(spata_obj, barcodes = bc_CNV$Barcode)
  n <- nrow(bc_CNV)
  Chr7 <- data.frame(spata_sub@fdata[[paste0("s_", x)]][["Chr7"]], rep("Chr7", n), rep(Healthy_list[i,3], n), rep(x, n))
  colnames(Chr7) <- c("Score", "Chr", "Entity", "Sample")
  Chr10 <-data.frame(spata_sub@fdata[[paste0("s_", x)]][["Chr10"]], rep("Chr10", n), rep(Healthy_list[i,3], n), rep(x, n))
  colnames(Chr10) <- c("Score", "Chr", "Entity", "Sample")
  Chr1p <- data.frame(spata_sub@fdata[[paste0("s_", x)]][["Chr1p"]], rep("1p", n), rep(Healthy_list[i,3], n), rep(x, n))
  colnames(Chr1p) <- c("Score", "Chr", "Entity", "Sample")
  Chr1q <- data.frame(spata_sub@fdata[[paste0("s_", x)]][["Chr1q"]], rep("1q", n), rep(Healthy_list[i,3], n), rep(x, n))
  colnames(Chr1q) <- c("Score", "Chr", "Entity", "Sample")
  Chr19p <- data.frame(spata_sub@fdata[[paste0("s_", x)]][["Chr19p"]], rep("19p", n), rep(Healthy_list[i,3], n), rep(x, n))
  colnames(Chr19p) <- c("Score", "Chr", "Entity", "Sample")
  Chr19q <- data.frame(spata_sub@fdata[[paste0("s_", x)]][["Chr19q"]], rep("19q", n), rep(Healthy_list[i,3], n), rep(x, n))
  colnames(Chr19q) <- c("Score", "Chr", "Entity", "Sample")
  
  ifelse(nrow(CNV_all)==0, 
         CNV_all <- rbind(Chr7, Chr10, Chr1p, Chr1q, Chr19p, Chr19q), 
         CNV_all <- rbind(CNV_all, Chr7, Chr10, Chr1p, Chr1q, Chr19p, Chr19q)
  )
  
  CNV_med[nrow(CNV_med)+1,] <- c(median(Chr7[,1]), "Chr7", Healthy_list[i,3], x)
  CNV_med[nrow(CNV_med)+1,] <- c(median(Chr10[,1]), "Chr10", Healthy_list[i,3], x)
  CNV_med[nrow(CNV_med)+1,] <- c(median(Chr1p[,1]), "1p", Healthy_list[i,3], x)
  CNV_med[nrow(CNV_med)+1,] <- c(median(Chr1q[,1]), "1q", Healthy_list[i,3], x)
  CNV_med[nrow(CNV_med)+1,] <- c(median(Chr19q[,1]), "19q", Healthy_list[i,3], x)
  CNV_med[nrow(CNV_med)+1,] <- c(median(Chr19p[,1]), "19p", Healthy_list[i,3], x)
  
}


colnames(CNV_med) <- c("Score", "Chr", "Entity", "Sample")
CNV_med <- na.omit(CNV_med)

#Reorder: colors=pal_npg("nrc")(5)
colors=c("#00A087FF", "#4DBBD5FF", "#F39B7FFF", "#E64B35FF", "#3C5488FF")

order=c("Healthy", "Astrocytoma", "Midline Glioma", "Glioblastoma", "Oligodendroglioma")

#plot median
#p1 <- ggplot(CNV_med, aes(x=factor(Chr, level=c("1p","1q","19p","19q", "Chr7", "Chr10")), y=as.numeric(CNV_med$Score), color = factor(Entity, level=order)))+
#  labs(x="Chromosome Arm", y="CNV Score", main="1p Deletion")+theme_classic()+
#  geom_jitter(position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.1), size = 0.5)+scale_color_manual(values=colors)+
#  geom_boxplot(outlier.shape = NA, alpha=0.5)+
#  labs(colour="Entity")

CNV_OG <- CNV_all[which(!(CNV_all$Chr %in% c("Chr7", "Chr10"))),]
CNV_GBM <- CNV_all[which(CNV_all$Chr %in% c("Chr7", "Chr10")),]

#plot each spot on its own
p2 <- ggplot(CNV_OG, aes(x=factor(Chr, level=c("1p","1q","19p","19q")), y=as.numeric(CNV_OG$Score), color = factor(Entity, level=order)))+
  labs(x="Chromosome Arm", y="CNV Score", main="1p Deletion")+theme_classic()+
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.1), size = 0.5)+scale_color_manual(values=colors)+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  labs(colour="Entity")
ggsave(filename = paste0("Comparison_1p_19q_across_STX.pdf"), plot = p2, height = 6, width = 9)


p3 <- ggplot(CNV_GBM, aes(x=factor(Chr, level=c("Chr7", "Chr10")), y=as.numeric(CNV_GBM$Score), color = factor(Entity, level=order)))+
  labs(x="Chromosome Arm", y="CNV Score", main="1p Deletion")+theme_classic()+
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.1), size = 0.5)+scale_color_manual(values=colors)+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  labs(colour="Entity")
ggsave(filename = paste0("Comparison_7_10_across_STX.pdf"), plot = p3, height = 6, width = 9)


#Run statistical testing


stat.test <- CNV_OG %>%
  group_by(Chr) %>%
  t_test(Score ~ Entity) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test

CNV_med[,1] <- as.numeric(CNV_med[,1])

stat.test1 <- CNV_med %>%
  group_by(Chr) %>%
  t_test(Score ~ Entity) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test1

