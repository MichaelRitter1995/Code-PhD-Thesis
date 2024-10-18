setwd("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/visium_Analysis/GBM/Merged/")

library(SPATA2)
library(dplyr)
library(Seurat)
library(ggplot2)
library(rstatix)

path <-  "/omics/groups/OE0146/internal/Christina/Visium_seq/analysis/GBM/"
pathCNVs <- "/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/yiheng/p16/infercnv/merged_ref/"
sample_list = list("195096","GB2","GB6","N20_2382","GB12", "GB9", "GB1", "N20_1260_IB", "GB3", "GB5", "GB7", "GB11",
                   "AC_NOS_STX1","GB_STX1","AC_STX6","GB_STX3_2972", "GB_STX4","GB13", "GB_STX2", "GB_STX6", "GB_STX7", "GB10")


Cell_Types<- list(
  c("MBP","OLIG2","OLIG3","CLDN11","MOG","SOX10","GJB1","MAG","CNP","PLP1"),
  c("ENO2","NEFL","NEFM","NEFH","INA","SLC17A7","RBFOX3","SNAP25","SYT1","NEUROD2", "MAP2", "RBFOX1", "SNAP25")
)

fin = data.frame()

for (i in 1:length(sample_list)) {
  sample_name = sample_list[[i]]
  print(sample_name)
  
  #load objects
  spata_obj <-
    updateSpataObject(readRDS(
      paste0(path, sample_name, "/", sample_name, "_SPATA_obj.rds")
    ))
  seurat_obj <-
    readRDS(paste0(path, sample_name, "/", sample_name, "_seurat_obj.rds"))
  
  #Add module scores to seurat object
  #counts= data.frame(seurat_obj@data[[sample_name]][["counts"]])
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = Cell_Types,
    name = c('Cell_Types'),
    assay = 'Spatial'
  )
  
  OF <- select(seurat_obj@meta.data, Cell_Types1)
  NF <- select(seurat_obj@meta.data, Cell_Types2)
  
  #Extract CNVs and calcualte 7/10 score
  cnv <- spata_obj@cnv[[sample_name]]$cnv_df
  cnv[, ncol(cnv) + 1] <- (cnv$chr7 - 1 + 1 - cnv$chr10) / 2 + 1
  stscore <- cbind(cnv$barcodes, cnv[, ncol(cnv)])
  
  #clean up dataframe
  colnames(stscore) <- c("barcodes", "7/10_Score")
  row.names(stscore) <- stscore[, 1]
  stscore <- data.frame(stscore[order(stscore[, 1]), ])
  
  #generate final dataframe
  stscore <- stscore[which(stscore[,1] %in% row.names(NF)),]
  OF <- OF[which(row.names(OF) %in% row.names(stscore)), ]
  NF <- NF[which(row.names(NF) %in% row.names(stscore)), ]
  
  merged <- cbind(stscore, OF, NF)
  
  
  #fin=rbind(temp[2,], MBP[1,1:ncol(MBP)-1], NEFL[1,1:ncol(NEFL)-1])
  fin <- rbind(fin, merged)
}

#Plotting V1

#plot(as.numeric(fin[1,]), as.numeric(fin[2,]), cex=0.2)
#plot(as.numeric(fin[1,]), as.numeric(fin[3,]), cex=0.2)

#p1 <-ggplot(fin)+theme_bw()
#p2 <- geom_point(aes(x= fin$X7.10_Score, y=fin[,3:4]), size=0.5, alpha=0.5, color="blue", stroke = 0)
#p3 <- geom_point(aes(x= fin$X7.10_Score, y=fin$NF), size=0.5, alpha=0.5, color="red", stroke = 0)
#l1 <- labs(color="Comparison Visium 850K", y="Marker Expression", x="+7/-10 Score")
#c1 <- scale_color_manual(breaks=c("Oligodendrocytes", "Neurons"), colour=c("blue", "red"))
#p1+p2+l1

#Plotting V2

fin$Olig <- "Oligo"
fin$Neuron <- "Neuron"

#generate dataframe
d1 <- select(fin, X7.10_Score, OF, Olig)
d2 <- select(fin, X7.10_Score, NF, Neuron)

#Normalize
d1$OF <- d1$OF/max(d1$OF)
d2$NF <- d2$NF/max(d2$NF)


###Approach driven by zeros
colnames(d2) <- colnames(d1)
stupid_dataframe <- bind_rows(d1, d2)

cutoff <- sd(d1$X7.10_Score)

p1 <- ggplot(stupid_dataframe,
             aes(
               x = X7.10_Score,
               y = OF,
               color = Olig
             )) + geom_hline(aes(yintercept = 0.00), color="lightgrey") + geom_vline(aes(xintercept = 1.00), color="lightgrey") +
  geom_point(size=1, alpha=0.5, stroke = 0) +
  scale_color_manual(values=c("#339900", "#CC6600"))+
  labs(
    x = "+7/-10 Score",
    y = "Module Score",
    color = "Cell Type"
  ) + theme_classic() + guides(colour = guide_legend(override.aes = list(size=3))) + ggtitle("Comparison Marker Expression +7/-10 Score") + ylim(c(-0.4,1))

p1+geom_vline(xintercept = 1+cutoff/2)
ggsave(plot = p1, filename = "Marker_Expression_7_10_Score.pdf", width = 24, height = 8)

#Test for significance
stupid_dataframe <- stupid_dataframe %>%
  mutate(category = ifelse(X7.10_Score < 1+cutoff/2, "Healthy", "Tumour"))

stupid_dataframe <- stupid_dataframe %>%
  mutate(Testing = ifelse(category == "Healthy" & Olig == "Neuron", "NH", 
                          ifelse(category == "Tumour" & Olig == "Neuron", "NT",
                                 ifelse(category == "Healthy" & Olig == "Oligo", "OH", "OT")
                                 )
                          )
         )

#Test for a normal distirbution
NH <- stupid_dataframe[which(stupid_dataframe$Testing=="NH"),]
NT <- stupid_dataframe[which(stupid_dataframe$Testing=="NT"),]
OT <- stupid_dataframe[which(stupid_dataframe$Testing=="OH"),]
OH <- stupid_dataframe[which(stupid_dataframe$Testing=="OT"),]
qqnorm(NH$OF)
qqline(NH$OF)
qqnorm(NT$OF)
qqline(NT$OF)
qqnorm(OH$OF)
qqline(OH$OF)
qqnorm(OT$OF)
qqline(OT$OF)



#Dunne's test
kruskal.test(OF~Testing, data=stupid_dataframe)
dunn_test(OF~Testing, data=stupid_dataframe, p.adjust.method = "bonferroni")

p11 <- ggplot(data=stupid_dataframe, aes(x=Testing, y=OF))+geom_violin()+geom_boxplot(aes(alpha=0))+theme_classic()
ggsave(plot = p11, filename = "Vln_Marker_Expression_7_10_Score.pdf", width = 24, height = 8)


median(NH$OF)
median(NT$OF)
median(OH$OF)
median(OT$OF)

#generate dataframe
d1 <- select(fin, X7.10_Score, OF, Olig)
d2 <- select(fin, X7.10_Score, NF, Neuron)

#Normalize
d1$OF <- d1$OF/max(d1$OF)
d2$NF <- d2$NF/max(d2$NF)
new_anno <- cbind(d1, d2)
new_anno <- new_anno[,c(1,2,5)] 
new_anno <- new_anno %>% mutate(
    Score = pmax(NF, OF),
    Annotation = ifelse(NF >= OF, "Cortex", "White Matter")
  )

p2 <- ggplot(new_anno,
             aes(
               x = X7.10_Score,
               y = Score,
               color = Annotation
             )) + geom_hline(aes(yintercept = 0.00), color="lightgrey") + geom_vline(aes(xintercept = 1.00), color="lightgrey") +
  geom_point(size=1, alpha=0.5, stroke = 0) +
  scale_color_manual(values=c("#339900", "#CC6600"))+
  labs(
    x = "+7/-10 Score",
    y = "Module Score",
    color = "Cell Type"
  ) + theme_classic() + guides(colour = guide_legend(override.aes = list(size=3))) + ggtitle("Comparison Marker Expression +7/-10 Score") + ylim(c(-0.4,1))
p2
ggsave(plot = p2, filename = "New_Annotation_Marker_Expression_7_10_Score.pdf", width = 24, height = 8)

p3 <-ggplot(data=new_anno, aes(x=Annotation, y=X7.10_Score))+geom_violin()+geom_boxplot(aes(alpha=0))+theme_classic()
ggsave(plot = p3, filename = "Vln_New_Annotation_Marker_Expression_7_10_Score.pdf", width = 24, height = 8)

#Test
wilcox.test(Score~Annotation, data=new_anno)

new_anno <- new_anno %>%
  mutate(category = ifelse(X7.10_Score < 1+cutoff/2, "Healthy", "Tumour"))

new_anno <- new_anno %>%
  mutate(Testing = ifelse(category == "Healthy" & Annotation == "Cortex", "NH", 
                          ifelse(category == "Tumour" & Annotation == "Cortex", "NT",
                                 ifelse(category == "Healthy", "OH", "OT")
                          )
  )
  )

NH <- nrow(new_anno[which(new_anno$Testing=="NH"),])
NT <- nrow(new_anno[which(new_anno$Testing=="NT"),])
OT <- nrow(new_anno[which(new_anno$Testing=="OT"),])
OH <- nrow(new_anno[which(new_anno$Testing=="OH"),])

testing <- cbind(c(NH, NT), c(OH, OT))
colnames(testing) <- c("Cortex", "White Matter")
row.names(testing) <- c("Healthy", "Tumour")
chisq.test(testing)
