library(Seurat)
library(ggplot2)
library(dplyr)
library(rstatix)

setwd("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")
seurat_path <- "/omics/groups/OE0146/internal/Christina/Visium_seq/analysis/GBM/"

#sample_name <- "GB_STX6"

ID_list=data.frame(c("GB_STX7", "GB_STX1", "GB_STX6", "GB_STX4","GB_STX2", "N20_1260_IB", "GB3", 
                     "195096", "GB2", "GB6", "GB11", "GB13",
                     "N20_2382", "GB1", "GB10", "GB12", "GB_STX3", "N21_2972", 
                     "GB9", "GB5", "GB14", "GB15", "GB16", "GB17_B", "GB17_A", "GB18", "GB19",
                     "GB20", "GB21", "GB22"
                     ))


list_of_dataframes=lapply(ID_list[,1],                           
                          function (x){
                            ifelse(x=="GB_STX3" | x=="N21_2972", sample_name <- "GB_STX3_2972", sample_name <- x)
                            print(x)
                            seurat_obj <- readRDS(paste0(seurat_path, sample_name,"/", sample_name, "_seurat_obj.rds"))
                            DefaultAssay(seurat_obj) <- "Spatial"

seurat_obj <- NormalizeData(seurat_obj)
counts <- seurat_obj@assays[["Spatial"]]@data
counts <- data.frame(counts[which(row.names(counts)=="EGFR"),])


Anno <- read.csv(paste0(x,"/Annotated.csv"))
Annot <- Anno[which(Anno$Annotated == "Main Tumour"),]
if(nrow(Annot)==0){print(paste0("No tumour found in sample ", x, ": ", unique(Anno$Annotated)))}
EGFR <- counts[which(row.names(counts) %in% Annot$Barcode),]
fin <- data.frame(cbind(EGFR, rep(x, length(EGFR)), rep(x, length(EGFR))))
colnames(fin) <- c("EGFR_count","Type","sample")

Ctrl_Anno <- Anno[which(Anno$Annotated == "Healthy"),]
Ctrl <- counts[which(row.names(counts) %in% Ctrl_Anno$Barcode),]
Ctrl_fin <- data.frame(cbind(Ctrl, rep("Healthy_Ctrl", length(Ctrl)), rep(x, length(Ctrl))))
colnames(Ctrl_fin) <- c("EGFR_count","Type", "sample")
fin <- rbind(fin, Ctrl_fin)
fin
}
)

fin <- dplyr::bind_rows(list_of_dataframes)
fin$EGFR_count <- as.numeric(fin$EGFR_count)
EGFR_status <- read.csv(file = "/omics/groups/OE0146/internal/Micha/Analysis/GBM/Batch_Analysis/EGFR/EGFR_Status.csv", sep=";")
colnames(EGFR_status) <- c("Type", "EGFR_Amplification", "tissue_type")
EGFR_status[nrow(EGFR_status),] <- c("Healthy_Ctrl", "Healthy_Ctrl", "unknown")
fin <- merge(fin, EGFR_status, by="Type")
fin <- fin[order(fin$EGFR_count),]

ggplot(fin)+theme_classic()+labs(x="Spot", y="LogNorm EGFR count", main="EGFR status")+
  geom_point(mapping = aes(x=1:nrow(fin),y=EGFR_count, color=sample))

ggplot(fin)+theme_classic()+labs(x="Spot", y="LogNorm EGFR count", main="EGFR status")+
  geom_point(mapping = aes(x=1:nrow(fin),y=EGFR_count, color=EGFR_Amplification))

p2=ggplot(fin, aes(x=EGFR_Amplification, y=EGFR_count))+
  labs(x="EGFR Amplification in 850K data", y="LogNorm EGFR count", main="EGFR status")+theme_classic()+
  geom_jitter(size=0.01, alpha=0.5)+
  geom_violin(aes(alpha=0.5, fill=EGFR_Amplification))+
  theme(plot.title=element_text(hjust=0.5))
  #scale_colour_gradient2(low= col[1], mid=col[5], high=col[9], midpoint=1.0, name="Combined Chr 7 and 10 score")+
  #labs(fill="Annotation", title = "Merged Chr. 7 and 10 scores")+
  #geom_boxplot(outlier.shape = NA, alpha=0.3, width=0.1)+
 p2 

fin_med <- data.frame(matrix(ncol=3))
 
for(i in 1:length(unique(fin$sample))){
  x <- unique(fin$sample)[i]
  med_Healthy  <- mean(fin[which(fin$sample==x & fin$EGFR_Amplification == "Healthy_Ctrl"),2])
  med_Tum <- mean(fin[which(fin$sample==x & fin$EGFR_Amplification != "Healthy_Ctrl"),2])
  
  if (!(is.na(med_Healthy))){
    fin_med <- rbind(fin_med, c(x, med_Healthy, "Healthy_Ctrl"))
  }
  
  if (!(is.na(med_Tum))){
  fin_med <- rbind(fin_med, c(x, med_Tum, fin[which(fin$Type==x),][1,4]))
  }
}

colnames(fin_med) <- c("sample", "median_EGFR_expression", "EGFR_Amplification")
fin_med$median_EGFR_expression <- as.numeric(fin_med$median_EGFR_expression)
fin_med <- na.omit(fin_med)

p3=ggplot(fin_med, aes(x=EGFR_Amplification, y=median_EGFR_expression))+
  labs(x="EGFR Amplification in 850K data", y="median EGFR expression", main="EGFR status")+theme_classic()+
  geom_jitter()+
  geom_boxplot(aes(fill=fin_med$EGFR_Amplification), outlier.shape = NA, alpha=0.3)+
  theme(plot.title=element_text(hjust=0.5))
#scale_colour_gradient2(low= col[1], mid=col[5], high=col[9], midpoint=1.0, name="Combined Chr 7 and 10 score")+
#labs(fill="Annotation", title = "Merged Chr. 7 and 10 scores")+
#geom_boxplot(outlier.shape = NA, alpha=0.3, width=0.1)+
p3 

shapiro.test(fin_med[which(fin_med$EGFR_Amplification==F),2])
shapiro.test(fin_med[which(fin_med$EGFR_Amplification==T),2])
shapiro.test(fin_med[which(fin_med$EGFR_Amplification=="Healthy_Ctrl"),2])

stat.test <- fin_med %>%
  #group_by(EGFR_Amplification) %>%
  t_test(median_EGFR_expression ~ EGFR_Amplification) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

#### Using Spata objects

spata_path <- "/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/"

list_of_dataframes=lapply(ID_list[,1],                           
                          function (x){
                            ifelse(x=="GB_STX3" | x=="N21_2972", sample_name <- "GB_STX3_2972", sample_name <- x)
                            print(x)
                            spata_obj <- updateSpataObject(readRDS(paste0(spata_path, x,"/SPATA_obj_", x, ".RDS")))
                            counts <- data.frame(spata_obj@data[[paste0("s_", sample_name)]][["scaled"]])
                            counts <- counts[which(row.names(counts)=="EGFR"),]
                            counts <- data.frame(t(counts))
                            counts[,1] <-(counts+abs(min(counts)))/(max(counts)+abs(min(counts)))
                            Anno <- read.csv(paste0(x,"/Annotated.csv"))
                            Annot <- Anno[which(Anno$Annotated == "Main Tumour"),]
                            Annot$Barcode <- gsub("-", ".", Annot$Barcode)
                            if(nrow(Annot)==0){print(paste0("No tumour found in sample ", x, ": ", unique(Anno$Annotated)))}
                            EGFR <- counts[which(row.names(counts) %in% Annot$Barcode),]
                            fin <- data.frame(cbind(EGFR, rep(x, length(EGFR))))
                            colnames(fin) <- c("EGFR_count", "sample")
                            
                            Ctrl_Anno <- Anno[which(Anno$Annotated == "Healthy"),]
                            Ctrl <- counts[which(row.names(counts) %in% Ctrl_Anno$Barcode),]
                            Ctrl_fin <- data.frame(cbind(Ctrl, rep("Healthy_Ctrl", length(Ctrl))))
                            colnames(Ctrl_fin) <- c("EGFR_count", "sample")
                            fin <- rbind(fin, Ctrl_fin)
                            
                            fin
                          }
)

fin <- dplyr::bind_rows(list_of_dataframes)
fin$EGFR_count <- as.numeric(fin$EGFR_count)
EGFR_status <- read.csv(file = "EGFR_Status.csv", sep=";")
colnames(EGFR_status) <- c("sample", "EGFR_Amplification", "tissue_type")
EGFR_status[nrow(EGFR_status),] <- c("Healthy_Ctrl", "Healthy_Ctrl", "unknown")
fin <- merge(fin, EGFR_status, by="sample")
fin <- fin[order(fin$EGFR_count),]
fin$EGFR_count <- as.numeric(fin$EGFR_count)

ggplot(fin)+theme_classic()+labs(x="Spot", y="LogNorm EGFR count", main="EGFR status")+
  geom_point(mapping = aes(x=1:nrow(fin),y=EGFR_count, color=sample))

ggplot(fin)+theme_classic()+labs(x="Spot", y="LogNorm EGFR count", main="EGFR status")+
  geom_point(mapping = aes(x=1:nrow(fin),y=EGFR_count, color=EGFR_Amplification))

p2=ggplot(fin, aes(x=EGFR_Amplification, y=EGFR_count))+
  labs(x="EGFR Amplification in 850K data", y="LogNorm EGFR count", main="EGFR status")+theme_classic()+
  geom_jitter(size=0.01, alpha=0.5)+
  geom_violin(aes(alpha=0.5, fill=EGFR_Amplification))+
  theme(plot.title=element_text(hjust=0.5))
#scale_colour_gradient2(low= col[1], mid=col[5], high=col[9], midpoint=1.0, name="Combined Chr 7 and 10 score")+
#labs(fill="Annotation", title = "Merged Chr. 7 and 10 scores")+
#geom_boxplot(outlier.shape = NA, alpha=0.3, width=0.1)+
p2 

fin_med <- data.frame(matrix(ncol=3))

for(i in 1:length(unique(fin$sample))){
  x <- unique(fin$sample)[i]
  med  <- median(fin[which(fin$sample==x),2])
  fin_med[i,] <- c(x, med, fin[which(fin$sample==x),][1,3])
}

colnames(fin_med) <- c("sample", "median_EGFR_expression", "EGFR_Amplification")
fin_med$median_EGFR_expression <- as.numeric(fin_med$median_EGFR_expression)

p3=ggplot(fin_med, aes(x=EGFR_Amplification, y=median_EGFR_expression))+
  labs(x="EGFR Amplification in 850K data", y="median EGFR expression", main="EGFR status")+theme_classic()+
  geom_jitter()+
  geom_boxplot(aes(fill=fin_med$EGFR_Amplification), outlier.shape = NA, alpha=0.3)+
  theme(plot.title=element_text(hjust=0.5))
#scale_colour_gradient2(low= col[1], mid=col[5], high=col[9], midpoint=1.0, name="Combined Chr 7 and 10 score")+
#labs(fill="Annotation", title = "Merged Chr. 7 and 10 scores")+
#geom_boxplot(outlier.shape = NA, alpha=0.3, width=0.1)+
p3 




###Figures Paper



x <-"GB_STX7" 
x <-"GB_STX6"

seurat_obj <- readRDS(paste0(seurat_path, x,"/", x, "_seurat_obj.rds"))
DefaultAssay(seurat_obj) <- "Spatial"
seurat_obj <- NormalizeData(seurat_obj)
spata_obj <- updateSpataObject(readRDS(paste0(spata_path, x,"/SPATA_obj_", x, ".RDS")))
counts <- seurat_obj@assays[["Spatial"]]@data
counts <- data.frame(counts[which(row.names(counts)=="EGFR"),])
spata_obj@fdata[[paste0("s_", x)]][["LogNorm_EGFR"]] <- counts[,1]

Frag <- read.csv(paste0(x,"/EGFR.csv"))
spata_obj <- subsetByBarcodes(spata_obj, Frag$Barcode, verbose = NULL)

plot1 <-
  #plotSurface(spata_obj, color_by="LogNorm_EGFR", pt_size=8)+
  plotSurface(spata_obj, color_by="LogNorm_EGFR", pt_size=5)+
  ggpLayerThemeCoords() + scale_color_gradientn(colours = viridis::inferno(5), limits=c(0,5))+
  ggpLayerAxesSI(test, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
  ggtitle("EGFR")
plot1
ggsave(filename = paste0(x,"_EGFR_Expression.pdf"), plot = plot1, path=paste0(x,"/"), height = 5, width = 7)
saveRDS(plot1$data, file=paste0(x,"/", x, "_EGFR_Expression.RDS"))
