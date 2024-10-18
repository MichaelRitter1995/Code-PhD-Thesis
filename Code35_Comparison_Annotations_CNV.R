library(ggplot2)
library(dplyr)
library(SPATA2)
library(RColorBrewer)
library(rstatix)

setwd("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")
pathVis <- ("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")

col=c(RColorBrewer::brewer.pal(9, "BrBG"))

ID_list=read.table("ID_list.txt")
ID_list=data.frame(c("GB1", "GB2", "GB3", "195096", "N20_1260_IB", "GB5", "GB6", "N20_2382",
                     "GB9", "GB10", "GB11", "GB12", "GB13",
                     "GB_STX1", "GB_STX2", "GB_STX3", "N21_2972", "GB_STX4", "GB_STX6", "GB_STX7",
                     "GB14", "GB15", "GB16", "GB17_B", "GB17_A", "GB18", "GB19",
                     "GB20", "GB21", "GB22"))
                    


#Not included for being astrocytoma -> clean cohort:
#"AC_NOS_STX1", "AC_STX6"
#Not included because of low UMI count -> noisy CNV calling:
#"GB7", "AC", "OG", "AC_OG", "MG"

list_of_dataframes=lapply(ID_list[,1],                           
                          function (x){
                            ifelse(x=="GB_STX3" | x=="N21_2972", sample_name <- "GB_STX3_2972", sample_name <- x)
                            print(x)
                            Annotation <- read.csv(paste0(x,"/Annotated.csv"), fill=TRUE)
                            Annotation[,ncol(Annotation)+1]=rep(x, nrow(Annotation)) 
                            colnames(Annotation)=c("barcodes", "Annotation", "Sample ID")
                            
                            spata_obj <- updateSpataObject(readRDS(paste0(pathVis, x, "/SPATA_obj_", x, ".RDS")))
                            CNVs <- cbind(spata_obj@fdata[[paste0("s_", sample_name)]][["barcodes"]], spata_obj@fdata[[paste0("s_", sample_name)]][["Chr7"]], spata_obj@fdata[[paste0("s_", sample_name)]][["Chr10"]])
                            colnames(CNVs)=c("barcodes", "Chr7", "Chr10")
                            CNV_merged <- merge(Annotation, CNVs, by="barcodes", all.x=TRUE)
                          })


CNV_all=dplyr::bind_rows(list_of_dataframes)

CNV_all=cbind(CNV_all, 1+(as.numeric(CNV_all$Chr7)-as.numeric(CNV_all$Chr10))/2)
colnames(CNV_all)[6]="7/10 Score"
CNV_all[,7]<-ifelse(CNV_all$Annotation %in% c("IFZ", "Necrosis", "Main Tumour", "Healthy"), CNV_all$Annotation, "Other")
#CNV_all[,30]<-ifelse(CNV_all$Annotation %in% c("Healthy"), "Mainly Healthy", CNV_all$V30)
CNV_all[,7] <- factor(CNV_all[,7] , levels=c("Healthy", "IFZ", "Main Tumour", "Necrosis", "Other"))

Necr=CNV_all[which(CNV_all[,7] =="Necrosis"),6]
MH=CNV_all[which(CNV_all[,7] =="Healthy"),6]
TM=CNV_all[which(CNV_all[,7] =="Main Tumour"),6]
IFZ=CNV_all[which(CNV_all[,7] =="IFZ"),6]
Other=CNV_all[which(CNV_all[,7] =="Other"),6]

CNV_median <- data.frame()

for(i in 1:nrow(ID_list)){
  x <- ID_list[i,1]
  CNV_samp <- CNV_all[which(CNV_all$`Sample ID` == x),]
  Annotations <- c("Necrosis", "Healthy", "Main Tumour", "IFZ", "Other")
  for(j in 1:length(Annotations)){
    Feat <- Annotations[j]
    Feat_samp=CNV_samp[which(CNV_samp[,7] ==Feat),]
    if(nrow(Feat_samp) != 0){
      Feat_samp[,4:5] <- lapply(Feat_samp[,4:5], as.numeric)
      Feat_samp <- na.omit( Feat_samp)
      CNV_median <- rbind(CNV_median, c(Feat, x, median(Feat_samp$Chr7), median(Feat_samp$Chr10),(1+(median(Feat_samp$Chr7)-median(Feat_samp$Chr10))/2)))
      }
    }
  }

colnames(CNV_median) <- c("Annotation", "Sample_ID", "Chr7", "Chr10", "7/10Score")

CNV_median$`7/10Score` <- as.numeric(CNV_median$`7/10Score`)

colors=c("IFZ"=brewer.pal(10, "Set3")[5], "Main Tumour"=brewer.pal(10, "Set3")[4], "Healthy"=brewer.pal(10, "Set3")[7], "Necrosis"="black", "Other"="grey")


p1=ggplot(CNV_all, aes(x=CNV_all[,7], y=CNV_all[,6]))+labs(x="Annotation", y="Combined Chr. 7 and 10 Alteration", main="Chr7 and 10 Changes according to Annotation")+theme_bw()+
   geom_jitter(aes(color=CNV_all[,6]), alpha=0.5, size=0.01)+
#t3=geom_boxplot(aes(color=CNV_all[,30]))+
   geom_boxplot(fill=colors, outlier.shape = NA, alpha=0.5)+
   scale_colour_gradient2(low= col[1], mid=col[5], high=col[9], midpoint=1.0, name="Combined Chr 7 and 10 score")+
   labs(fill="Annotation", title = "Merged Chr. 7 and 10 scores")+
   theme(plot.title=element_text(hjust=0.5))
p1

p2=ggplot(CNV_all, aes(x=CNV_all[,7], y=CNV_all[,6]))+labs(x="Annotation", y="Combined Chr. 7 and 10 Alteration", main="Chr7 and 10 Changes according to Annotation")+theme_bw()+
  geom_jitter(aes(color=CNV_all[,6]), size=0.01)+
  #t3=geom_boxplot(aes(color=CNV_all[,30]))+
  geom_violin(aes(fill=CNV_all[,7]), alpha=0.5)+
  scale_fill_manual(values=colors)+
  scale_colour_gradient2(low= col[1], mid=col[5], high=col[9], midpoint=1.0, name="Combined Chr 7 and 10 score")+
  labs(fill="Annotation", title = "Merged Chr. 7 and 10 scores")+
  geom_boxplot(outlier.shape = NA, alpha=0.3, width=0.1)+
  theme(plot.title=element_text(hjust=0.5))
p2

p3=ggplot(CNV_median, aes(x=Annotation, y=`7/10Score`))+labs(x="Annotation", y="Combined Median chr. 7 and 10 Alteration", main="Median Chr7 and 10 Changes according to Annotation")+theme_bw()+
  geom_jitter(aes(color=`7/10Score`), size=1)+
  #t3=geom_boxplot(aes(color=CNV_all[,30]))+
  geom_boxplot(aes(fill=Annotation), alpha=0.5, width=0.5)+
  scale_fill_manual(values=colors)+
  scale_colour_gradient2(low= col[9], mid="white", high=col[1], midpoint=(min(CNV_median$`7/10Score`)+(max(CNV_median$`7/10Score`)-min(CNV_median$`7/10Score`))/2), name="Combined Chr 7 and 10 score")+
  #scale_colour_gradientn(colours = viridis::inferno(n = 100))+
  labs(fill="Annotation", title = "Merged Chr. 7 and 10 scores")+
  theme(plot.title=element_text(hjust=0.5))
p3

sd_chr7_healty <- CNV_all %>% filter(Annotation=="Healthy") %>% pull(Chr7) %>% sd()
median_chr7__healty <- CNV_all %>% filter(Annotation=="Healthy") %>% pull(Chr7) %>% median()

sd_chr10_healty <- CNV_all %>% filter(Annotation=="Healthy") %>% pull(Chr10) %>% sd()
median_chr10__healty <- CNV_all %>% filter(Annotation=="Healthy") %>% pull(Chr10) %>% median()
  
sd_merged_healty <- CNV_all %>% filter(Annotation=="Healthy") %>% pull("7/10 Score") %>% sd()
median_merged__healty <- CNV_all %>% filter(Annotation=="Healthy") %>% pull("7/10 Score") %>% median()

CNV_all <- CNV_all %>% arrange(CNV_all$`7/10 Score`)

p4 <- ggplot(CNV_all %>% arrange(CNV_all$`7/10 Score`))+
  geom_point(mapping = aes(x=1:nrow(CNV_all), y=as.numeric(CNV_all$Chr10)), alpha=0.1, size=0.1, color=scales::muted("blue"))+
             geom_point(mapping = aes(x=1:nrow(CNV_all), y=as.numeric(CNV_all$Chr7)), alpha=0.1,size=0.1, color=scales::muted("red"))+
               geom_point(mapping = aes(x=1:nrow(CNV_all), y=CNV_all$`7/10 Score`, alpha="7/10 Score"), size=0.8, color=scales::muted("pink"))+
               geom_hline(yintercept = median_merged__healty, color="black")+
               geom_hline(yintercept = median_merged__healty+sd_merged_healty, color="black", linetype="dashed")+
                          geom_hline(yintercept = median_merged__healty+(2*sd_merged_healty), color="black", linetype="dashed")+
                                     theme_classic()+
                                       xlab("Index Spots")+
                                       ylab("Inferred CNV Score")+
                                       Seurat::NoLegend()
p4

p5 <- ggplot(CNV_all %>% arrange(CNV_all$`7/10 Score`))+
  geom_point(mapping = aes(x=1:nrow(CNV_all), y=CNV_all$`7/10 Score`, alpha="7/10 Score"), size=0.8, color=scales::muted("pink"))+
  geom_hline(yintercept = median_merged__healty, color="black")+
  geom_hline(yintercept = median_merged__healty+sd_merged_healty, color="black", linetype="dashed")+
  geom_hline(yintercept = median_merged__healty+(2*sd_merged_healty), color="black", linetype="dashed")+
  theme_classic()+
  xlab("Index Spots")+
  ylab("Inferred CNV Score")+
  Seurat::NoLegend()
p5

CNV_all$Annotation <- 
CNV_all[2,!(which(CNV_all[2,] %in% c("IFZ", "Main Tumour", "Healthy")))]<-"Other"

shapiro.test(CNV_median[which(CNV_median$Annotation=="Healthy"),]$`7/10Score`)
shapiro.test(CNV_median[which(CNV_median$Annotation=="Healthy"),]$`7/10Score`)
shapiro.test(CNV_median[which(CNV_median$Annotation=="Healthy"),]$`7/10Score`)

CNV_median1 <- CNV_median[which(CNV_median$Annotation %in% c("IFZ", "Main Tumour", "Healthy")),]

stat.test <- CNV_median1 %>%
  t_test(`7/10Score` ~ Annotation) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
