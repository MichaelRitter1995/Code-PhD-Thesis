library(ggplot2)
library(dplyr)
library(reshape)
library(rstatix)

setwd("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")
ID_list=read.table("ID_list.txt", header = T)

list_of_dataframes=lapply(ID_list[,1],
                          function (x){ 
                            Annotation <- read.csv(paste0(x,"/Annotated.csv"), fill=TRUE)
                            Annotation[,ncol(Annotation)+1]=rep(x, nrow(Annotation)) 
                            colnames(Annotation)=c("barcodes", "Annotation", "Sample ID")
                            
                            Subtypes <- t(read.csv(paste0(x,"/",x,"_predictions_Neftel.csv"), fill=TRUE, header=FALSE, stringsAsFactors = FALSE))
                            colnames(Subtypes)=c("barcodes", Subtypes[1,2:ncol(Subtypes)])
                            Subtypes <- Subtypes[2:nrow(Subtypes),]
                            Subtypes_merged <- merge(Annotation, Subtypes, by="barcodes", all.x=TRUE)
                          })

Subtypes_all=dplyr::bind_rows(list_of_dataframes)
Subtypes_all[which(Subtypes_all$Annotation == "Necrose"),2]<- "Necrosis"
Subtypes_all[,ncol(Subtypes_all)+1]<-ifelse(Subtypes_all$Annotation %in% c("IFZ", "Necrosis", "Main Tumour", "Healthy"), Subtypes_all$Annotation, "Other")

Subtypes_per_sample <- data.frame()

for(i in 1:length(unique(Subtypes_all$`Sample ID`))){
  samp <- unique(Subtypes_all$`Sample ID`)[[i]]
  subset <- Subtypes_all[which(Subtypes_all$`Sample ID`==samp),]
  for(j in 1:length(unique(subset$Annotation))){
    st <- unique(subset$Annotation)[[j]]
    Hist <- subset[which(subset$Annotation == st),]
    Hist <- na.omit(Hist)
    if(nrow(Hist)>1){
      Mean_hist <- colMeans(sapply(Hist[,4:17], as.numeric))
      Mean_fin <-cbind(Hist[1,2:3], t(Mean_hist), Hist[1,c(19)])
      colnames(Mean_fin)[[17]] <- "V19"
    } else {
      Mean_fin <-Hist[1,c(2:17, 19)]
      colnames(Mean_fin) <- colnames(Subtypes_per_sample)
    }
    Subtypes_per_sample <- rbind(Subtypes_per_sample, Mean_fin)
    }
}

Subtypes_per_sample_only <- melt(Subtypes_per_sample, id.vars='V19' , measure.vars=c("NPC1","OPC","AC","NPC2", "MES1","MES2") )


t1 <- ggplot(Subtypes_per_sample_only)+theme_bw()+labs(x="Mapping Score", y= NA, main="Mapping of GBM Subtypes")
t2 <- geom_boxplot(aes(x=V19, y=as.numeric(value), color=variable), outlier.shape = NA)
t3 <- geom_jitter(aes(x=V19, y=as.numeric(value), color=variable), position = position_jitterdodge())

t1+t3+t2+theme_classic()

#Only analyse necrotic spots

Subtypes_only <- melt(Subtypes_all, id.vars=c("barcodes",'V19') , measure.vars=c("NPC1","OPC","AC","NPC2", "MES1","MES2"))
#Necr=Subtypes_only[which(Subtypes_only[,2]=="Necrosis" & Subtypes_only$variable %in% c("OPC", "AC", "MES1", "MES2")),]
Necr=Subtypes_only[which(Subtypes_only[,2]=="Necrosis"),]


t4 <- ggplot(Necr)+theme_bw()+labs(y="Mapping Score", x= "Necrosis", main="Mapping of GBM Subtypes")+
  geom_jitter(aes(x=variable, y=as.numeric(value), color=variable), alpha=0.8)+
  geom_violin(aes(x=variable, y=as.numeric(value), color=variable))+
  theme_bw()
t4

NecrS <- Necr[which(Necr$variable %in% c("MES1", "MES2", "AC", "OPC")),]
pairwise.wilcox.test(as.numeric(NecrS$value), NecrS$variable, p.adjust.method="bonferroni")

Necr1 <- Necr[which(Necr$variable %in% c("MES1", "AC")),]
Necr1$value <- as.numeric(Necr1$value)

ggwithinstats( # paired samples
data = Necr1,
x = variable,
y = value,
type = "nonparametric",
centrality.plotting = FALSE
)
