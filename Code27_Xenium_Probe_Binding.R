setwd("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/Combined_Analysis/visium/")

library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)

path="/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/"

####Prepare probeset####

vis.probes <- read.csv("/software/spaceranger/2.0.1/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv", sep="|")
vis.probes <- matrix(vis.probes[6:nrow(vis.probes),1], ncol=3, byrow = T)

#Extract all probes for all genes used (including control probes)
Overlapping.probes <- read.csv(paste0(path,"reference_seq/costumized_probes/Xenium_CytAssist_probes/Pre_designed_panel/Visium_Xenium_Alignment/All_probes.csv"), sep="|")
Overlapping.probes <- matrix(Overlapping.probes[6:nrow(Overlapping.probes),1], ncol=3, byrow = T)
all.probes <- vis.probes[which(vis.probes[,2] %in% Overlapping.probes[,2]),]
#Extract all probes of genes with overlapping probes
Overlapping.probes <- read.csv(paste0(path,"reference_seq/costumized_probes/Xenium_CytAssist_probes/Pre_designed_panel/Visium_Xenium_Alignment/All_Overlapping_probes.csv"), sep="|", header = F)
ov.probes <- vis.probes[which(vis.probes[,2] %in% Overlapping.probes[,2]),]

#Extracting the total counts of all genes
list_raw_counts <- lapply(c("4405", "4137"), function(x, probeset=all.probes)
{
  Seu_obj1.data <- Read10X(data.dir = paste0(path, "/visium/aligned/Xenium/GB/",x,"/outs/filtered_feature_bc_matrix/"))
  Seu_obj1 <- CreateSeuratObject(counts = Seu_obj1.data, min.cells = 0, min.features = 0)
  
  Counts <- data.frame(Seu_obj1@assays[["RNA"]]@counts[which(rownames(Seu_obj1@assays[["RNA"]]@counts) %in% probeset[,2]),])
  Counts <- cbind(row.names(Counts), rowSums(Counts), rep(x, nrow(Counts)), rep("all_probes", nrow(Counts)))
  colnames(Counts) <- c("Gene", "Counts", "Sample", "Probe_set")
  probe.number <- data.frame(table(all.probes[,2]))
  colnames(probe.number) <- c("Gene", "Probe_Number")
  Counts <- merge(Counts, probe.number, by.x="Gene")
  
}
)

raw_counts <- bind_rows(list_raw_counts, .id = "sample")


####Analyse Probe distribution#######

#Create Reference
Ref_probes <- ov.probes
Ref_probes <- data.frame(table(Ref_probes[,2]))
colnames(Ref_probes) <- c("Gene", "Total_probes")

All_probes <- all.probes
All_probes <- data.frame(table(All_probes[,2]))
colnames(Ref_probes) <- c("Gene", "Total_probes")

PE3 <- read.csv(paste0(path, "reference_seq/costumized_probes/Xenium_CytAssist_probes/Pre_designed_panel/Visium_Xenium_Alignment/5_prime_end_seq.txt"), sep="|")
PE3 <- cbind(matrix(unlist(strsplit(PE3[,1], ","), recursive = FALSE), ncol=3, byrow=T)[,2], PE3[,2])

PE5 <-read.csv(paste0(path, "reference_seq/costumized_probes/Xenium_CytAssist_probes/Pre_designed_panel/Visium_Xenium_Alignment/3_prime_end_seq.txt"), sep="|")
PE5 <- cbind(matrix(unlist(strsplit(PE5[,1], ","), recursive = FALSE), ncol=3, byrow=T)[,2], PE5[,2])



fin <- data.frame()

for(i in c(0,5:20)){
  #Extract all genes of interest form the probe set and also store number of probes in a data.frame
  probes <- read.csv(paste0(path, "reference_seq/costumized_probes/Xenium_CytAssist_probes/Pre_designed_panel/Visium_Xenium_Alignment/Revised_Probe_set_",i,"_BP_Binding.csv"), sep="|")
  probes <- data.frame(matrix(probes[6:nrow(probes),1], ncol=3, byrow = T))
  #Save sequence for GC content
  sequence <- data.frame(matrix(unlist(strsplit(probes[,1], ","), recursive = FALSE), ncol=3, byrow=T))
  sequence <- data.frame(cbind(probes[,2], sequence[,2]))
  
  #Procede with probe clean-up
  probes <- data.frame(table(probes[,2]))
  if(i !=0){probes <- probes[which(probes[,1] %in% Ref_probes[,1]),]}
  probes <- merge(probes, All_probes[which(All_probes[,1] %in% probes[,1]),], by="Var1", all.y = T)
  probes[,ncol(probes)+1] <- (probes[,3]-probes[,2])
  colnames(probes) <- c("Gene", "Probes_number", "Total_probes", "Difference")
  probes[,ncol(probes)+1] <- round((probes$Probes_number)/(probes$Total_probes), digits=2)*100
  
  #Analyse CG content
  
  sequence <- sequence[which(sequence[,1] %in% probes[which(probes$V5=="33"),1]),]
  colnames(sequence) <- c("Gene","Sequence")
  
  if(i !=0 & i !=20 ){
    sequencePE3 <- sequence[which(sequence[,2] %in% PE3[,1]),]
    sequencePE3[,3] <- substring(sequencePE3[,2],26,25+i)
    sequencePE3[,4] <- (str_count(sequencePE3[,3], "G")+str_count(sequencePE3[,3], "C"))/(i)*100
    sequencePE5 <- sequence[which(sequence[,2] %in% PE5[,1]),]
    sequencePE5[,3] <- substring(sequencePE5[,2],26-i,25)
    sequencePE5[,4] <- (str_count(sequencePE5[,3], "G")+str_count(sequencePE5[,3], "C"))/(i)*100
    sequence <- rbind(sequencePE3, sequencePE5)
    probes <- merge(probes, sequence, by = "Gene", all.x = T)
    probes <- probes[,c(1:5,8)]
    
  } else if (i ==20) {
    sequence[,3] <- substring(sequence[,2],6,45)
    sequence[,4] <- (str_count(sequence[,3], "G")+str_count(sequence[,3], "C"))/(i*2)*100
    sequence <- rbind(sequence, sequence)
    probes <- merge(probes, sequence, by = "Gene", all.x = T)
    probes <- probes[,c(1:5,8)]
  } else{
    sequence[,3] <- (str_count(sequence[,2], "G")+str_count(sequence[,2], "C"))/50*100
    probes <- merge(probes, sequence, by = "Gene", all.x = T)
    probes <- probes[,c(1:5,7)]
  }
  colnames(probes)[c(5,6)] <- c("Percentage_removed", "GC-content")
  
  #Add 3' and 5' Analysis to it
  
  #if probes are in the list, use the generated list to extract the appropriate reads
  if(nrow(probes)!=0){
    
    list_counts <- lapply(c("4405", "4137"), function(x, probeset=probes, y=i, Ref=raw_counts)
    {
      #Create Seurat object
      Seu_obj1.data <- Read10X(data.dir = paste0(path, "visium/aligned/Xenium/Probe_Binding/Brain_panel/", y,"_BP_Binding/",x,"/outs/filtered_feature_bc_matrix/"))
      Seu_obj1 <- CreateSeuratObject(counts = Seu_obj1.data, min.cells = 0, min.features = 0)
      
      #Extract reads
      Counts1 <- data.frame(Seu_obj1@assays[["RNA"]]@counts[which(rownames(Seu_obj1@assays[["RNA"]]@counts) %in% probes[,1]),])
      Counts1 <- cbind(row.names(Counts1), rowSums(Counts1), rep(x, nrow(Counts1)), rep(paste0(y,"_BP_Binding"), nrow(Counts1)))
      colnames(Counts1) <- c("Gene", "Counts", "Sample", "Probe_set")
      fin_counts <- list(data.frame(Counts1), Ref[which(Ref$Gene %in% Counts1[,1] & Ref$Sample == x),]) %>% reduce(full_join, by="Gene")
      fin_counts <- fin_counts[,c(1:4,6)]
      fin_counts[,c(2,5)] <- apply(fin_counts[,c(2,5)],2, function(x) {as.numeric(x)})
      fin_counts[,ncol(fin_counts)+1] <- (fin_counts$Counts.x)/(fin_counts$Counts.y)*100
      colnames(fin_counts) <- c("Gene", "Counts", "Sample", "Probe_set", "Total_Counts", "Percentage_Lost")
      fin_counts
    }
    )
    fin_counts <- bind_rows(list_counts)
    fin_i <- merge(fin_counts, probes, by="Gene")
    fin <- rbind(fin, fin_i)
  }
}

fin <- cbind(fin, as.numeric(data.frame(str_split(fin[,4], "_"))[1,]))
colnames(fin)[ncol(fin)] <- "BP_Binding"

fin1 <- fin[which(fin$Percentage_removed == 33),]
fin1 <- rbind(fin1, c(rep(NA, 3), "6_BP_Binding", rep(NA, 5), 33, NA, 6))
fin1$Percentage_Lost <- as.numeric(fin1$Percentage_Lost)
fin1$Counts <- as.numeric(fin1$Counts)
fin1$`GC-content` <- as.numeric(fin1$`GC-content`)

p1 <- ggplot(fin1, aes(x=factor(Percentage_removed, level=c(33, 50, 67,100)), y=Percentage_Lost, colour=Probe_set))+theme_classic()+
  geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single"))+    
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.1), size = 0.5)+
  geom_hline(yintercept = 66.66, linetype=2)+geom_hline(yintercept = 33.33, linetype=2)+geom_hline(yintercept = 0, linetype=2)+geom_hline(yintercept = 50, linetype=2)
p1

plot_aes=c("0_BP_Binding", "5_BP_Binding","6_BP_Binding","7_BP_Binding","8_BP_Binding","9_BP_Binding","10_BP_Binding", "11_BP_Binding", "12_BP_Binding","13_BP_Binding","14_BP_Binding","15_BP_Binding", "16_BP_Binding","17_BP_Binding","18_BP_Binding","19_BP_Binding","20_BP_Binding")

p2 <- ggplot(fin1, aes(x=factor(BP_Binding, levels = c(0:20)), y=Percentage_Lost, colour=factor(Probe_set, levels=plot_aes)))+theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  ggtitle("Determining binding of probes (only 1 probe per gene included)")+
  labs(y = "Fraction of total reads", x = "number of BP of Visium probe binding")+
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.1), size = 1)+
  geom_hline(yintercept = 66.66, linetype=2)+geom_hline(yintercept = 33.33, linetype=2)+geom_hline(yintercept = 0, linetype=2)+geom_hline(yintercept = 50, linetype=2)
p2

barplot(sort(fin$Counts))

p3 <- ggplot(fin1, aes(x=factor(BP_Binding, levels = c(0:20)), y=Percentage_Lost, colour=Counts))+theme_classic()+
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.1), size = 1)+
  scale_colour_gradientn(limits = c(0,25000), colours=c("navyblue", "white", "darkred"))+
  geom_hline(yintercept = 66.66, linetype=2)+geom_hline(yintercept = 33.33, linetype=2)+geom_hline(yintercept = 0, linetype=2)+geom_hline(yintercept = 50, linetype=2)
p3

p4 <- ggplot(fin1, aes(x=factor(BP_Binding, levels = c(0:20)), y=Percentage_Lost, colour=fin1$`GC-content`))+theme_classic()+
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.1), size = 2)+
  scale_colour_gradientn(limits=c(30,70), colours=c("navyblue",
    "#ffffbf",
    "darkred"))+
  geom_hline(yintercept = 66.66, linetype=2)+geom_hline(yintercept = 33.33, linetype=2)+geom_hline(yintercept = 0, linetype=2)+geom_hline(yintercept = 50, linetype=2)
p4


