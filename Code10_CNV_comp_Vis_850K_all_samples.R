library(ggplot2)
library(dplyr)
library(infercnv)
library(SPATA2)

setwd("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/visium_Analysis/GBM/")
path850K <- ("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")
pathCNV_Chr <- ("/omics/groups/OE0146/internal/Christina/Visium_seq/analysis/GBM/")
pathCNV <- ("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/yiheng/p16/infercnv/merged_ref/")

ID_list=data.frame(c("GB1", "GB2", "GB3", "195096", "GB6", "N20_2382",
                     "GB10", "GB11", "GB12", "GB13", "GB5",
                     "GB_STX1", "GB_STX2", "GB_STX3", "N21_2972", "GB_STX4", "GB_STX6", "GB_STX7"))
ID_list_cryo=data.frame(c("GB15", "GB17_A", "GB18", "GB19",
                     "GB20", "GB22"))


#Not included for being astrocytoma
#"AC_NOS_STX1", "AC_STX6"

#Not included because of low UMI count -> noisy CNV calling:
#"GB7", "AC", "OG", "AC_OG", "MG", "GB14", 

#Not included: Mainly Healthy or no CNVs detected
# "GB9", "N21_1260_IB", "GB17_B", "GB16", "GB21"



list_alias <- read.csv(file="/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/alias.txt")
colnames(list_alias) <- c("ID", "Alias")


list_of_dataframes = lapply(ID_list[, 1],
                            function (x) {
                              #Load 850K data
                              print(x)
                              ifelse(x=="GB_STX3" | x=="N21_2972", sample_name <- "GB_STX3_2972", sample_name <- x)
                              CNV_850K <-
                                read.csv(
                                  paste0(path850K, x, "/", x, "_CNV.segments.seg"),
                                  fill = TRUE,
                                  sep = "\t" ,
                                  stringsAsFactors = FALSE
                                )
                              #Make data.frame of 850K data containing information about length, chormosome, scores, ...
                              CNV_850K <-
                                data.frame(
                                  cbind(
                                    CNV_850K$chrom,
                                    as.numeric(CNV_850K$loc.start),
                                    as.numeric(CNV_850K$loc.end),
                                    as.numeric(CNV_850K$loc.end) - as.numeric(CNV_850K$loc.start),
                                    as.numeric(CNV_850K$seg.mean),
                                    as.numeric(CNV_850K$seg.median)
                                  )
                                )
                              
                              #Rename
                              CNV_850K <-
                                cbind(CNV_850K, as.numeric(sub("chr", "", CNV_850K[, 1])))
                              colnames(CNV_850K) <-
                                c("chrom", "start", "end", "size", "mean", "median", "chr")
                              
                              #get annotation of spots with high tumour content (usually only around 50-100)
                              Tumour <-
                                read.csv(
                                  paste0(path850K, x, "/CNV.csv"),
                                  fill = TRUE,
                                  sep = "," ,
                                  stringsAsFactors = FALSE
                                )
                              #Rename for downstream analysis
                              Tumour[, 1] <- gsub("-", ".", Tumour[, 1])
                              
                              #import infercnv object
                              spata_obj <- updateSpataObject(readRDS(paste0(path850K, x, "/SPATA_obj_", x, ".RDS")))
                              x1 <- readRDS(paste0(pathCNV, sample_name, "/22_denoiseHMMi6.rand_trees.NF_NA.SD_1.5.NL_FALSE.infercnv_obj"))
                              #Extract CNV data for all genes
                              all_chr <-
                                cbind(data.frame(x1@gene_order), data.frame(x1@expr.data))
                              #Extract tumour spots
                              all_chr <-
                                cbind(all_chr[, 1:4], all_chr[, which(colnames(all_chr) %in% Tumour[, 1])])
                              
                              CNV_scores <- data.frame()
                              Gene_number <- data.frame()
                              
                              #Get 850K segements, by cycling through 850K file and using chromosome+start+end information to get all 
                              #infercnv scores
                              for (i in 1:nrow(CNV_850K)) {
                                #select genes on chromosome within stop and start range
                                #Extract genes on a certain chromsosome
                                Genes <- all_chr[which(all_chr$chr == CNV_850K[i, 1]), ]
                                #Select the start and stop position of each segment
                                Genes <- Genes[which(Genes$start >= as.numeric(CNV_850K[i, 2])), ]
                                Genes <- Genes[which(Genes$stop <= as.numeric(CNV_850K[i, 3])), ]
                                #calculate average
                                Genes[, 5:ncol(Genes)] <- lapply(Genes[, 5:ncol(Genes)], as.numeric)
                                Genes[nrow(Genes) + 1, 5:ncol(Genes)] = colMeans(Genes[, 5:ncol(Genes)])
                                CNV_scores <- rbind(CNV_scores, Genes[nrow(Genes), 5:ncol(Genes)])
                                Gene_number <- rbind(Gene_number, nrow(Genes)-1)
                              }
                              #Make Final data.frame
                              Average_score = rowMeans(CNV_scores)
                              merged = cbind(CNV_850K, CNV_scores)
                              fin = cbind(CNV_850K, Average_score, Gene_number, rep(x, nrow(CNV_850K)))
                              colnames(fin) <- c("Chrom_850K", "start", "stop", "size", "850K_mean", "850K_median", "Chrom_Vis", "Vis_mean", "Number_of_Genes", "ID")
                              fin
                            })




fin <- dplyr::bind_rows(list_of_dataframes)

#Run again for cryo samples, they used the cryo reference for CNV calling.

list_of_dataframes_cryo = lapply(ID_list_cryo[, 1],
                            function (x) {
                              #Load 850K data
                              print(x)
                              ifelse(x=="GB_STX3" | x=="N21_2972", sample_name <- "GB_STX3_2972", sample_name <- x)
                              CNV_850K <-
                                read.csv(
                                  paste0(path850K, x, "/", x, "_CNV.segments.seg"),
                                  fill = TRUE,
                                  sep = "\t" ,
                                  stringsAsFactors = FALSE
                                )
                              #Make data.frame of 850K data containing information about length, chormosome, scores, ...
                              CNV_850K <-
                                data.frame(
                                  cbind(
                                    CNV_850K$chrom,
                                    as.numeric(CNV_850K$loc.start),
                                    as.numeric(CNV_850K$loc.end),
                                    as.numeric(CNV_850K$loc.end) - as.numeric(CNV_850K$loc.start),
                                    as.numeric(CNV_850K$seg.mean),
                                    as.numeric(CNV_850K$seg.median)
                                  )
                                )
                              
                              #Rename
                              CNV_850K <-
                                cbind(CNV_850K, as.numeric(sub("chr", "", CNV_850K[, 1])))
                              colnames(CNV_850K) <-
                                c("chrom", "start", "end", "size", "mean", "median", "chr")
                              
                              #get annotation of spots with high tumour content (usually only around 50-100)
                              Tumour <-
                                read.csv(
                                  paste0(path850K, x, "/CNV.csv"),
                                  fill = TRUE,
                                  sep = "," ,
                                  stringsAsFactors = FALSE
                                )
                              #Rename for downstream analysis
                              Tumour[, 1] <- gsub("-", ".", Tumour[, 1])
                              
                              #import infercnv object
                              spata_obj <- updateSpataObject(readRDS(paste0(path850K, x, "/SPATA_obj_", x, ".RDS")))
                              x1 <- readRDS(paste0(pathCNV_Chr, sample_name, "/CNV_SPATA2/infercnv-obj.RDS"))
                              #Extract CNV data for all genes
                              all_chr <-
                                cbind(data.frame(x1@gene_order), data.frame(x1@expr.data))
                              #Extract tumour spots
                              all_chr <-
                                cbind(all_chr[, 1:4], all_chr[, which(colnames(all_chr) %in% Tumour[, 1])])
                              
                              CNV_scores <- data.frame()
                              Gene_number <- data.frame()
                              
                              #Get 850K segements, by cycling through 850K file and using chromosome+start+end information to get all 
                              #infercnv scores
                              for (i in 1:nrow(CNV_850K)) {
                                #select genes on chromosome within stop and start range
                                #Extract genes on a certain chromsosome
                                Genes <- all_chr[which(all_chr$chr == unlist(strsplit(CNV_850K[i, 1], "chr"))[[2]]), ]
                                #Select the start and stop position of each segment
                                Genes <- Genes[which(Genes$start >= as.numeric(CNV_850K[i, 2])), ]
                                Genes <- Genes[which(Genes$stop <= as.numeric(CNV_850K[i, 3])), ]
                                #calculate average
                                Genes[, 5:ncol(Genes)] <- lapply(Genes[, 5:ncol(Genes)], as.numeric)
                                Genes[nrow(Genes) + 1, 5:ncol(Genes)] = colMeans(Genes[, 5:ncol(Genes)])
                                CNV_scores <- rbind(CNV_scores, Genes[nrow(Genes), 5:ncol(Genes)])
                                Gene_number <- rbind(Gene_number, nrow(Genes)-1)
                              }
                              #Make Final data.frame
                              Average_score = rowMeans(CNV_scores)
                              merged = cbind(CNV_850K, CNV_scores)
                              fin = cbind(CNV_850K, Average_score, Gene_number, rep(x, nrow(CNV_850K)))
                              colnames(fin) <- c("Chrom_850K", "start", "stop", "size", "850K_mean", "850K_median", "Chrom_Vis", "Vis_mean", "Number_of_Genes", "ID")
                              fin
                            })


fin_cryo <- dplyr::bind_rows(list_of_dataframes_cryo)

#fin <- rbind(fin, fin_cryo)

#remove rows with NAs
fin <- na.omit(fin)
fin_cryo <- na.omit(fin_cryo)

#save raw file
fin_raw <- fin
fin_cryo_raw <- fin_cryo


#remove segements with less than 5 genes
#print number of rows removed
nrow(fin[which(fin$Number_of_Genes < 5),])
fin <- fin[which(fin$Number_of_Genes >= 5),]
fin <- merge(fin, list_alias, by= "ID")

nrow(fin_cryo[which(fin_cryo$Number_of_Genes < 5),])
fin_cryo <- fin_cryo[which(fin_cryo$Number_of_Genes >= 5),]
fin_cryo <- merge(fin_cryo, list_alias, by= "ID")

#add colours:
#merged[, ncol(merged)+1]=ifelse(as.numeric(merged$size) <=500000, "500.000",ifelse(as.numeric(merged$size) <=5000000,"5.000.000", ifelse(as.numeric(merged$size) <=50000000,"50.000.000", ">50.000.000")))
fin[, ncol(fin)+1]=ifelse(as.numeric(fin$size) <=500000, "500.000",ifelse(as.numeric(fin$size) <=5000000,"5.000.000", ifelse(as.numeric(fin$size) <=50000000,"50.000.000", ">50.000.000")))
fin[, ncol(fin)+1]=ifelse(as.numeric(fin$Number_of_Genes) <5, "<5",ifelse(as.numeric(fin$Number_of_Genes) <=10,"5-10", ifelse(as.numeric(fin$Number_of_Genes) <=25,"10-25", ">25")))
colnames(fin)[12:13] <- c("Segment_Size", "Genes_per_seg_colour")

fin_raw[, ncol(fin_raw)+1]=ifelse(as.numeric(fin_raw$size) <=500000, "500.000",ifelse(as.numeric(fin_raw$size) <=5000000,"5.000.000", ifelse(as.numeric(fin_raw$size) <=50000000,"50.000.000", ">50.000.000")))
fin_raw[, ncol(fin_raw)+1]=ifelse(as.numeric(fin_raw$Number_of_Genes) <5, "<5",ifelse(as.numeric(fin_raw$Number_of_Genes) <=10,"5-10", ifelse(as.numeric(fin_raw$Number_of_Genes) <=25,"10-25", ">25")))
colnames(fin_raw)[11:12] <- c("Segment_Size", "Genes_per_seg_colour")

fin_cryo[, ncol(fin_cryo)+1]=ifelse(as.numeric(fin_cryo$size) <=500000, "500.000",ifelse(as.numeric(fin_cryo$size) <=5000000,"5.000.000", ifelse(as.numeric(fin_cryo$size) <=50000000,"50.000.000", ">50.000.000")))
fin_cryo[, ncol(fin_cryo)+1]=ifelse(as.numeric(fin_cryo$Number_of_Genes) <5, "<5",ifelse(as.numeric(fin_cryo$Number_of_Genes) <=10,"5-10", ifelse(as.numeric(fin_cryo$Number_of_Genes) <=25,"10-25", ">25")))
colnames(fin_cryo)[12:13] <- c("Segment_Size", "Genes_per_seg_colour")

fin_cryo_raw[, ncol(fin_cryo_raw)+1]=ifelse(as.numeric(fin_cryo_raw$size) <=500000, "500.000",ifelse(as.numeric(fin_cryo_raw$size) <=5000000,"5.000.000", ifelse(as.numeric(fin_cryo_raw$size) <=50000000,"50.000.000", ">50.000.000")))
fin_cryo_raw[, ncol(fin_cryo_raw)+1]=ifelse(as.numeric(fin_cryo_raw$Number_of_Genes) <5, "<5",ifelse(as.numeric(fin_cryo_raw$Number_of_Genes) <=10,"5-10", ifelse(as.numeric(fin_cryo_raw$Number_of_Genes) <=25,"10-25", ">25")))
colnames(fin_cryo_raw)[11:12] <- c("Segment_Size", "Genes_per_seg_colour")

#Claculate standard deviation
No_CNV <- fin[which((as.numeric(fin[,6]) >= -0.075) & as.numeric(fin[,6]) <= 0.075),]
lim <- sd(fin[,9])/2

No_CNV_cryo <- fin_cryo[which(fin_cryo[,6] >= -0.075 & fin_cryo[,6] <= 0.075),]
lim_cryo <- sd(fin_cryo[,9])/2

#Basic functions to plot:
#ggplot(merged, aes(x=as.numeric(mean), y=as.numeric(merged[,8]), color= merged[,ncol(merged)]))+geom_point()
#p1 <- ggplot(fin, aes(x=as.numeric(fin$'850K_mean')+1, y=as.numeric(fin$'Vis_mean')))+geom_point()+labs(x="CNV Scores 850K", y="CNV Scores Visium")+theme_bw()

#Colour according to sample, shape according to segemnt size
p1 <- ggplot(fin, aes(x=as.numeric(fin$'850K_mean')+1, y=as.numeric(fin$'Vis_mean'), color=Alias))+
  geom_point(aes(shape=fin$'Segment_Size'))+labs(x="CNV Scores 850K", y="CNV Scores Visium")+
  theme_bw()+
  geom_hline(yintercept=c(1+lim,1-lim), linetype="dotted")+
  geom_vline(xintercept=c(1.075, 0.925), linetype="dotted")+
  labs(color="Genes per Segment", shape ="Segment Size", title = "CNV Comparison Median")+
  ylab("Visium CNV scores") + theme(plot.title=element_text(hjust=0.5))+
  xlim(0.25, 1.75)
p1

#Colour according to genes per segment
p2 <- ggplot(fin_raw, aes(x=as.numeric(fin_raw$'850K_mean')+1, y=as.numeric(fin_raw$'Vis_mean'), color=Genes_per_seg_colour))+
  geom_point(size=1)+scale_color_manual(values = c("red", "blue", "green3", "orange"))+
  theme_bw()+
  geom_hline(yintercept=c(1+lim,1-lim), linetype="dotted")+
  geom_vline(xintercept=c(1.075, 0.925), linetype="dotted")+
  labs(color="Genes per Segment", shape ="Segment Size", title = "CNV Comparison Median")+
  ylab("Visium CNV scores") + theme(plot.title=element_text(hjust=0.5))+
  xlim(0.0, 2.0)
#p1 <- ggplot(fin, aes(x=as.numeric(fin$'850K_mean')+1, y=as.numeric(fin$'Vis_mean'), color=Genes_per_seg_colour))+geom_point(size=1)+scale_color_manual(values = c("blue", "green3", "orange"))+theme_bw()
p2

p1 <- ggplot(fin_cryo, aes(x=as.numeric(fin_cryo$'850K_mean')+1, y=as.numeric(fin_cryo$'Vis_mean'), color=Alias))+
  geom_point(aes(shape=fin_cryo$'Segment_Size'))+labs(x="CNV Scores 850K", y="CNV Scores Visium")+
  theme_bw()+
  geom_hline(yintercept=c(1+lim_cryo,1-lim_cryo), linetype="dotted")+
  geom_vline(xintercept=c(1.075, 0.925), linetype="dotted")+
  labs(color="Genes per Segment", shape ="Segment Size", title = "CNV Comparison Median")+
  ylab("Visium CNV scores") + theme(plot.title=element_text(hjust=0.5))+
  xlim(0.25, 1.75)
p1

#Colour according to genes per segment
p2 <- ggplot(fin_cryo_raw, aes(x=as.numeric(fin_cryo_raw$'850K_mean')+1, y=as.numeric(fin_cryo_raw$'Vis_mean'), color=Genes_per_seg_colour))+
  geom_point(size=1)+scale_color_manual(values = c("red", "blue", "green3", "orange"))+
  theme_bw()+
  geom_hline(yintercept=c(1+lim_cryo,1-lim_cryo), linetype="dotted")+
  geom_vline(xintercept=c(1.075, 0.925), linetype="dotted")+
  labs(color="Genes per Segment", shape ="Segment Size", title = "CNV Comparison Median")+
  ylab("Visium CNV scores") + theme(plot.title=element_text(hjust=0.5))+
  xlim(0.0, 2.0)
p2


# Calculating specifity and sensitivity
#All correct predicitons
CNV_T <- fin[which((as.numeric(fin[,9]) <= 1-lim & as.numeric(fin[,6]) <= -0.075) | (as.numeric(fin[,9]) > 1-lim & as.numeric(fin[,9]) < 1+lim & as.numeric(fin[,6]) > -0.075 & as.numeric(fin[,6]) < 0.075) | (as.numeric(fin[,9]) >= 1+lim & as.numeric(fin[,6]) >= 0.075)),]
#All wrong predictions
CNV_F <- fin[which((as.numeric(fin[,9]) > 1-lim & as.numeric(fin[,6]) <= -0.075) | (as.numeric(fin[,9]) < 1-lim  & as.numeric(fin[,6]) > -0.075 & as.numeric(fin[,6]) < 0.075) | (as.numeric(fin[,9]) > 1+lim  & as.numeric(fin[,6]) > -0.075 & as.numeric(fin[,6]) < 0.075) | (as.numeric(fin[,9]) <= 1+lim & as.numeric(fin[,6]) >= 0.075)),]
#How many predicitons were correct:
print(nrow(CNV_T)/(nrow(CNV_T)+nrow(CNV_F))*100)
#True positives
CNV_1 <- fin[which((as.numeric(fin[,9]) <= 1-lim & as.numeric(fin[,6]) <= -0.075) | (as.numeric(fin[,9]) >= 1+lim & as.numeric(fin[,6]) >= 0.075)),]
#True negatives
CNV_2 <- fin[which((as.numeric(fin[,9]) > 1-lim & as.numeric(fin[,9]) < 1+lim & as.numeric(fin[,6]) > -0.075 & as.numeric(fin[,6]) < 0.075)),]
#False positives
CNV_3 <- fin[which((as.numeric(fin[,9]) > 1-lim & as.numeric(fin[,6]) <= -0.075) | (as.numeric(fin[,9]) <= 1+lim & as.numeric(fin[,6]) >= 0.075)),]
#False negatives
CNV_4 <- fin[which((as.numeric(fin[,9]) < 1-lim  & as.numeric(fin[,6]) > -0.075 & as.numeric(fin[,6]) < 0.075) | (as.numeric(fin[,9]) > 1+lim  & as.numeric(fin[,6]) > -0.075 & as.numeric(fin[,6]) < 0.075)),]
#Specificity:
print(nrow(CNV_2)/(nrow(CNV_2)+nrow(CNV_3))*100)
#Sensitivity:
print(nrow(CNV_1)/(nrow(CNV_1)+nrow(CNV_3))*100)

#CNV_TV <- fin[which(as.numeric(fin[,9]) <= 0.95 | as.numeric(fin[,9]) >= 1.05),]
#CNV_FV <- fin[which(as.numeric(fin[,9]) > 0.95 & as.numeric(fin[,9]) < 1.05),]
#CNV_T8 <- fin[which(as.numeric(fin[,6]) <= -0.1 | as.numeric(fin[,6]) >= 0.1),]
#CNV_F8 <- fin[which(as.numeric(fin[,6]) > -0.1 & as.numeric(fin[,6]) < 0.1),]

# Calculating specifity and sensitivity
#All correct predicitons
CNV_T <- fin_cryo[which((as.numeric(fin_cryo[,9]) <= 1-lim_cryo & as.numeric(fin_cryo[,6]) <= -0.075) | (as.numeric(fin_cryo[,9]) > 1-lim_cryo & as.numeric(fin_cryo[,9]) < 1+lim_cryo & as.numeric(fin_cryo[,6]) > -0.075 & as.numeric(fin_cryo[,6]) < 0.075) | (as.numeric(fin_cryo[,9]) >= 1+lim_cryo & as.numeric(fin_cryo[,6]) >= 0.075)),]
#All wrong predictions
CNV_F <- fin_cryo[which((as.numeric(fin_cryo[,9]) > 1-lim_cryo & as.numeric(fin_cryo[,6]) <= -0.075) | (as.numeric(fin_cryo[,9]) < 1-lim_cryo  & as.numeric(fin_cryo[,6]) > -0.075 & as.numeric(fin_cryo[,6]) < 0.075) | (as.numeric(fin_cryo[,9]) > 1+lim_cryo  & as.numeric(fin_cryo[,6]) > -0.075 & as.numeric(fin_cryo[,6]) < 0.075) | (as.numeric(fin_cryo[,9]) <= 1+lim_cryo & as.numeric(fin_cryo[,6]) >= 0.075)),]
#How many predicitons were correct:
print(nrow(CNV_T)/(nrow(CNV_T)+nrow(CNV_F))*100)
#True positives
CNV_1 <- fin_cryo[which((as.numeric(fin_cryo[,9]) <= 1-lim_cryo & as.numeric(fin_cryo[,6]) <= -0.075) | (as.numeric(fin_cryo[,9]) >= 1+lim_cryo & as.numeric(fin_cryo[,6]) >= 0.075)),]
#True negatives
CNV_2 <- fin_cryo[which((as.numeric(fin_cryo[,9]) > 1-lim_cryo & as.numeric(fin_cryo[,9]) < 1+lim_cryo & as.numeric(fin_cryo[,6]) > -0.075 & as.numeric(fin_cryo[,6]) < 0.075)),]
#False positives
CNV_3 <- fin_cryo[which((as.numeric(fin_cryo[,9]) > 1-lim_cryo & as.numeric(fin_cryo[,6]) <= -0.075) | (as.numeric(fin_cryo[,9]) <= 1+lim_cryo & as.numeric(fin_cryo[,6]) >= 0.075)),]
#False negatives
CNV_4 <- fin_cryo[which((as.numeric(fin_cryo[,9]) < 1-lim_cryo  & as.numeric(fin_cryo[,6]) > -0.075 & as.numeric(fin_cryo[,6]) < 0.075) | (as.numeric(fin_cryo[,9]) > 1+lim_cryo  & as.numeric(fin_cryo[,6]) > -0.075 & as.numeric(fin_cryo[,6]) < 0.075)),]
#Specificity:
print(nrow(CNV_2)/(nrow(CNV_2)+nrow(CNV_3))*100)
#Sensitivity:
print(nrow(CNV_1)/(nrow(CNV_1)+nrow(CNV_3))*100)

#CNV_TV <- fin_cryo[which(as.numeric(fin_cryo[,9]) <= 0.99 | as.numeric(fin_cryo[,9]) >= 1.01),]
#CNV_FV <- fin_cryo[which(as.numeric(fin_cryo[,9]) > 0.99 & as.numeric(fin_cryo[,9]) < 1.01),]
#CNV_T8 <- fin_cryo[which(as.numeric(fin_cryo[,6]) <= -0.1 | as.numeric(fin_cryo[,6]) >= 0.1),]
#CNV_F8 <- fin_cryo[which(as.numeric(fin_cryo[,6]) > -0.1 & as.numeric(fin_cryo[,6]) < 0.1),]



