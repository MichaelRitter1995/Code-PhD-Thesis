setwd("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")
#setwd("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Normal/")

library(ggplot2)
library(dplyr)
library(infercnv)
library(SPATA2)

list_alias <- read.csv(file="/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/alias.txt")

path <- "/omics/groups/OE0146/internal/Christina/Visium_seq/analysis/GBM/AC/CNV_SPATA2"
x1 <- readRDS(paste0(path, "/infercnv-obj.RDS"))
ID_list <- data.frame(c("GB107", "AC_STX5", "AC_STX4", "AC_STX3"))
list_of_dataframes = lapply(ID_list[, 1],
                            function (x) {
                              print(x)
                              CNV_850K <-
                                read.csv(
                                  paste0("AC/",x, "/", x, "_CNV.segments.seg"),
                                  fill = TRUE,
                                  sep = "\t" ,
                                  stringsAsFactors = FALSE
                                )
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
                              
                              #Read850K file
                              CNV_850K <-
                                cbind(CNV_850K, as.numeric(sub("chr", "", CNV_850K[, 1])))
                              colnames(CNV_850K) <-
                                c("chrom", "start", "end", "size", "mean", "median", "chr")
                              
                              #get Tumour Spots
                              Tumour <-
                                read.csv(
                                  paste0("AC/",x, "/CNV.csv"),
                                  fill = TRUE,
                                  sep = "," ,
                                  stringsAsFactors = FALSE
                                )
                              Tumour[, 1] <- gsub("-", ".", Tumour[, 1])
                              
                              #Get genes with chromosome information
                              all_chr <-
                                cbind(data.frame(x1@gene_order), data.frame(x1@expr.data))
                              all_chr <-
                                cbind(all_chr[, 1:4], all_chr[, which(colnames(all_chr) %in% Tumour[, 1])])
                              
                              CNV_scores <- data.frame()
                              Gene_number <- data.frame()
                              
                              for (i in 1:nrow(CNV_850K)) {
                                Genes <- all_chr[which(all_chr$chr == CNV_850K[i, 7]), ]
                                Genes <- Genes[which(Genes$start >= as.numeric(CNV_850K[i, 2])), ]
                                Genes <- Genes[which(Genes$stop <= as.numeric(CNV_850K[i, 3])), ]
                                Genes[, 5:ncol(Genes)] <- lapply(Genes[, 5:ncol(Genes)], as.numeric)
                                Genes[nrow(Genes) + 1, 5:ncol(Genes)] = colMeans(Genes[, 5:ncol(Genes)])
                                CNV_scores <- rbind(CNV_scores, Genes[nrow(Genes), 5:ncol(Genes)])
                                Gene_number <- rbind(Gene_number, nrow(Genes)-1)
                              }
                              Average_score = rowMeans(CNV_scores)
                              merged = cbind(CNV_850K, CNV_scores)
                              fin = cbind(CNV_850K, Average_score, Gene_number, rep(x, nrow(CNV_850K)))
                              colnames(fin) <- c("Chrom_850K", "start", "stop", "size", "850K_mean", "850K_median", "Chrom_Vis", "Vis_mean", "Number_of_Genes", "ID")
                              fin
                            })
fin_AC <- dplyr::bind_rows(list_of_dataframes)


path850K <- ("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")
pathCNV_Chr <- ("/omics/groups/OE0146/internal/Christina/Visium_seq/analysis/GBM/")
pathCNV <- ("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/yiheng/p16/infercnv/merged_ref/")
ID_list <- data.frame(c("AC_NOS_STX1", "AC_STX6"))
list_of_dataframes = lapply(ID_list[, 1],
                            function (x) {
                              #Load 850K data
                              print(x)
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
                              
                              #import infercnv object; different possibilities: Yiheng (FFPE ref); Christina (cryo ref); Henrik (850K bins)
                              x1 <- readRDS(paste0("/omics/groups/OE0146/internal/Christina/Visium_seq/analysis/GBM/", x, "/CNV_SPATA2/infercnv-obj.RDS"))
                              
                              #spata_obj <- updateSpataObject(readRDS(paste0(path850K, x, "/SPATA_obj_", x, ".RDS")))
                              #x1 <- readRDS(paste0(pathCNV, sample_name, "/22_denoiseHMMi6.rand_trees.NF_NA.SD_1.5.NL_FALSE.infercnv_obj"))
                              
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
fin_AC1 <- dplyr::bind_rows(list_of_dataframes)


x1 <- readRDS(paste0(path, "/infercnv-obj.RDS"))
ID_list <- data.frame(c("AC_STX2", "AC_STX1", "OG_STX4"))
list_of_dataframes = lapply(ID_list[, 1],
                            function (x) {
                              print(x)
                              CNV_850K <-
                                read.csv(
                                  paste0("AC_OG/",x, "/", x, "_CNV.segments.seg"),
                                  fill = TRUE,
                                  sep = "\t" ,
                                  stringsAsFactors = FALSE
                                )
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
                              CNV_850K <-
                                cbind(CNV_850K, as.numeric(sub("chr", "", CNV_850K[, 1])))
                              colnames(CNV_850K) <-
                                c("chrom", "start", "end", "size", "mean", "median", "chr")
                              
                              Tumour <-
                                read.csv(
                                  paste0("AC_OG/",x, "/CNV.csv"),
                                  fill = TRUE,
                                  sep = "," ,
                                  stringsAsFactors = FALSE
                                )
                              Tumour[, 1] <- gsub("-", ".", Tumour[, 1])
                              
                              all_chr <-
                                cbind(data.frame(x1@gene_order), data.frame(x1@expr.data))
                              all_chr <-
                                cbind(all_chr[, 1:4], all_chr[, which(colnames(all_chr) %in% Tumour[, 1])])
                              
                              CNV_scores <- data.frame()
                              Gene_number <- data.frame()
                              
                              for (i in 1:nrow(CNV_850K)) {
                                Genes <- all_chr[which(all_chr$chr == CNV_850K[i, 7]), ]
                                Genes <- Genes[which(Genes$start >= as.numeric(CNV_850K[i, 2])), ]
                                Genes <- Genes[which(Genes$stop <= as.numeric(CNV_850K[i, 3])), ]
                                Genes[, 5:ncol(Genes)] <- lapply(Genes[, 5:ncol(Genes)], as.numeric)
                                Genes[nrow(Genes) + 1, 5:ncol(Genes)] = colMeans(Genes[, 5:ncol(Genes)])
                                CNV_scores <- rbind(CNV_scores, Genes[nrow(Genes), 5:ncol(Genes)])
                                Gene_number <- rbind(Gene_number, nrow(Genes)-1)
                              }
                              Average_score = rowMeans(CNV_scores)
                              merged = cbind(CNV_850K, CNV_scores)
                              fin = cbind(CNV_850K, Average_score, Gene_number, rep(x, nrow(CNV_850K)))
                              colnames(fin) <- c("Chrom_850K", "start", "stop", "size", "850K_mean", "850K_median", "Chrom_Vis", "Vis_mean", "Number_of_Genes", "ID")
                              fin
                            })
fin_AC_OG <- dplyr::bind_rows(list_of_dataframes)


x1 <- readRDS(paste0(path, "/infercnv-obj.RDS"))
ID_list <- data.frame(c("OG_STX3", "OG_STX1"))
list_of_dataframes = lapply(ID_list[, 1],
                            function (x) {
                              print(x)
                              CNV_850K <-
                                read.csv(
                                  paste0("OG/",x, "/", x, "_CNV.segments.seg"),
                                  fill = TRUE,
                                  sep = "\t" ,
                                  stringsAsFactors = FALSE
                                )
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
                              CNV_850K <-
                                cbind(CNV_850K, as.numeric(sub("chr", "", CNV_850K[, 1])))
                              colnames(CNV_850K) <-
                                c("chrom", "start", "end", "size", "mean", "median", "chr")
                              
                              Tumour <-
                                read.csv(
                                  paste0("OG/",x, "/CNV.csv"),
                                  fill = TRUE,
                                  sep = "," ,
                                  stringsAsFactors = FALSE
                                )
                              Tumour[, 1] <- gsub("-", ".", Tumour[, 1])
                              
                              all_chr <-
                                cbind(data.frame(x1@gene_order), data.frame(x1@expr.data))
                              all_chr <-
                                cbind(all_chr[, 1:4], all_chr[, which(colnames(all_chr) %in% Tumour[, 1])])
                              
                              CNV_scores <- data.frame()
                              Gene_number <- data.frame()
                              
                              for (i in 1:nrow(CNV_850K)) {
                                Genes <- all_chr[which(all_chr$chr == CNV_850K[i, 7]), ]
                                Genes <- Genes[which(Genes$start >= as.numeric(CNV_850K[i, 2])), ]
                                Genes <- Genes[which(Genes$stop <= as.numeric(CNV_850K[i, 3])), ]
                                Genes[, 5] <-  as.numeric(Genes[,5])
                                CNV_scores <- rbind(CNV_scores, colMeans(Genes[, 5, drop = FALSE]))
                                Gene_number <- rbind(Gene_number, nrow(Genes)-1)
                              }
                              Average_score = rowMeans(CNV_scores)
                              merged = cbind(CNV_850K, CNV_scores)
                              fin = cbind(CNV_850K, Average_score, Gene_number, rep(x, nrow(CNV_850K)))
                              colnames(fin) <- c("Chrom_850K", "start", "stop", "size", "850K_mean", "850K_median", "Chrom_Vis", "Vis_mean", "Number_of_Genes", "ID")
                              fin
                            })
fin_OG <- dplyr::bind_rows(list_of_dataframes)




fin <- rbind(fin_AC, fin_AC_OG, fin_OG, fin_AC1)

#remove rows with NAs
fin <- na.omit(fin)
fin_AC1 <- na.omit(fin_AC1)

#save raw file
fin_raw <- fin
fin_AC1_raw <- fin_AC1

#remove segements with less than 5 genes
#print number of rows removed
nrow(fin[which(fin$Number_of_Genes < 5),])
fin <- fin[which(fin$Number_of_Genes >= 5),]

nrow(fin_AC1[which(fin_AC1$Number_of_Genes < 5),])
fin_AC1 <- fin_AC1[which(fin_AC1$Number_of_Genes >= 5),]


#add colours:
#merged[, ncol(merged)+1]=ifelse(as.numeric(merged$size) <=500000, "500.000",ifelse(as.numeric(merged$size) <=5000000,"5.000.000", ifelse(as.numeric(merged$size) <=50000000,"50.000.000", ">50.000.000")))
fin[, ncol(fin)+1]=ifelse(as.numeric(fin$size) <=500000, "500.000",ifelse(as.numeric(fin$size) <=5000000,"5.000.000", ifelse(as.numeric(fin$size) <=50000000,"50.000.000", ">50.000.000")))
fin[, ncol(fin)+1]=ifelse(as.numeric(fin$Number_of_Genes) <5, "<5",ifelse(as.numeric(fin$Number_of_Genes) <=10,"5-10", ifelse(as.numeric(fin$Number_of_Genes) <=25,"10-25", ">25")))
colnames(fin)[11:12] <- c("Segment_Size", "Genes_per_seg_colour")

fin_raw[, ncol(fin_raw)+1]=ifelse(as.numeric(fin_raw$size) <=500000, "500.000",ifelse(as.numeric(fin_raw$size) <=5000000,"5.000.000", ifelse(as.numeric(fin_raw$size) <=50000000,"50.000.000", ">50.000.000")))
fin_raw[, ncol(fin_raw)+1]=ifelse(as.numeric(fin_raw$Number_of_Genes) <5, "<5",ifelse(as.numeric(fin_raw$Number_of_Genes) <=10,"5-10", ifelse(as.numeric(fin_raw$Number_of_Genes) <=25,"10-25", ">25")))
colnames(fin_raw)[11:12] <- c("Segment_Size", "Genes_per_seg_colour")

fin_AC1[, ncol(fin_AC1)+1]=ifelse(as.numeric(fin_AC1$size) <=500000, "500.000",ifelse(as.numeric(fin_AC1$size) <=5000000,"5.000.000", ifelse(as.numeric(fin_AC1$size) <=50000000,"50.000.000", ">50.000.000")))
fin_AC1[, ncol(fin_AC1)+1]=ifelse(as.numeric(fin_AC1$Number_of_Genes) <5, "<5",ifelse(as.numeric(fin_AC1$Number_of_Genes) <=10,"5-10", ifelse(as.numeric(fin_AC1$Number_of_Genes) <=25,"10-25", ">25")))
colnames(fin_AC1)[11:12] <- c("Segment_Size", "Genes_per_seg_colour")

fin_AC1_raw[, ncol(fin_AC1_raw)+1]=ifelse(as.numeric(fin_AC1_raw$size) <=500000, "500.000",ifelse(as.numeric(fin_AC1_raw$size) <=5000000,"5.000.000", ifelse(as.numeric(fin_AC1_raw$size) <=50000000,"50.000.000", ">50.000.000")))
fin_AC1_raw[, ncol(fin_AC1_raw)+1]=ifelse(as.numeric(fin_AC1_raw$Number_of_Genes) <5, "<5",ifelse(as.numeric(fin_AC1_raw$Number_of_Genes) <=10,"5-10", ifelse(as.numeric(fin_AC1_raw$Number_of_Genes) <=25,"10-25", ">25")))
colnames(fin_AC1_raw)[11:12] <- c("Segment_Size", "Genes_per_seg_colour")

#Claculate standard deviation
No_CNV <- fin[which((as.numeric(fin[,6]) >= -0.075) & as.numeric(fin[,6]) <= 0.075),]
lim <- sd(fin[,8])/2


#Basic functions to plot:
#ggplot(merged, aes(x=as.numeric(mean), y=as.numeric(merged[,8]), color= merged[,ncol(merged)]))+geom_point()
#p1 <- ggplot(fin, aes(x=as.numeric(fin$'850K_mean')+1, y=as.numeric(fin$'Vis_mean')))+geom_point()+labs(x="CNV Scores 850K", y="CNV Scores Visium")+theme_bw()

#Colour according to sample, shape according to segemnt size
p1 <- ggplot(fin, aes(x=as.numeric(fin$'850K_mean')+1, y=as.numeric(fin$'Vis_mean'), color=ID))+
  geom_point(aes(shape=fin$'Segment_Size'))+labs(x="CNV Scores 850K", y="CNV Scores Visium")+
  theme_bw()+
  geom_hline(yintercept=c(1+lim,1-lim), linetype="dotted")+
  geom_vline(xintercept=c(1.075, 0.925), linetype="dotted")+
  labs(color="Genes per Segment", shape ="Segment Size", title = "CNV Comparison Median")+
  ylab("Visium CNV scores") + theme(plot.title=element_text(hjust=0.5))
p1

#Colour according to genes per segment
p2 <- ggplot(fin_raw, aes(x=as.numeric(fin_raw$'850K_mean')+1, y=as.numeric(fin_raw$'Vis_mean'), color=Genes_per_seg_colour))+
  geom_point(size=1)+scale_color_manual(values = c("red", "blue", "green3", "orange"))+
  theme_bw()+
  geom_hline(yintercept=c(1+lim,1-lim), linetype="dotted")+
  geom_vline(xintercept=c(1.075, 0.925), linetype="dotted")+
  labs(color="Genes per Segment", shape ="Segment Size", title = "CNV Comparison Median")+
  ylab("Visium CNV scores") + theme(plot.title=element_text(hjust=0.5))
#p1 <- ggplot(fin, aes(x=as.numeric(fin$'850K_mean')+1, y=as.numeric(fin$'Vis_mean'), color=Genes_per_seg_colour))+geom_point(size=1)+scale_color_manual(values = c("blue", "green3", "orange"))+theme_bw()
p2


p1 <- ggplot(fin_AC1, aes(x=as.numeric(fin_AC1$'850K_mean')+1, y=as.numeric(fin_AC1$'Vis_mean'), color=ID))+
  geom_point(aes(shape=fin_AC1$'Segment_Size'))+labs(x="CNV Scores 850K", y="CNV Scores Visium")+
  theme_bw()+
  geom_hline(yintercept=c(1.01,0.99), linetype="dotted")+
  geom_vline(xintercept=c(1.1, 0.9), linetype="dotted")+
  labs(color="Genes per Segment", shape ="Segment Size", title = "CNV Comparison Median")+
  ylab("Visium CNV scores") + theme(plot.title=element_text(hjust=0.5))
p1

#Colour according to genes per segment
p2 <- ggplot(fin_AC1_raw, aes(x=as.numeric(fin_AC1_raw$'850K_mean')+1, y=as.numeric(fin_AC1_raw$'Vis_mean'), color=Genes_per_seg_colour))+
  geom_point(size=1)+scale_color_manual(values = c("red", "blue", "green3", "orange"))+
  theme_bw()+
  geom_hline(yintercept=c(1.01,0.99), linetype="dotted")+
  geom_vline(xintercept=c(1.1, 0.9), linetype="dotted")+
  labs(color="Genes per Segment", shape ="Segment Size", title = "CNV Comparison Median")+
  ylab("Visium CNV scores") + theme(plot.title=element_text(hjust=0.5))
#p1 <- ggplot(fin_AC1, aes(x=as.numeric(fin_AC1$'850K_mean')+1, y=as.numeric(fin_AC1$'Vis_mean'), color=Genes_per_seg_colour))+geom_point(size=1)+scale_color_manual(values = c("blue", "green3", "orange"))+theme_bw()
p2

CNV_T <- fin[which((as.numeric(fin[,8]) <= 1-lim & as.numeric(fin[,6]) <= -0.075) | (as.numeric(fin[,8]) > 1-lim & as.numeric(fin[,8]) < 1+lim & as.numeric(fin[,6]) > -0.075 & as.numeric(fin[,6]) < 0.075) | (as.numeric(fin[,8]) >= 1+lim & as.numeric(fin[,6]) >= 0.075)),]
#All wrong predictions
CNV_F <- fin[which((as.numeric(fin[,8]) > 1-lim & as.numeric(fin[,6]) <= -0.075) | (as.numeric(fin[,8]) < 1-lim  & as.numeric(fin[,6]) > -0.075 & as.numeric(fin[,6]) < 0.075) | (as.numeric(fin[,8]) > 1+lim  & as.numeric(fin[,6]) > -0.075 & as.numeric(fin[,6]) < 0.075) | (as.numeric(fin[,8]) <= 1+lim & as.numeric(fin[,6]) >= 0.075)),]
#How many predicitons were correct:
print(nrow(CNV_T)/(nrow(CNV_T)+nrow(CNV_F))*100)
#True positives
CNV_1 <- fin[which((as.numeric(fin[,8]) <= 1-lim & as.numeric(fin[,6]) <= -0.075) | (as.numeric(fin[,8]) >= 1+lim & as.numeric(fin[,6]) >= 0.075)),]
#True negatives
CNV_2 <- fin[which((as.numeric(fin[,8]) > 1-lim & as.numeric(fin[,8]) < 1+lim & as.numeric(fin[,6]) > -0.075 & as.numeric(fin[,6]) < 0.075)),]
#False positives
CNV_3 <- fin[which((as.numeric(fin[,8]) > 1-lim & as.numeric(fin[,6]) <= -0.075) | (as.numeric(fin[,8]) <= 1+lim & as.numeric(fin[,6]) >= 0.075)),]
#False negatives
CNV_4 <- fin[which((as.numeric(fin[,8]) < 1-lim  & as.numeric(fin[,6]) > -0.075 & as.numeric(fin[,6]) < 0.075) | (as.numeric(fin[,8]) > 1+lim  & as.numeric(fin[,6]) > -0.075 & as.numeric(fin[,6]) < 0.075)),]
#Specificity:
print(nrow(CNV_2)/(nrow(CNV_2)+nrow(CNV_3))*100)
#Sensitivity:
print(nrow(CNV_1)/(nrow(CNV_1)+nrow(CNV_3))*100)



CNV_T <- fin_AC1[which((as.numeric(fin_AC1[,8]) <= 0.99 & as.numeric(fin_AC1[,6]) <= -0.1) | (as.numeric(fin_AC1[,8]) > 0.99 & as.numeric(fin_AC1[,8]) < 1.05 & as.numeric(fin_AC1[,6]) > -0.1 & as.numeric(fin_AC1[,6]) < 0.1) | (as.numeric(fin_AC1[,8]) >= 1.05 & as.numeric(fin_AC1[,6]) >= 0.1)),]
#All wrong predictions
CNV_F <- fin_AC1[which((as.numeric(fin_AC1[,8]) > 0.99 & as.numeric(fin_AC1[,6]) <= -0.1) | (as.numeric(fin_AC1[,8]) < 0.99  & as.numeric(fin_AC1[,6]) > -0.1 & as.numeric(fin_AC1[,6]) < 0.1) | (as.numeric(fin_AC1[,8]) > 1.05  & as.numeric(fin_AC1[,6]) > -0.1 & as.numeric(fin_AC1[,6]) < 0.1) | (as.numeric(fin_AC1[,8]) <= 1.05 & as.numeric(fin_AC1[,6]) >= 0.1)),]
#How many predicitons were correct:
print(nrow(CNV_T)/(nrow(CNV_T)+nrow(CNV_F))*100)
#True positives
CNV_1 <- fin_AC1[which((as.numeric(fin_AC1[,8]) <= 0.99 & as.numeric(fin_AC1[,6]) <= -0.1) | (as.numeric(fin_AC1[,8]) >= 1.05 & as.numeric(fin_AC1[,6]) >= 0.1)),]
#True negatives
CNV_2 <- fin_AC1[which((as.numeric(fin_AC1[,8]) > 0.99 & as.numeric(fin_AC1[,8]) < 1.05 & as.numeric(fin_AC1[,6]) > -0.1 & as.numeric(fin_AC1[,6]) < 0.1)),]
#False positives
CNV_3 <- fin_AC1[which((as.numeric(fin_AC1[,8]) > 0.99 & as.numeric(fin_AC1[,6]) <= -0.1) | (as.numeric(fin_AC1[,8]) <= 1.05 & as.numeric(fin_AC1[,6]) >= 0.1)),]
#False negatives
CNV_4 <- fin_AC1[which((as.numeric(fin_AC1[,8]) < 0.99  & as.numeric(fin_AC1[,6]) > -0.1 & as.numeric(fin_AC1[,6]) < 0.1) | (as.numeric(fin_AC1[,8]) > 1.05  & as.numeric(fin_AC1[,6]) > -0.1 & as.numeric(fin_AC1[,6]) < 0.1)),]
#Specificity:
print(nrow(CNV_2)/(nrow(CNV_2)+nrow(CNV_3))*100)
#Sensitivity:
print(nrow(CNV_1)/(nrow(CNV_1)+nrow(CNV_3))*100)

