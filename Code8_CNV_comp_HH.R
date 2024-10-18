library(tidyverse)
library(SPATA2)
library(igraph)
library(reticulate)
library(ggpubr)
library(dplyr)
library(infercnv)
library(RColorBrewer)


setwd("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")

#set colorPalette to match infercnv colours
col=colorRampPalette(rev(c(RColorBrewer::brewer.pal(9, "RdBu"))))

#path to 850K data
path850K <- ("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/Methylation_data/")
#path to Visium file
pathVis <- ("/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/")


#Cryo samples are missing
ID_list=data.frame(c("GB_STX2", "GB_STX7", "GB_STX1", "GB_STX6", "GB_STX4",
                     "195096", "GB2", "GB6", "GB11", "GB13",
                     "N20_2382", "GB1", "GB10", "GB12", "GB_STX3", "N21_2972",
                      "GB15", "GB17_A", "GB18", "GB19", 
                     "GB20","GB22", "GB5"))

#Missing inferred CNVs in SPATA object:
#"GB3", "GB7"

#Not included for being astrocytoma
#"AC_NOS_STX1", "AC_STX6"

#Not included because of low UMI count -> noisy CNV calling:
#"GB7", "AC", "OG", "AC_OG", "MG", "GB14", 

#Not included: Mainly Healthy or no CNVs detected
# "GB9", "N21_1260_IB", "GB17_B", "GB16", "GB21"




#load list of alias
list_alias <- read.csv(file="/omics/groups/OE0146/internal/Micha/Analysis/GBM/Per_Sample/alias.txt")
colnames(list_alias) <- c("ID", "Alias")

list_of_dataframes = lapply(ID_list[,1] , function (x) {
  
  #adapt variables for samples with 2 patients
  ifelse(x=="GB_STX3" | x=="N21_2972", sample_name <- "GB_STX3_2972", sample_name <- x)
  print(x)
  #load spata object
  spata_obj <- updateSpataObject(readRDS(paste0(pathVis, x, "/SPATA_obj_", x, ".RDS")))
  #load 850K data (.igv file)
  Meth_bins <- read.csv(paste0(path850K,x,"/",x,".igv"), sep="\t")
  colnames(Meth_bins) <-c("Chr","start","end","bin","value")
  
  #####Plot inferred CNVs like 850K
  
  #load barcodes with high tumour content used for CNV comparison
  #roughly 50-400 spots per sample -> less pooling of subclones, same number of spots per samples (stereotactic vs. bulk)
  bc_CNV <- read.csv(paste0(x,"/CNV.csv"))
  bc_CNV[,3] <- rep("Spots selected for CNV analysis", nrow(bc_CNV))
  Missing_barcodes <- spata_obj@fdata[[paste0("s_",sample_name)]][["barcodes"]]
  Missing_barcodes <- Missing_barcodes[which(!(spata_obj@fdata[[paste0("s_",sample_name)]][["barcodes"]] %in% bc_CNV$Barcode))]
  Missing_barcodes <- data.frame(cbind(Missing_barcodes, rep(NA, length(Missing_barcodes)), rep(NA, length(Missing_barcodes))))
  colnames(Missing_barcodes) <- colnames(bc_CNV)
  bc_CNV1 <- rbind(bc_CNV, Missing_barcodes)
  spata_obj@fdata[[paste0("s_",sample_name)]][["CNV_calling"]] <- bc_CNV1[match(spata_obj@fdata[[paste0("s_",sample_name)]]$barcodes, bc_CNV1$Barcode),3]
  
  p3 <-
    plotSurface(spata_obj, color_by="CNV_calling", pt_alpha=0.5)+
    scale_color_manual(values=RColorBrewer::brewer.pal(10, "Set3")[4], na.value="grey", name="")+
    ggpLayerThemeCoords() +
    ggpLayerAxesSI(spata_obj, unit = "mm", breaks = str_c(1:7, "mm"), add_labs = TRUE)+
    ggtitle("Spots used for CNV Inference")                  
  p3
  ggsave(filename = paste0("XX_", x, "_inferCNV_spots.pdf"), plot = p3, path= paste0(x, "/"), height = 5, width = 7)
  saveRDS(p3$data, file = paste0(x, "/", x, "_Spots_CNV_inference.RDS"))
  
  #Subset CNVs according to barcodes containing the tumour tissue
  # => there is probably an easier way
  spata_sub <- subsetByBarcodes(spata_obj, barcodes = bc_CNV$Barcode)
  spata_sub@cnv[[paste0("s_",sample_name)]][["cnv_df"]] <- spata_sub@cnv[[paste0("s_",sample_name)]][["cnv_df"]][which(spata_sub@cnv[[paste0("s_",sample_name)]][["cnv_df"]][["barcodes"]] %in% bc_CNV$Barcode),]
  spata_sub@cnv[[paste0("s_",sample_name)]][["Normalized.bin.matrix"]] <- spata_sub@cnv[[paste0("s_",sample_name)]][["Normalized.bin.matrix"]][,which(colnames(spata_sub@cnv[[paste0("s_",sample_name)]][["Normalized.bin.matrix"]]) %in% bc_CNV$Barcode)]
  spata_sub@cnv[[paste0("s_",sample_name)]][["cnv_mtr"]] <- spata_sub@cnv[[paste0("s_",sample_name)]][["cnv_mtr"]][,which(colnames(spata_sub@cnv[[paste0("s_",sample_name)]][["cnv_mtr"]]) %in% bc_CNV$Barcode)]
  
  #plot CNVs per sample
  p_dat <- SPATAwrappers::plotCNV.RefMode(spata_sub)
  
  p <- ggplot(p_dat$data)+
    geom_point(mapping=aes(x=xaxsis, y=cnv_mean, color=cnv_mean), size=1)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5),
          axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))+
    scale_colour_gradientn(colours = col(50),limit=c(1, 3), oob = scales::squish, na.value ="white")+
    ylab("CNV")+xlab("Bin Index")+
    guides(color = guide_colourbar(barwidth = 0.3, barheight = 8, ticks =F, frame.colour="black", label=T))+
    geom_vline(xintercept=0, alpha=0.2)+
    ylim(-4, 12)+
    ggtitle("Inferred CNV")
  
  #add horizontal lines
  for(z in unique(p_dat$data$Chr)){
    p=p+geom_vline(xintercept=p_dat$data %>% filter(Chr==z) %>% tail(1) %>% pull(xaxsis), alpha=0.2)
  }
  
  saveRDS(p$data, file = paste0(x, "/", x, "_spatial_CNVs.RDS"))
  
  #extract CNV information of spatial data
  st_data <- p_dat$data
  
  
  #Prepare to match bins
  st_data$Chr <- paste0("chr", st_data$Chr)
  #only segemtns with more than 5 genes, probably not needed
  #st_data <- st_data[which(st_data$n>=5),]
  
  #match bins
  CNV_mean <- map(.x=1:length(st_data$bin), .f=function(i){
    chr <- st_data$Chr[i]
    s <- st_data$start[i]
    e <- st_data$end[i]
    mean_CNV <- Meth_bins %>% filter(Chr==chr) %>% filter(start>s & end<e) %>% pull(value) %>% mean()
    return(mean_CNV)
  }, .progress=T) %>% unlist()
  st_data$cnv_mean_meth <- CNV_mean
  
  #plot EPIC CNVs
  p2 <- ggplot(st_data)+
    geom_point(mapping=aes(x=xaxsis,y=cnv_mean_meth, color=cnv_mean_meth), size=1)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5),
          axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))+
    scale_colour_gradientn(colours = col(50),limit=c(-0.5, 0.5), oob = scales::squish, na.value ="white")+
    ylab("CNV")+xlab("Bin Index")+
    guides(color = guide_colourbar(barwidth = 0.3, barheight = 8, ticks =F, frame.colour="black", label=T))+
    geom_vline(xintercept=0, alpha=0.2)+
    ggtitle("EPIC CNV")
  
  #add vertical lines:
  for(z in unique(p_dat$data$Chr)){
    p2=p2+geom_vline(xintercept=p_dat$data %>% filter(Chr==z) %>% tail(1) %>% pull(xaxsis), alpha=0.2)
  }
  ggsave(filename = paste0("1_", x, "_850K_vs_Vis_CNVs.pdf"), plot = p2+p, path= paste0(x, "/"), height = 4, width = 15)
  saveRDS(p2$data, file = paste0(x, "/", x, "_EPIC_CNVs.RDS"))
  
  #plot comparison of 850K and Visium per bin:
  p3 <- ggplot(st_data)+
    geom_point(mapping=aes(x=cnv_mean_meth, y=cnv_mean), size=1, stroke=0)+
    coord_fixed(ratio=0.2)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5),
          axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))+
    ylab("CNV Scores Visium")+xlab("CNV Scores Methylome")
  
  #Add lines to separate gains and losses
  p3 <- p3 + geom_hline(yintercept=c(1.5,2.5), linetype="dotted")+geom_vline(xintercept=c(-0.1, 0.1), linetype="dotted")
  ggsave(filename = paste0("2_", x, "_comp_850K_vs_Vis.pdf"), plot = p3, path= paste0(x, "/"), height = 6, width = 6)
  
  #Add all CNV values into one dataframe
  st_data <- cbind(st_data, rep(x,nrow(st_data)))
}
)

#create one big dataframe containing all information of all samples
fin_raw <- dplyr::bind_rows(list_of_dataframes)
fin_raw <- na.omit(fin_raw)

#remove segments with less than 5 genes, needed?
fin <- fin_raw[which(fin_raw$n >= 5),]
#Annotate number of genes per segment
#fin_raw[, ncol(fin_raw)+1]=ifelse(fin_raw$n <5, "<5",ifelse(fin_raw$n <=10,"5-10", ifelse(fin_raw$n <=25,"10-25", ">25")))
fin_raw[, ncol(fin_raw)+1]=ifelse(fin_raw$n ==0, "0",ifelse(fin_raw$n <=2,"1-2", ifelse(fin_raw$n <=5,"3-5", ">5")))

colnames(fin_raw) <- c("bin", "start", "end", "Chr", "Chr.arm", "n", "xaxsis", "Arm", "cnv_mean", "cnv_sd", "cnv_mean_meth", "ID", "colour")
colnames(fin) <- c("bin", "start", "end", "Chr", "Chr.arm", "n", "xaxsis", "Arm", "cnv_mean", "cnv_sd", "cnv_mean_meth", "ID")
fin <- merge(fin, list_alias, by= "ID")
fin_raw <- merge(fin_raw, list_alias, by= "ID")

#Claculate standard deviation
No_CNV <- fin[which((as.numeric(fin$cnv_mean_meth) >= -0.075) & as.numeric(fin$cnv_mean_meth <= 0.075)),]
lim <- sd(fin$cnv_mean)/2

No_CNV <- fin[which((as.numeric(fin_raw$cnv_mean_meth) >= -0.075) & as.numeric(fin_raw$cnv_mean_meth <= 0.075)),]
lim_raw <- sd(fin_raw$cnv_mean)/2


#Colour according to number of genes per segment
p6 <- ggplot(fin_raw)+
  geom_point(mapping=aes(x=cnv_mean_meth, y=cnv_mean, colour=colour), size=1, alpha=0.5, stroke=0)+
  coord_fixed(ratio=0.2)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))+
  ylab("CNV Scores Visium")+xlab("CNV Scores Methylome")+scale_color_manual(values = c("red", "blue", "green3", "orange"))
p6 <- p6 + geom_hline(yintercept=c(1.5,2.5), linetype="dotted")+geom_vline(xintercept=c(-0.1, 0.1), linetype="dotted")
p6
ggsave(filename = "2_comp_850K_vs_Vis_genes_per_segment.pdf", plot = p6, path=getwd(), height = 6, width = 6)


#Colour according to sample

p7 <- ggplot(fin)+
  geom_point(mapping=aes(x=ifelse(cnv_mean_meth>-1, cnv_mean_meth, -1.1), y=cnv_mean, colour=Alias), size=1, alpha=0.5, stroke=0)+
  coord_fixed(ratio=0.2)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))+
  xlim(c(-1.1,1.1))+
  ylab("CNV Scores Visium")+xlab("CNV Scores Methylome")
p7 <- p7 + geom_hline(yintercept=c(2-lim,2+lim), linetype="dotted")+geom_vline(xintercept=c(-0.075, 0.075), linetype="dotted")
p7
ggsave(filename = "1_comp_850K_vs_Vis_all_samples.pdf", plot = p7, path=getwd(), height = 6, width = 6)



#Extract true positives (CNV)
CNV_1 <- fin[which((fin$cnv_mean_meth <= -0.075 & fin$cnv_mean <= 2-lim) | (fin$cnv_mean_meth >= 0.075 & fin$cnv_mean >= 2+lim)),]
#Extract true negatives (no CNV)
CNV_2 <- fin[which((fin$cnv_mean_meth > -0.075 & fin$cnv_mean_meth < 0.075 & fin$cnv_mean > 2-lim & fin$cnv_mean < 2+lim)),]
#Extract false negatives
CNV_3 <- fin[which((fin$cnv_mean_meth > -0.075 & fin$cnv_mean <= 2-lim) | (fin$cnv_mean_meth <= 0.075 & fin$cnv_mean >= 2+lim)),]
#Extract false positives
CNV_4 <- fin[which((fin$cnv_mean_meth < -0.075  & fin$cnv_mean > 2-lim & fin$cnv_mean < 2+lim) | (fin$cnv_mean_meth > 0.075  & fin$cnv_mean > 2-lim & fin$cnv_mean < 2+lim)),]

#Sensitivity:
nrow(CNV_1)/(nrow(CNV_1)+nrow(CNV_3))
#Specificity:
nrow(CNV_2)/(nrow(CNV_2)+nrow(CNV_4))

#All correct predictions
CNV_T <- fin[which((fin$cnv_mean <= 2-lim & fin$cnv_mean_meth <= -0.075) | (fin$cnv_mean > 2-lim & fin$cnv_mean < 2+lim & fin$cnv_mean_meth > -0.075 & fin$cnv_mean_meth < 0.075) | (fin$cnv_mean >= 2+lim & fin$cnv_mean_meth >= 0.075)),]
#All wrong predictions
CNV_F <- fin[which((fin$cnv_mean > 2-lim & fin$cnv_mean_meth <= -0.075) | (fin$cnv_mean < 2-lim  & fin$cnv_mean_meth > -0.075 & fin$cnv_mean_meth < 0.075) | (fin$cnv_mean > 2+lim  & fin$cnv_mean_meth > -0.075 & fin$cnv_mean_meth < 0.075) | (fin$cnv_mean <= 2+lim & fin$cnv_mean_meth >= 0.075)),]
#How many predicitons were correct:
print(nrow(CNV_T)/(nrow(CNV_T)+nrow(CNV_F))*100)

##For the raw data

#Extract true positives (CNV)
CNV_1 <- fin_raw[which((fin_raw$cnv_mean_meth <= -0.075 & fin_raw$cnv_mean <= 2-lim) | (fin_raw$cnv_mean_meth >= 0.075 & fin_raw$cnv_mean >= 2+lim)),]
#Extract true negatives (no CNV)
CNV_2 <- fin_raw[which((fin_raw$cnv_mean_meth > -0.075 & fin_raw$cnv_mean_meth < 0.075 & fin_raw$cnv_mean > 2-lim & fin_raw$cnv_mean < 2+lim)),]
#Extract false negatives
CNV_3 <- fin_raw[which((fin_raw$cnv_mean_meth > -0.075 & fin_raw$cnv_mean <= 2-lim) | (fin_raw$cnv_mean_meth <= 0.075 & fin_raw$cnv_mean >= 2+lim)),]
#Extract false positives
CNV_4 <- fin_raw[which((fin_raw$cnv_mean_meth < -0.075  & fin_raw$cnv_mean > 2-lim & fin_raw$cnv_mean < 2+lim) | (fin_raw$cnv_mean_meth > 0.075  & fin_raw$cnv_mean > 2-lim & fin_raw$cnv_mean < 2+lim)),]

#Sensitivity:
nrow(CNV_1)/(nrow(CNV_1)+nrow(CNV_3))
#Specificity:
nrow(CNV_2)/(nrow(CNV_2)+nrow(CNV_4))

#All correct predictions
CNV_T <- fin_raw[which((fin_raw$cnv_mean <= 2-lim & fin_raw$cnv_mean_meth <= -0.075) | (fin_raw$cnv_mean > 2-lim & fin_raw$cnv_mean < 2+lim & fin_raw$cnv_mean_meth > -0.075 & fin_raw$cnv_mean_meth < 0.075) | (fin_raw$cnv_mean >= 2+lim & fin_raw$cnv_mean_meth >= 0.075)),]
#All wrong predictions
CNV_F <- fin_raw[which((fin_raw$cnv_mean > 2-lim & fin_raw$cnv_mean_meth <= -0.075) | (fin_raw$cnv_mean < 2-lim  & fin_raw$cnv_mean_meth > -0.075 & fin_raw$cnv_mean_meth < 0.075) | (fin_raw$cnv_mean > 2+lim  & fin_raw$cnv_mean_meth > -0.075 & fin_raw$cnv_mean_meth < 0.075) | (fin_raw$cnv_mean <= 2+lim & fin_raw$cnv_mean_meth >= 0.075)),]
#How many predicitons were correct:
print(nrow(CNV_T)/(nrow(CNV_T)+nrow(CNV_F))*100)



