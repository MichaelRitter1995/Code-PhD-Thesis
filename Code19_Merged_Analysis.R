# integrate all samples into one Seurat object

library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(stringr)
library(SCpubr)

setwd("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/Combined_Analysis/sn/analysis_merged/")

### read in data
#sample.names <- dir("analysis_per_sample")
#sample.names <- dir("analysis_per_sample", pattern= "GB[1-9]+(_[A, B])?.rds", recursive = T, full.names = T)
sample.names <- dir("analysis_per_sample", pattern= "DR_0.3*", recursive = T, full.names = T)
all.data <- lapply(sample.names, readRDS)

merged <- merge(all.data[[1]], all.data[-1], merge.data = F)
merged <- SCTransform(merged, vars.to.regress = c("nCount_RNA"), verbose = F)

print("Dimensional reduction ...")
# merged <- FindVariableFeatures(merged)
VariableFeatures(merged) <- rownames(merged@assays$SCT@scale.data)
merged <- RunPCA(merged, features = VariableFeatures(object = merged), verbose = F)
merged <- FindNeighbors(merged, dims = 1:30, verbose = F)
merged <- FindClusters(merged, resolution = .6, verbose = F)
merged <- RunUMAP(merged, dims = 1:30, verbose = F)
DimPlot(merged, reduction = "umap", label = T)

saveRDS(merged, "Combined_Analysis/merged_strict_DR.rds")

values=cbind(c(0,5,10), c(0.1,2,4), c(0.1,1,2))

for(i in 1:3){
  for(j in 1:3){
    for(k in 1:3){
      t <- values[i,1]
      l <- values[j,2]
      s <- values[k,3]
      harmon <- RunHarmony(merged, c("orig.ident"), assay.use = "SCT", theta=t, lamda=l, sigma=s)
      harmon <- RunUMAP(harmon, reduction = "harmony", dims = 1:20, verbose = F)
      harmon <- FindNeighbors(harmon, reduction = "harmony", dims = 1:20, verbose = F)
      harmon <- FindClusters(harmon, reduction = "harmony", resolution = .6, verbose = F)
      pdf(paste0("UMAP_merged_GB_sn_t",t,"_l",l,"_s",s,".pdf"), width = 6.5, height = 6)
      print(DimPlot(harmon, group.by = "orig.ident"))
      dev.off()
      saveRDS(harmon, paste0("harmony_strictt_DR_t",t,"_l",l,"_s",s,".rds"))
      }
  }
}

#Re-Analyze lambda

values1=cbind(c(0,5,10), c(0.1,10,30), c(0.1,1,2))

for(j in 1:3){
  l <- values1[j,2]
  harmon <- RunHarmony(merged, c("orig.ident"), assay.use = "SCT", theta=0, lambda=l, sigma=0.1)
  harmon <- RunUMAP(harmon, reduction = "harmony", dims = 1:20, verbose = F)
  harmon <- FindNeighbors(harmon, reduction = "harmony", dims = 1:20, verbose = F)
  harmon <- FindClusters(harmon, reduction = "harmony", resolution = .6, verbose = F)
  pdf(paste0("UMAP_merged_GB_sn_t0_l",l,"_s0.1.pdf"), width = 6.5, height = 6)
  print(DimPlot(harmon, group.by = "orig.ident"))
  dev.off()
  saveRDS(harmon, paste0("harmony_strictt_DR_t0_l",l,"_s0.1.rds"))
}

harmon <- RunHarmony(merged, c("orig.ident"), assay.use = "SCT", theta=0, lambda=NULL, sigma=0.1)
harmon <- RunUMAP(harmon, reduction = "harmony", dims = 1:20, verbose = F)
harmon <- FindNeighbors(harmon, reduction = "harmony", dims = 1:20, verbose = F)
harmon <- FindClusters(harmon, reduction = "harmony", resolution = .6, verbose = F)
pdf(paste0("UMAP_merged_GB_sn_t0_l_s0.1.pdf"), width = 6.5, height = 6)
print(DimPlot(harmon, group.by = "orig.ident"))
dev.off()
saveRDS(harmon, paste0("harmony_strictt_DR_t0_l_s0.1.rds"))

#harmon1 <- readRDS("Combined_Analysis/merged_harmony.rds")

harmon <- readRDS("harmony_strictt_DR_t0_l0.1_s0.1.rds")

#Col20A1 and PDGFRA for OPC; AQP4 for astrocytes
FeaturePlot(harmon, features=c("GFAP", "SYNPR", "CD163", "MBP", "PDGFRA", "COL20A1", "AQP4", "CLDN5"))
FeaturePlot(harmon, features=c("nCount_RNA"))

#####Annotate Cells#####
Cell_Types<- list(
  c("MBP","OLIG2","OLIG3","CLDN11","MOG","SOX10","GJB1","MAG","CNP","PLP1"),
  c("ENO2","NEFL","NEFM","NEFH","INA","SLC17A7","RBFOX3","SNAP25","SYT1","NEUROD2", "MAP2", "RBFOX1", "SNAP25"),
  c("NCAN", "NTSR2", "ADGRV1", "SLC1A2", "SLCO1C1", "ACSBG1", "GFAP", "MLC1", "PIRT", "SLC7A10", "AQP4"),
  c("TNR", "LHFPL3", "KLHL1", "PTPRZ1", "PLPP4", "FBN3", "MMP16", "TLL1", "COL11A1", "COL9A1"),
  c("FGF20", "P2RY12", "SMIM35", "ABCC4", "MUC22", "RASGEF1C"),
  c("RNASE2", "CEBPE", "CD163", "RETN", "CD207", "MSR1", "AIF1"),
  c("LCN6", "SELE", "CLDN5", "NOSTRIN", "ICAM2"),
  c("CDH6", "NR2F2", "PLEC1", "RFTN1"),
  c("CEMIP", "COL12A1", "FBLN1", "IGFBP5", "TRPC6"),
  c("EREG", "LIPN", "RIPK2", "MAPK7", "LILRA5", "IRAK3", "IL1B","C5AR1"),
  c("SERPINA9", "BLK", "FCRL1", "BANK1", "CD79B", "PAX5", "CNR2", "CXCR5"),
  c("CD3D", "CD8A", "CD8B"),
  c("CD160","KLRD1","NCR1","NKG7","KRT37","FGFBP2")
)

Cell_type_names <- c("Oligodendrocytes","Neurons", "Astrocytes", "OPC", "Microglia", "Macrophages", "Endothelia", "VLMC", "Monocytes", "B cells", "T cells", "NK")
names(Cell_Types) <- Cell_type_names
harmon <- AddModuleScore(object =harmon, features= Cell_Types, name = Cell_type_names, assay ='SCT')
f1 <- function(x){
  paste0(x, seq_along(x))
  }
l1<- f1(Cell_type_names)

FeaturePlot(harmon, features = unlist(l1))
DimPlot(harmon, reduction = "umap", label = T)

harmon$Annot_v1 <- "Tumour"
harmon$Annot_v1[which(harmon$seurat_clusters %in% c(8,5))] <- "Oligodendrocytes"
harmon$Annot_v1[which(harmon$seurat_clusters %in% c(0,18,13,26,6))] <- "Neurons"
harmon$Annot_v1[which(harmon$seurat_clusters %in% c(22,23))] <- "Lymphocytes"
harmon$Annot_v1[which(harmon$seurat_clusters %in% c(11,25))] <- "Endothelia"
harmon$Annot_v1[which(harmon$seurat_clusters %in% c(9))] <- "Astrocytes"
harmon$Annot_v1[which(harmon$seurat_clusters %in% c(16))] <- "OPC"
harmon$Annot_v1[which(harmon$seurat_clusters %in% c(2,14))] <- "Microglia"
#harmon$Annot_v1[which(harmon$seurat_clusters %in% c(22,23))] <- "T cells"
#harmon$Annot_v1[which(harmon$seurat_clusters %in% c(22))] <- "Monocytes"
harmon$Annot_v1[which(harmon$seurat_clusters %in% c(24))] <- "VLMC"
harmon$Annot_v1[which(harmon$seurat_clusters %in% c(19))] <- "Pericytes"



###Annotate Lymphocytes

Lymphocytes <- subset(x = harmon, subset= Annot_v1 == "Lymphocytes")
Lymphocytes <- RunUMAP(Lymphocytes, dims = 1:20, verbose = F)
Lymphocytes <- FindNeighbors(Lymphocytes, dims = 1:20, verbose = F)
Lymphocytes <- FindClusters(Lymphocytes, resolution = 2.0, verbose = F)
Lymphocyte_Types<- list(
c("EREG", "LIPN", "RIPK2", "MAPK7", "LILRA5", "IRAK3", "IL1B","C5AR1"),
c("SERPINA9", "BLK", "FCRL1", "BANK1", "CD79B", "PAX5", "CNR2", "CXCR5"),
c("CD3D", "CD8A", "CD8B"),
c("CD160","KLRD1","NCR1","NKG7","KRT37","FGFBP2"),
c("CCR7", "SELL")
)

Cell_type_names <- c("Monocytes", "B cells", "T cells", "NK")
names(Lymphocyte_Types) <- Cell_type_names
Lymphocytes <- AddModuleScore(object =Lymphocytes, features= Lymphocyte_Types, name = Cell_type_names, assay ='SCT')
f1 <- function(x){
  paste0(x, seq_along(x))
}
l1<- f1(Cell_type_names)

FeaturePlot(Lymphocytes, features = unlist(l1))
DimPlot(Lymphocytes, reduction = "umap", label = T)
DimPlot(Lymphocytes, group.by="orig.ident")
Lymphocytes <- PrepSCTFindMarkers(Lymphocytes)
Lymphocytes1 <- FindAllMarkers(object = Lymphocytes, recorrect_umi=FALSE) 
Markers10 <- Lymphocytes1[which(Lymphocytes1$p_val_adj <= 0.001 & Lymphocytes1$cluster == 10),]
Markers8 <- Lymphocytes1[which(Lymphocytes1$p_val_adj <= 0.001 & Lymphocytes1$cluster == 8),]

Lymphocytes$Annot_Lymph <- "Unknown"
Lymphocytes$Annot_Lymph[which(Lymphocytes$seurat_clusters %in% c(9))] <- "B cells"
Lymphocytes$Annot_Lymph[which(Lymphocytes$seurat_clusters %in% c(4,7))] <- "classical monocytes"
Lymphocytes$Annot_Lymph[which(Lymphocytes$seurat_clusters %in% c(6,3))] <- "non-classical monocytes"
Lymphocytes$Annot_Lymph[which(Lymphocytes$seurat_clusters %in% c(0,1))] <- "T helper"
Lymphocytes$Annot_Lymph[which(Lymphocytes$seurat_clusters %in% c(2))] <- "T regs"
Lymphocytes$Annot_Lymph[which(Lymphocytes$seurat_clusters %in% c(10))] <- "T cells"
DimPlot(Lymphocytes, group.by="Annot_Lymph")
DimPlot(Lymphocytes, group.by="orig.ident")

saveRDS(Lymphocytes, paste0("Lymphocyte_subset.rds"))


###Annotate Microglia

Microglia <- subset(x = harmon, subset= Annot_v1 == "Microglia")
Microglia  <- RunUMAP(Microglia , dims = 1:20, verbose = F)
Microglia  <- FindNeighbors(Microglia , dims = 1:20, verbose = F)
Microglia  <- FindClusters(Microglia , resolution = 1.0, verbose = F)
Microglia_Types<- list(
  c("ATP8B4", "MUC22", "RASGEF1C", "P2RY12", "ABCC4", "FGF20", "ABCC4", "KCNK13", "PRDM11", "SMIM35", "TBX18"),
  c("GPR183", "RNASET2", "NRP1", "TGFBI", "GPNMB"),
  c("RNASE2", "CEBPE", "MMP9", "CD163", "FCN1", "S100A12", "CCL18", "CCL17", "CD207", "RETN"),
  c("CD36", "RLF2", "LY86", "PECAM", "ITGAM", "CTSS", "CD48")
)

Cell_type_names <- c("Microglia","Macrophages", "Macrophages M1", "Macrophages M2")
names(Microglia_Types) <- Cell_type_names
Microglia <- AddModuleScore(object =Microglia, features= Microglia_Types, name = Cell_type_names, assay ='SCT')
f1 <- function(x){
  paste0(x, seq_along(x))
}
l1<- f1(Cell_type_names)
FeaturePlot(Microglia, features = unlist(l1))

DimPlot(Microglia)
DimPlot(Microglia, group.by="orig.ident")

Microglia$Annot_Micr <- "Microglia"
Microglia$Annot_Micr[which(Microglia$seurat_clusters %in% c(11,7,12,9,3))] <- "Macrophages"

DimPlot(Microglia, group.by="Annot_Micr")

###Annotate Neftel Subtypes

Cell_Types<- list(
  c("CHI3L1", "ANXA1", "ANXA2", "VIM", "HILPDA", "ADM", "ENO2", "MT1X", "ATF3"),
  c("CST3", "S100B", "SLC1A3", "HEPN1", "HOPX", "MT3", "SPARRCL1", "MCL1", "GFAP", "BCAN"),
  c("BCAN", "PLP1", "FIBIN", "OMG", "APOD", "SIRT2", "PLPPR1", "PTPRZ1", "VCAN", "TTYH1"),
  c("DLL3", "DLL1", "SOX4", "NEU4", "PTPRS", "DCX", "CD24", "CD200", "PAK3", "LBH", "NFIB", "TCF4")
)

Cell_type_names <- c("MES_like", "AC_like", "OPC_like","NPC_like")
names(Cell_Types) <- Cell_type_names
harmon <- AddModuleScore(object =harmon, features= Cell_Types, name = Cell_type_names, assay ='SCT')
FeaturePlot(harmon, features = c("MES_like1", "AC_like2", "OPC_like3","NPC_like4"))
DimPlot(harmon, label=T)

Lymph_types <- Lymphocytes$Annot_Lymph
Micr <- Microglia$Annot_Micr

#Excitatory Neuron Markers
FeaturePlot(harmon, features = c("NEUROD6", "GAP43", "ANO3", "CDH22", "NELL2", "FSTL4"))

harmon$Annot_v2 <- harmon$Annot_v1
harmon$Annot_v2[which(harmon$seurat_clusters %in% c(12))] <- "Proliferating Tumour"
harmon$Annot_v2[which(harmon$seurat_clusters %in% c(6))] <- "Inhibitory Neurons"
harmon$Annot_v2[which(harmon$seurat_clusters %in% c(0,18,13,26))] <- "Excitatory Neurons"
harmon$Annot_v2[which(harmon$seurat_clusters %in% c(4))] <- "NPC/OPC-like tumour"
harmon$Annot_v2[which(harmon$seurat_clusters %in% c(20,3,7,1,10,21,15))] <- "AC-like tumour"
harmon$Annot_v2[which(harmon$seurat_clusters %in% c(20))] <- "MES1-like tumour"
harmon$Annot_v2[which(harmon$seurat_clusters %in% c(17))] <- "MES2-like tumour"

Annot_v2 <- data.frame(harmon$Annot_v2, colnames(harmon))
colnames(Annot_v2) <- c("Annot", "barcode")
Lymph_v2 <- data.frame(Lymphocytes$Annot_Lymph, colnames(Lymphocytes))
colnames(Lymph_v2) <- c("Annot", "barcode")
#Annot_v2 <- dplyr::left_join(Annot_v2, Lymph_v2, by="barcode")
Annot_v2$Annot[match(Lymph_v2$barcode, Annot_v2$barcode)] <- Lymph_v2$Annot
Micro_v2 <- data.frame(Microglia$Annot_Micr, colnames(Microglia))
colnames(Micro_v2) <- c("Annot", "barcode")
Annot_v2$Annot[match(Micro_v2$barcode, Annot_v2$barcode)] <- Micro_v2$Annot

harmon$Annot_v2 <- Annot_v2$Annot
SCpubr::do_DimPlot(harmon, reduction = "umap", group.by= "Annot_v2")


saveRDS(harmon, paste0("harmon_strict_DR_Annot.rds"))

p1 <- SCpubr::do_FeaturePlot(harmon, features = c("FGFR3", "EGFR", "PCDH15"), ncol=3)
ggsave(plot = p1,filename = "Marker_Expression.pdf", width = 24, height = 8)

p2 <- SCpubr::do_DimPlot(harmon, group.by = "Annot_v1")
ggsave(plot = p2,filename = "Cell_Annotation.pdf", width = 8, height = 8)

p3 <- SCpubr::do_ViolinPlot(harmon, features = c("FGFR3", "EGFR","PCDH15"), group.by="Annot_v1")
ggsave(plot = p3,filename = "Marker_Expression_Annotation.pdf", width = 24, height = 8)

####CNV Analysis####

#Create Reference for CNV calling

ref <- subset(x = harmon, subset= Annot_v1 != "Tumour" & Annot_v1 != "OPC" & Annot_v1 !="Astrocytes")
keep_Olig <- sample(colnames(subset(ref, subset= Annot_v1=="Oligodendrocytes")), 2000)
keep_Neurons <- sample(colnames(subset(ref, subset= Annot_v1=="Neurons")), 2000)
keep_Microglia <- sample(colnames(subset(ref, subset= Annot_v1=="Microglia")), 2000)
keep_other <- colnames(subset(x = ref, subset= Annot_v1 != "Oligodendrocytes" & Annot_v1 != "Neurons" & Annot_v1 !="Microglia"))
keep <- c(keep_Olig, keep_Neurons, keep_Microglia, keep_other)
ref <- subset(x = ref, cells= keep)
ref$Barcodes <- colnames(ref)
names_ref <- ref$Annot_v1
t <- f1(names_ref)
colnames(ref) <- t

saveRDS(ref, paste0("Reference_CNV_inference.rds"))

#Run infercnv
#plot results

#load infercnv object
harmon <- readRDS("harmon_strict_DR_Annot.rds")
infer.cnv <- readRDS("CNVcalling/outs/GB_FFPE_MP_i6_wind201/run.final.infercnv_obj")

cnvs <- infer.cnv@expr.data[,unlist(infer.cnv@observation_grouped_cell_indices)]
cnvs <- cnvs[,match(colnames(cnvs),colnames(harmon))]


#Add all chromosomes
for(i in 1:length(unique(infer.cnv@gene_order$chr))){
  chr <- as.character(unique(infer.cnv@gene_order$chr)[[i]])
  temp <- colMeans(cnvs[which(infer.cnv@gene_order$chr == chr), ])
  harmon[[chr]] <- temp
}

#FeaturePlot(harmon, "chr7",
#            min.cutoff = "q5", max.cutoff = "q95")+ coord_fixed() +
#  scale_color_gradientn(colors = c( rev(RColorBrewer::brewer.pal(name = "RdBu",n=11))), limits = c(0.92,1.08)) + NoAxes()

p1 <- SCpubr::do_FeaturePlot(harmon, "chr7")+ coord_fixed() +
  scale_color_gradientn(colors = c( rev(RColorBrewer::brewer.pal(name = "RdBu",n=11))), limits = c(0.92,1.08)) + NoAxes()
ggsave(plot = p1,width = 10,height = 8, filename ="Harmon_Chr_7.pdf")

p2 <- SCpubr::do_FeaturePlot(harmon, "chr10")+ coord_fixed() +
  scale_color_gradientn(colors = c( rev(RColorBrewer::brewer.pal(name = "RdBu",n=11))), limits = c(0.92,1.08)) + NoAxes()
ggsave(plot = p2, width = 10,height = 8, filename = "Harmon_Chr_10.pdf")

p2 <- SCpubr::do_FeaturePlot(harmon, "chr8")+ coord_fixed() +
  scale_color_gradientn(colors = c( rev(RColorBrewer::brewer.pal(name = "RdBu",n=11))), limits = c(0.92,1.08)) + NoAxes()
ggsave(plot = p2, width = 10,height = 8, filename = "Harmon_Chr_10.pdf")

#plot changes according to sample
p <- FeaturePlot(harmon, "chr10", split.by = "orig.ident", combine = F)
plot <- list()
for(i in 1:length(p)){
  temp <- p[[i]]+ coord_fixed() + scale_color_gradientn(colors = c( rev(RColorBrewer::brewer.pal(name = "RdBu",n=11))), limits = c(0.92,1.08)) + NoAxes()
  plot[[length(plot)+1]] = temp 
}
p3 <- CombinePlots(plot, ncol = 4)
ggsave(plot = p3, width = 40,height = 32, filename = "Harmon_Chr_10_by_sample.pdf")

p <- FeaturePlot(harmon, "chr7", split.by = "orig.ident", combine = F)
plot <- list()
for(i in 1:length(p)){
  temp <- p[[i]]+ coord_fixed() + scale_color_gradientn(colors = c( rev(RColorBrewer::brewer.pal(name = "RdBu",n=11))), limits = c(0.92,1.08)) + NoAxes()
  plot[[length(plot)+1]] = temp 
}
p4 <- CombinePlots(plot, ncol = 4)
ggsave(plot = p4, width = 40,height = 32, filename = "Harmon_Chr_7_by_sample.pdf")

p <- FeaturePlot(harmon, "FGFR3", split.by = "orig.ident", combine = F)
p5 <- CombinePlots(p, ncol = 4)
ggsave(plot = p5, width = 40,height = 32, filename = "Harmon_FGFR3_by_sample.pdf")

p <- FeaturePlot(harmon, "PCDH15", split.by = "orig.ident", combine = F)
p6 <- CombinePlots(p, ncol = 4)
ggsave(plot = p6, width = 40,height = 32, filename = "Harmon_PCDH15_by_sample.pdf")


p <- DimPlot(harmon, group.by="Annot_v1", split.by = "orig.ident", combine = F, ncol=4)
p7 <- CombinePlots(p, ncol = 4)

#####Analyse Neurons#####

Neurons <- subset(x = harmon, subset = Annot_v1 == "Neurons")
Neurons <- RunHarmony(Neurons, c("orig.ident"), assay.use = "SCT")
Neurons <- RunUMAP(Neurons, reduction = "harmony", dims = 1:20, verbose = F)
Neurons <- FindNeighbors(Neurons, reduction = "harmony", dims = 1:20, verbose = F)
Neurons <- FindClusters(Neurons, reduction = "harmony", resolution = .6, verbose = F)
DimPlot(Neurons, reduction = "umap", label = T)

FeaturePlot(Neurons, features = c("VIP", "LAMP5", "NTNG2", "SST", "KIT", "HTR2C", "KCNH5"))


