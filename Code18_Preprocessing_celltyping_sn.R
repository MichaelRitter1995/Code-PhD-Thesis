#####Clean-up per sample#####

library(scater)
library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(stringr)
library(DoubletFinder)
source("/omics/groups/OE0146/internal/Micha/Unlimited_Power/R/functions/Doublet_Finder.R")

setwd("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/Combined_Analysis/sn")

####GB1####

#sample.names <- dir("analysis_per_sample")
#sample.names <- dir("analysis_per_sample", pattern= "GB[1-9]+(_[A, B])?.rds", recursive = T, full.names = T)
sample <- "GB1"

sampledir <- "/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/aligned/GB/P2_3/outs/per_sample_outs/"
data <- Read10X(data.dir=paste0(sampledir,sample,"/count/sample_filtered_feature_bc_matrix"))
data <- CreateSeuratObject(counts = data, project = sample, min.cells = 0, min.features = 500)

#Get mitochondrial Reads
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2

#Filter for UMI and mt-counts

mito.keep <- data$percent.mt <= 3
feat.keep <- !(isOutlier(data$nCount_RNA, nmads=3, type="higher"))

p1 <- ggplot() +
  geom_histogram(aes(data$nFeature_RNA), colour="black", fill="white", binwidth = 100) +
  geom_vline(aes(xintercept=500), color="blue", linetype="dashed") +
  labs(x="Detected features per cell")
p1
p2 <- ggplot() +
  geom_histogram(aes(data$nCount_RNA), colour="black", fill="white", binwidth = 1000) +
  geom_vline(aes(xintercept=200), color="blue", linetype="dashed") +
  geom_vline(aes(xintercept=max(data$nFeature_RNA[feat.keep])), color="blue", linetype="dashed") +
  labs(x="Total counts per cell")
p2
p3 <- ggplot() +
  geom_histogram(aes(data$percent.mt), colour="black", fill="white", binwidth = .5) +
  geom_vline(aes(xintercept=max(data$percent.mt[mito.keep])), color="blue", linetype="dashed") +
  labs(x="Percentage of mtRNA") +
  xlim(c(0,30))
p3

keep <- mito.keep & feat.keep
data <- data[, keep]

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2

### normalization and scaling with SCTransform ###
print("Performing normalization, regression, scaling and dimensionality reduction...")
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F) # default options
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = F) # default options
data <- ScaleData(data, vars.to.regress = "nCount_RNA", verbose = F)
### dimensionality reduction
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = F)
### clustering and generating UMAP
data <- FindNeighbors(data, dims = 1:20, verbose = F)
data <- FindClusters(data, resolution = .6, verbose = F)
data <- RunUMAP(data, dims = 1:30, verbose = F)


###Run Doublet Finder###

print("Doublet detection...")
## pK Identification (no ground-truth)
sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk <- as.numeric(levels(bcmvn$pK))[which.max(bcmvn$BCmetric)]

## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
exp_dbl <- (ncol(data)*.8/1000)/100 ## according to 10X's expected multiplet rate depnding on the no. of cells
nExp_poi <- round(exp_dbl*nrow(data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
data.dbl <- doubletFinder_v3(data, PCs = 1:15, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

## create output files
dbl.scores <- colnames(data.dbl@meta.data)[grepl("pANN",colnames(data.dbl@meta.data), fixed = T)]
data.dbl@meta.data[["Doublets"]] <- ifelse( data.dbl@meta.data[[dbl.scores]] <= 0.35, "Singlet", "Doublet")
#dbl.binary <- colnames(data.dbl@meta.data)[grepl("DF.classifications",colnames(data.dbl@meta.data), fixed = T)]
DimPlot(data.dbl, group.by = dbl.binary) + ggtitle(paste0(table(data.dbl[[dbl.binary]])[[1]], " doublets detected"))

FeaturePlot(data.dbl, dbl.scores) + ggtitle("Doublet score") + scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits = c(0,0.35))
FeaturePlot(data.dbl, features = "nCount_RNA")+ scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits=c(2000,7000))
FeaturePlot(data.dbl, features = "nFeature_RNA")+ scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits=c(1000,5000))
FeaturePlot(data.dbl, features = "percent.mt")+ scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits=c(0,2))


data <- data[,which(data.dbl[["Doublets"]]=="Singlet")]
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = F)
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:15, verbose = F)
data <- FindClusters(data, resolution = .6, verbose = F)
data <- RunUMAP(data, dims = 1:15, verbose = F)
DimPlot(data, reduction = "umap", label = T)

#Using Feature vs RNA count matrix does not work...
#test <- data[,which(data[["seurat_clusters"]]=="0" | data[["seurat_clusters"]]=="1"| data[["seurat_clusters"]]=="6"|data[["seurat_clusters"]]=="11"| data[["seurat_clusters"]]=="3")]
#test <- data[,which(data[["seurat_clusters"]]=="0" | data[["seurat_clusters"]]=="1"| data[["seurat_clusters"]]=="6"|data[["seurat_clusters"]]=="11"| data[["seurat_clusters"]]=="3")]
#plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", split.by = "seurat_clusters", ncol = 4)
#plot2 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", ncol = 3)
#plot2

Cell_Types<- list(
  c("MBP","OLIG2","OLIG3","CLDN11","MOG","SOX10","GJB1","MAG","CNP","PLP1"),
  c("ENO2","NEFL","NEFM","NEFH","INA","SLC17A7","RBFOX3","SNAP25","SYT1","NEUROD2", "MAP2", "RBFOX1", "SNAP25"),
  c("NCAN", "NTSR2", "ADGRV1", "SLC1A2", "SLCO1C1", "ACSBG1", "GFAP", "MLC1", "PIRT", "SLC7A10", "AQP4"),
  c("TNR", "LHFPL3", "KLHL1", "PTPRZ1", "PLPP4", "FBN3", "MMP16", "TLL1", "COL11A1", "COL9A1"),
  c("FGF20", "P2RY12", "SMIM35", "ABCC4", "MUC22", "RASGEF1C"),
  c("RNASE2", "CEBPE", "CD163", "RETN", "CD207", "MSR1", "AIF1"),
  c("LCN6", "SELE", "CLDN5", "NOSTRIN", "ICAM2")
)

Cell_type_names <- c("Oligodendrocytes","Neurons", "Astrocytes", "OPC", "Microglia", "Macrophages", "Endothelia")
names(Cell_Types) <- Cell_type_names
data <- AddModuleScore(object = data, features= Cell_Types, name = Cell_type_names, assay ='RNA')
f1 <- function(x){
  paste0(x, seq_along(x))
}
l1<- f1(Cell_type_names)
FeaturePlot(data, features = unlist(l1))
DimPlot(data, reduction = "umap", label = T)

GB1 <- data



#####GB9_A#####

#sample.names <- dir("analysis_per_sample")
#sample.names <- dir("analysis_per_sample", pattern= "GB[1-9]+(_[A, B])?.rds", recursive = T, full.names = T)
sample <- "GB9_A"

sampledir <- "/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/aligned/GB/P2_3/outs/per_sample_outs/"
data <- Read10X(data.dir=paste0(sampledir,sample,"/count/sample_filtered_feature_bc_matrix"))
data <- CreateSeuratObject(counts = data, project = sample, min.cells = 0, min.features = 500)

#Get mitochondrial Reads
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2

#Filter for UMI and mt-counts

mito.keep <- data$percent.mt <= 4
feat.keep <- !(isOutlier(data$nCount_RNA, nmads=3, type="higher"))

p1 <- ggplot() +
  geom_histogram(aes(data$nFeature_RNA), colour="black", fill="white", binwidth = 100) +
  geom_vline(aes(xintercept=500), color="blue", linetype="dashed") +
  labs(x="Detected features per cell")
p1
p2 <- ggplot() +
  geom_histogram(aes(data$nCount_RNA), colour="black", fill="white", binwidth = 1000) +
  geom_vline(aes(xintercept=200), color="blue", linetype="dashed") +
  geom_vline(aes(xintercept=max(data$nFeature_RNA[feat.keep])), color="blue", linetype="dashed") +
  labs(x="Total counts per cell")
p2
p3 <- ggplot() +
  geom_histogram(aes(data$percent.mt), colour="black", fill="white", binwidth = .5) +
  geom_vline(aes(xintercept=max(data$percent.mt[mito.keep])), color="blue", linetype="dashed") +
  labs(x="Percentage of mtRNA") +
  xlim(c(0,30))
p3

keep <- mito.keep & feat.keep
data <- data[, keep]

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2




#####GB4_A#####
#sample.names <- dir("analysis_per_sample")
#sample.names <- dir("analysis_per_sample", pattern= "GB[1-9]+(_[A, B])?.rds", recursive = T, full.names = T)
sample <- "GB4_A"

sampledir <- "/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/aligned/GB/P2_3/outs/per_sample_outs/"
data <- Read10X(data.dir=paste0(sampledir,sample,"/count/sample_filtered_feature_bc_matrix"))
data <- CreateSeuratObject(counts = data, project = sample, min.cells = 0, min.features = 500)

#Get mitochondrial Reads
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2

#Filter for UMI and mt-counts

mito.keep <- data$percent.mt <= 3
feat.keep <- !(isOutlier(data$nCount_RNA, nmads=3, type="higher"))

p1 <- ggplot() +
  geom_histogram(aes(data$nFeature_RNA), colour="black", fill="white", binwidth = 100) +
  geom_vline(aes(xintercept=500), color="blue", linetype="dashed") +
  labs(x="Detected features per cell")
p1
p2 <- ggplot() +
  geom_histogram(aes(data$nCount_RNA), colour="black", fill="white", binwidth = 1000) +
  geom_vline(aes(xintercept=200), color="blue", linetype="dashed") +
  geom_vline(aes(xintercept=max(data$nFeature_RNA[feat.keep])), color="blue", linetype="dashed") +
  labs(x="Total counts per cell")
p2
p3 <- ggplot() +
  geom_histogram(aes(data$percent.mt), colour="black", fill="white", binwidth = .5) +
  geom_vline(aes(xintercept=max(data$percent.mt[mito.keep])), color="blue", linetype="dashed") +
  labs(x="Percentage of mtRNA") +
  xlim(c(0,30))
p3

keep <- mito.keep & feat.keep
data <- data[, keep]

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2



### normalization and scaling with SCTransform
print("Performing normalization, regression, scaling and dimensionality reduction...")
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F) # default options
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = F) # default options
data <- ScaleData(data, vars.to.regress = "nCount_RNA", verbose = F)
### dimensionality reduction
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = F)
### clustering and generating UMAP
data <- FindNeighbors(data, dims = 1:20, verbose = F)
data <- FindClusters(data, resolution = .6, verbose = F)
data <- RunUMAP(data, dims = 1:30, verbose = F)
DimPlot(data)

###Run Doublet Finder###

print("Doublet detection...")
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk <- as.numeric(levels(bcmvn$pK))[which.max(bcmvn$BCmetric)]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
exp_dbl <- (ncol(data)*.8/1000)/100 ## according to 10X's expected multiplet rate depnding on the no. of cells
nExp_poi <- round(exp_dbl*nrow(data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
data.dbl <- doubletFinder_v3(data, PCs = 1:15, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

## create output files -------------------------------------------------------------------------------------------------------
dbl.scores <- colnames(data.dbl@meta.data)[grepl("pANN",colnames(data.dbl@meta.data), fixed = T)]
data.dbl@meta.data[["Doublets"]] <- ifelse( data.dbl@meta.data[[dbl.scores]] <= 0.35, "Singlet", "Doublet")
#dbl.binary <- colnames(data.dbl@meta.data)[grepl("DF.classifications",colnames(data.dbl@meta.data), fixed = T)]
#DimPlot(data.dbl, group.by = dbl.binary) + ggtitle(paste0(table(data.dbl[[dbl.binary]])[[1]], " doublets detected"))

FeaturePlot(data.dbl, dbl.scores) + ggtitle("Doublet score") + scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits = c(0,0.35))
FeaturePlot(data.dbl, features = "nCount_RNA")+ scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits=c(2000,7000))
FeaturePlot(data.dbl, features = "nFeature_RNA")+ scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits=c(1000,5000))
FeaturePlot(data.dbl, features = "percent.mt")+ scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits=c(0,2))


data <- data[,which(data.dbl[["Doublets"]]=="Singlet")]
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = F)
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:15, verbose = F)
data <- FindClusters(data, resolution = .6, verbose = F)
data <- RunUMAP(data, dims = 1:15, verbose = F)
DimPlot(data, reduction = "umap", label = T)

#Using Feature vs RNA count matrix does not work...
#test <- data[,which(data[["seurat_clusters"]]=="0" | data[["seurat_clusters"]]=="1"| data[["seurat_clusters"]]=="6"|data[["seurat_clusters"]]=="11"| data[["seurat_clusters"]]=="3")]
#test <- data[,which(data[["seurat_clusters"]]=="0" | data[["seurat_clusters"]]=="1"| data[["seurat_clusters"]]=="6"|data[["seurat_clusters"]]=="11"| data[["seurat_clusters"]]=="3")]
#plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", split.by = "seurat_clusters", ncol = 4)
#plot2 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", ncol = 3)
#plot2

Cell_Types<- list(
  c("MBP","OLIG2","OLIG3","CLDN11","MOG","SOX10","GJB1","MAG","CNP","PLP1"),
  c("ENO2","NEFL","NEFM","NEFH","INA","SLC17A7","RBFOX3","SNAP25","SYT1","NEUROD2", "MAP2", "RBFOX1", "SNAP25"),
  c("NCAN", "NTSR2", "ADGRV1", "SLC1A2", "SLCO1C1", "ACSBG1", "GFAP", "MLC1", "PIRT", "SLC7A10", "AQP4"),
  c("TNR", "LHFPL3", "KLHL1", "PTPRZ1", "PLPP4", "FBN3", "MMP16", "TLL1", "COL11A1", "COL9A1"),
  c("FGF20", "P2RY12", "SMIM35", "ABCC4", "MUC22", "RASGEF1C"),
  c("RNASE2", "CEBPE", "CD163", "RETN", "CD207", "MSR1", "AIF1"),
  c("LCN6", "SELE", "CLDN5", "NOSTRIN", "ICAM2")
)

Cell_type_names <- c("Oligodendrocytes","Neurons", "Astrocytes", "OPC", "Microglia", "Macrophages", "Endothelia")
names(Cell_Types) <- Cell_type_names
data <- AddModuleScore(object = data, features= Cell_Types, name = Cell_type_names, assay ='RNA')
f1 <- function(x){
  paste0(x, seq_along(x))
}
l1<- f1(Cell_type_names)
FeaturePlot(data, features = unlist(l1))
FeaturePlot(data, features = c("LAMP5", "HTR2C", "NTNG2", "VIP", "KIT", "SST"))
DimPlot(data, reduction = "umap", label = T)

GB4_A <- data




#####GB8####

#sample.names <- dir("analysis_per_sample")
#sample.names <- dir("analysis_per_sample", pattern= "GB[1-9]+(_[A, B])?.rds", recursive = T, full.names = T)
sample <- "GB8"

sampledir <- "/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/aligned/GB/P2_3/outs/per_sample_outs/"
data <- Read10X(data.dir=paste0(sampledir,sample,"/count/sample_filtered_feature_bc_matrix"))
data <- CreateSeuratObject(counts = data, project = sample, min.cells = 0, min.features = 500)

#Get mitochondrial Reads
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2

#Filter for UMI and mt-counts

mito.keep <- data$percent.mt <= 3
feat.keep <- !(isOutlier(data$nCount_RNA, nmads=3, type="higher"))

p1 <- ggplot() +
  geom_histogram(aes(data$nFeature_RNA), colour="black", fill="white", binwidth = 100) +
  geom_vline(aes(xintercept=500), color="blue", linetype="dashed") +
  labs(x="Detected features per cell")
p1
p2 <- ggplot() +
  geom_histogram(aes(data$nCount_RNA), colour="black", fill="white", binwidth = 1000) +
  geom_vline(aes(xintercept=200), color="blue", linetype="dashed") +
  geom_vline(aes(xintercept=max(data$nFeature_RNA[feat.keep])), color="blue", linetype="dashed") +
  labs(x="Total counts per cell")
p2
p3 <- ggplot() +
  geom_histogram(aes(data$percent.mt), colour="black", fill="white", binwidth = .5) +
  geom_vline(aes(xintercept=max(data$percent.mt[mito.keep])), color="blue", linetype="dashed") +
  labs(x="Percentage of mtRNA") +
  xlim(c(0,30))
p3

keep <- mito.keep & feat.keep
data <- data[, keep]

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2



### normalization and scaling with SCTransform ###############################################################################
print("Performing normalization, regression, scaling and dimensionality reduction...")
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F) # default options
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = F) # default options
data <- ScaleData(data, vars.to.regress = "nCount_RNA", verbose = F)
### dimensionality reduction
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = F)
### clustering and generating UMAP
data <- FindNeighbors(data, dims = 1:20, verbose = F)
data <- FindClusters(data, resolution = .6, verbose = F)
data <- RunUMAP(data, dims = 1:30, verbose = F)
DimPlot(data)

###Run Doublet Finder###

print("Doublet detection...")
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk <- as.numeric(levels(bcmvn$pK))[which.max(bcmvn$BCmetric)]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
exp_dbl <- (ncol(data)*.8/1000)/100 ## according to 10X's expected multiplet rate depnding on the no. of cells
nExp_poi <- round(exp_dbl*nrow(data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
data.dbl <- doubletFinder_v3(data, PCs = 1:15, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

## create output files -------------------------------------------------------------------------------------------------------
dbl.scores <- colnames(data.dbl@meta.data)[grepl("pANN",colnames(data.dbl@meta.data), fixed = T)]
data.dbl@meta.data[["Doublets"]] <- ifelse( data.dbl@meta.data[[dbl.scores]] <= 0.35, "Singlet", "Doublet")
#dbl.binary <- colnames(data.dbl@meta.data)[grepl("DF.classifications",colnames(data.dbl@meta.data), fixed = T)]
#DimPlot(data.dbl, group.by = dbl.binary) + ggtitle(paste0(table(data.dbl[[dbl.binary]])[[1]], " doublets detected"))

FeaturePlot(data.dbl, dbl.scores) + ggtitle("Doublet score") + scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits = c(0,0.35))
FeaturePlot(data.dbl, features = "nCount_RNA")+ scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits=c(2000,7000))
FeaturePlot(data.dbl, features = "nFeature_RNA")+ scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits=c(1000,5000))
FeaturePlot(data.dbl, features = "percent.mt")+ scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits=c(0,2))


data <- data[,which(data.dbl[["Doublets"]]=="Singlet")]
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = F)
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:15, verbose = F)
data <- FindClusters(data, resolution = .6, verbose = F)
data <- RunUMAP(data, dims = 1:15, verbose = F)
DimPlot(data, reduction = "umap", label = T)

#Using Feature vs RNA count matrix does not work...
#test <- data[,which(data[["seurat_clusters"]]=="0" | data[["seurat_clusters"]]=="1"| data[["seurat_clusters"]]=="6"|data[["seurat_clusters"]]=="11"| data[["seurat_clusters"]]=="3")]
#test <- data[,which(data[["seurat_clusters"]]=="0" | data[["seurat_clusters"]]=="1"| data[["seurat_clusters"]]=="6"|data[["seurat_clusters"]]=="11"| data[["seurat_clusters"]]=="3")]
#plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", split.by = "seurat_clusters", ncol = 4)
#plot2 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", ncol = 3)
#plot2

Cell_Types<- list(
  c("MBP","OLIG2","OLIG3","CLDN11","MOG","SOX10","GJB1","MAG","CNP","PLP1"),
  c("ENO2","NEFL","NEFM","NEFH","INA","SLC17A7","RBFOX3","SNAP25","SYT1","NEUROD2", "MAP2", "RBFOX1", "SNAP25"),
  c("NCAN", "NTSR2", "ADGRV1", "SLC1A2", "SLCO1C1", "ACSBG1", "GFAP", "MLC1", "PIRT", "SLC7A10", "AQP4"),
  c("TNR", "LHFPL3", "KLHL1", "PTPRZ1", "PLPP4", "FBN3", "MMP16", "TLL1", "COL11A1", "COL9A1"),
  c("FGF20", "P2RY12", "SMIM35", "ABCC4", "MUC22", "RASGEF1C"),
  c("RNASE2", "CEBPE", "CD163", "RETN", "CD207", "MSR1", "AIF1"),
  c("LCN6", "SELE", "CLDN5", "NOSTRIN", "ICAM2")
)

Cell_type_names <- c("Oligodendrocytes","Neurons", "Astrocytes", "OPC", "Microglia", "Macrophages", "Endothelia")
names(Cell_Types) <- Cell_type_names
data <- AddModuleScore(object = data, features= Cell_Types, name = Cell_type_names, assay ='RNA')
f1 <- function(x){
  paste0(x, seq_along(x))
}
l1<- f1(Cell_type_names)
FeaturePlot(data, features = unlist(l1))
FeaturePlot(data, features = c("LAMP5", "HTR2C", "NTNG2", "VIP", "KIT", "SST"))
DimPlot(data, reduction = "umap", label = T)

GB8<- data





##### #T4 ####

#sample.names <- dir("analysis_per_sample")
#sample.names <- dir("analysis_per_sample", pattern= "GB[1-9]+(_[A, B])?.rds", recursive = T, full.names = T)
#sample <- "T4"

sampledir <- "/omics/groups/OE0146/internal/Christina/10X_SC_RNA/FFPE/cellranger_align/GBM_MNG/C9/outs/per_sample_outs/T4_Pepsin/count/"
data <- Read10X(data.dir=paste0(sampledir,"/sample_filtered_feature_bc_matrix"))
data <- CreateSeuratObject(counts = data, project = sample, min.cells = 0, min.features = 500)

#Get mitochondrial Reads
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2

#Filter for UMI and mt-counts

mito.keep <- data$percent.mt <= 5
feat.keep <- !(isOutlier(data$nCount_RNA, nmads=3, type="higher"))

p1 <- ggplot() +
  geom_histogram(aes(data$nFeature_RNA), colour="black", fill="white", binwidth = 100) +
  geom_vline(aes(xintercept=500), color="blue", linetype="dashed") +
  labs(x="Detected features per cell")
p1
p2 <- ggplot() +
  geom_histogram(aes(data$nCount_RNA), colour="black", fill="white", binwidth = 1000) +
  geom_vline(aes(xintercept=200), color="blue", linetype="dashed") +
  geom_vline(aes(xintercept=max(data$nFeature_RNA[feat.keep])), color="blue", linetype="dashed") +
  labs(x="Total counts per cell")
p2
p3 <- ggplot() +
  geom_histogram(aes(data$percent.mt), colour="black", fill="white", binwidth = .5) +
  geom_vline(aes(xintercept=max(data$percent.mt[mito.keep])), color="blue", linetype="dashed") +
  labs(x="Percentage of mtRNA") +
  xlim(c(0,30))
p3

keep <- mito.keep & feat.keep
data <- data[, keep]

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2

data_raw <- data

### normalization and scaling with SCTransform ###############################################################################
print("Performing normalization, regression, scaling and dimensionality reduction...")
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F) # default options
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = F) # default options
data <- ScaleData(data, vars.to.regress = "nCount_RNA", verbose = F)
### dimensionality reduction
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = F)
### clustering and generating UMAP
data <- FindNeighbors(data, dims = 1:20, verbose = F)
data <- FindClusters(data, resolution = .6, verbose = F)
data <- RunUMAP(data, dims = 1:30, verbose = F)
DimPlot(data)

###Run Doublet Finder###

print("Doublet detection...")
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk <- as.numeric(levels(bcmvn$pK))[which.max(bcmvn$BCmetric)]

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
exp_dbl <- (ncol(data)*.8/1000)/100 ## according to 10X's expected multiplet rate depnding on the no. of cells
nExp_poi <- round(exp_dbl*nrow(data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
data.dbl <- doubletFinder_v3(data, PCs = 1:15, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

## create output files -------------------------------------------------------------------------------------------------------
dbl.scores <- colnames(data.dbl@meta.data)[grepl("pANN",colnames(data.dbl@meta.data), fixed = T)]
data.dbl@meta.data[["Doublets"]] <- ifelse( data.dbl@meta.data[[dbl.scores]] <= 0.35, "Singlet", "Doublet")
#dbl.binary <- colnames(data.dbl@meta.data)[grepl("DF.classifications",colnames(data.dbl@meta.data), fixed = T)]
#DimPlot(data.dbl, group.by = dbl.binary) + ggtitle(paste0(table(data.dbl[[dbl.binary]])[[1]], " doublets detected"))

FeaturePlot(data.dbl, dbl.scores) + ggtitle("Doublet score") + scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits = c(0,0.35))
FeaturePlot(data.dbl, features = "nCount_RNA")+ scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits=c(2000,7000))
FeaturePlot(data.dbl, features = "nFeature_RNA")+ scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits=c(1000,5000))
FeaturePlot(data.dbl, features = "percent.mt")+ scale_color_gradientn(colours = c("green","lightgrey"), na.value = "red", limits=c(0,2))


data <- data[,which(data.dbl[["Doublets"]]=="Singlet")]
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = F)
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:15, verbose = F)
data <- FindClusters(data, resolution = .6, verbose = F)
data <- RunUMAP(data, dims = 1:15, verbose = F)
DimPlot(data, reduction = "umap", label = T)

Cell_Types<- list(
  c("MBP","OLIG2","OLIG3","CLDN11","MOG","SOX10","GJB1","MAG","CNP","PLP1"),
  c("ENO2","NEFL","NEFM","NEFH","INA","SLC17A7","RBFOX3","SNAP25","SYT1","NEUROD2", "MAP2", "RBFOX1", "SNAP25"),
  c("NCAN", "NTSR2", "ADGRV1", "SLC1A2", "SLCO1C1", "ACSBG1", "GFAP", "MLC1", "PIRT", "SLC7A10", "AQP4"),
  c("TNR", "LHFPL3", "KLHL1", "PTPRZ1", "PLPP4", "FBN3", "MMP16", "TLL1", "COL11A1", "COL9A1"),
  c("FGF20", "P2RY12", "SMIM35", "ABCC4", "MUC22", "RASGEF1C"),
  c("RNASE2", "CEBPE", "CD163", "RETN", "CD207", "MSR1", "AIF1"),
  c("LCN6", "SELE", "CLDN5", "NOSTRIN", "ICAM2")
)

Cell_type_names <- c("Oligodendrocytes","Neurons", "Astrocytes", "OPC", "Microglia", "Macrophages", "Endothelia")
names(Cell_Types) <- Cell_type_names
data <- AddModuleScore(object = data, features= Cell_Types, name = Cell_type_names, assay ='RNA')
f1 <- function(x){
  paste0(x, seq_along(x))
}
l1<- f1(Cell_type_names)
FeaturePlot(data, features = unlist(l1))
FeaturePlot(data, features = c("LAMP5", "HTR2C", "NTNG2", "VIP", "KIT", "SST"))
DimPlot(data, reduction = "umap", label = T)

T4 <- data
