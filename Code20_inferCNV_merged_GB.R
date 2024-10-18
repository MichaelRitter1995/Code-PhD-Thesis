#!/usr/bin/env Rscript
#Submit this script with the following command - make sure to laod the R module:
#bsub -q "verylong" -n 32 -R "rusage[mem=200G]" "./inferCNV_merged_GB.R"

library(Seurat)
library(rjags)
library(infercnv)


setwd("/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Micha/Combined_Analysis/sn/analysis_merged/")

annot.path <- "CNVcalling/annot/"
out.path <- "CNVcalling/outs/"



### create annotation sets ###################################################################################################
samp <-'GB_FFPE_MP'
data <- readRDS("harmon_strict_DR_Annot.rds")
ref <- readRDS("Reference_CNV_inference.rds")
#merge reference and data
data <- merge(ref, data)
data <- SCTransform(data, vars.to.regress = c("nCount_RNA"), verbose = F)

dir.create(annot.path, recursive = T)
dir.create(out.path, recursive = T)

#extracts a table with the raw counts
write.table(as.matrix(data@assays$SCT@counts), paste0(annot.path,samp,"_counts.txt"), sep = "\t", col.names = T, row.names = T, quote = F)

annot <- data.frame(colnames(data), data$Annot_v1)

#Rename the Unknown cells according to sample, in this step, the groups used for the output are selected

#annot$data.Annot_v1[which(annot$data.Annot_v1== "Unknown")] <- data$orig.ident[which(annot$data.Annot_v1== "Unknown")]
#annot$data.annot_v1[is.na(annot$data.annot_v1)] <- data$orig.ident[is.na(annot$data.annot_v1)]

write.table(annot, paste0(annot.path,samp,"_annot.txt"), sep = "\t", col.names = F, row.names = F, quote = F)


### create inferCNV object ###################################################################################################

#print(samp)

### addition ###
annot <- read.table(paste0(annot.path, samp, "_annot.txt"), sep = "\t")
#ref <- unique(annot[,2])[!grepl("^([0-9])", unique(annot[,2]))]
ref <- unique(annot[,2])[grepl("Ref_", unique(annot[,2]))]
#ref <- ref[!ref %in% c("Tumor")]
#ref <- ref[which(ref %in% annot[,2])]
###

infcnv <- CreateInfercnvObject(raw_counts_matrix = paste0(annot.path,samp,"_counts.txt"),
                               annotations_file = paste0(annot.path,samp,"_annot.txt"),
                               gene_order_file = '/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Dilan/FFPEvsCryo/singelCell/codeshg38_gencode_v27.txt',
                               ref_group_names = ref)


infcnv <- infercnv::run(infcnv, cutoff = 0.1, 
                        out_dir = paste0(out.path,samp,"_i6_wind201"),
                        denoise = T, HMM = T, HMM_type='i6', cluster_by_groups = T, num_threads = 8, window_length = 201,
                        cluster_references = T, analysis_mode = "samples")

infercnv.obj.medianfilt = infercnv::apply_median_filtering(infcnv)

infercnv::plot_cnv(infercnv.obj.medianfilt,
                   out_dir = paste0(out.path,samp,"_i6_wind201"),
                   output_filename = 'infercnv.median_filtered',
                   x.range = "auto",
                   x.center = 1,
                   title = "infercnv",
                   output_format = "pdf")

quit(save="no")
