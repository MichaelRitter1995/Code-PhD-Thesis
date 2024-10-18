#!/bin/bash
# 1. copy sequencing data to /omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/raw_data/
# 2. make sample_sheet.csv according to template
# 3. copy sample_sheet.csv to /omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/raw_data/Illumina_folder_name
# 4. Run spaceranger_mkfastq.sh and use the full name of the sequencing folder as input argument


if [ $# != 1 ]
  then
    echo "Wrong number of arguments supplied"
  else

    #data_dir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium/raw_data/$1
    #outdir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium/fastq

    #data_dir=/omics/groups/OE0146/internal/Helin/10x_snRNA/SC00012
    data_dir=/omics/groups/OE0146/internal/Christina/Visium_seq/sequencing/$1
    outdir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium/fastq


    cd $outdir

    bsub -q "long" -n 32 -R "rusage [mem=200G]" spaceranger mkfastq --run=$data_dir \
                                                                       --csv=$data_dir/sample_sheet.csv.save.1 \
                                                                       --localcores=32 \
                                                                       --localmem=200
fi


