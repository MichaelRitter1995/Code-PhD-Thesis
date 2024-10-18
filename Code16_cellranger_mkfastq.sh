#copy sequencing data to /omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/raw_data/
#make sample_sheet.csv according to template
#copy sample_sheet.csv to /omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/raw_data/Illumina_folder_name

cellranger_mkfastq () {
if [ $# != 1 ]
  then
    echo "Wrong number of arguments supplied"
  else
    data_dir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/raw_data/$1
    outdir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/fastq

    cd $outdir

    bsub -q "long" -n 32 -R "rusage [mem=200G]" cellranger mkfastq --run=$data_dir \
                                                                       --csv=$data_dir/sample_sheet.csv \
                                                                       --localcores=32 \
                                                                       --localmem=200
fi
}

