#$3 - project (check with "projects" command)
#Create a samp_config_{your_sample_name}.csv file for each sample
#"gene-expression" and libraries" must be in there always
#"samples" only for multiplexed or multibarcode samples
#Make sure there is a samp_config.csv file in the /omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/fastq/<flow cell ID> directory
# $1=flow cell
# $2=sample
# $3=project


if [ $# != 3 ]
  then
    echo "Wrong number of arguments supplied"
  else
    datadir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/fastq/$1
    outdir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/singlecell/aligned/$3

    samples=$(cut -f2 -d, $datadir/outs/input_samplesheet.csv | tail -n $2)

    cd $outdir
    #for samp in $samples
     #do
      echo $2
      if [ -d $2/ ]
       then
	echo "Error: Sample Directory already exists"
       else
	bsub -q "verylong" -n 32 -R "rusage[mem=200G]" "cellranger multi --id=$2 \
                                                                     --csv=$datadir/samp_config_${2}.csv \
                                                                     --localcores=32 \
                                                                     --localmem=200"
    fi
  #done
fi

