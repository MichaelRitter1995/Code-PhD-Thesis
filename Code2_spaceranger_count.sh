#$1 - project name (folder with this name has to be created)
#$2 - flow cell ID
#Make sure there is a samp_config.csv file in the <project_name>/fastq/<flow cell ID> directory

if [ $# != 2 ]
  then
    echo "Wrong number of arguments supplied"
  else
    module load spaceranger/2.0.1

    prodir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium
    outdir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium/aligned/mouse
    datadir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium/fastq/$1
    probeset=/software/spaceranger/2.0.1/probe_sets/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv
    #probeset=/software/spaceranger/2.0.1/probe_sets/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv
    #probeset=/software/spaceranger/2.0.1/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv
    #probeset=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/reference_seq/costumized_probes/Xenium_CytAssist_probes/Revised_Probe_set.csv
    #probeset=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/reference_seq/costumized_probes/Xenium_CytAssist_probes/Revised_Probe_set_strict.csv
    #probeset=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/reference_seq/costumized_probes/Xenium_CytAssist_probes/Ctrl_Set.csv

    samples=$(cut -f2 -d, $datadir/outs/input_samplesheet.csv | tail -n $2)
    refpath=/omics/groups/OE0146/internal/shared/Ref/refdata-gex-GRCh38-2020-A

    cd $outdir
    for samp in $samples
    do  
       echo $samp
       bsub -q "long" -n 32 -R "rusage[mem=200G]" "spaceranger count --id=$samp \
                                                                         --fastqs=$datadir \
                                                                         --sample=$samp \
                                                                         --transcriptome=$refpath \
                                                                         --probe-set=$probeset \
                                                                         --cytaimage=$prodir/cytassist_images/$samp.tif \
                                                                         --image=$prodir/orig_images/$samp.tif \
                                                                         --unknown-slide=visium-2-large \
                                                                         --localcores=32 \
									 --no-bam \
                                                                         --localmem=200"

   done
fi

