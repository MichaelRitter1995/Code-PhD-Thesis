#$1 - flow cell ID
#$2 - number of samples in flow cell

if [ $# != 2 ]
  then
    echo "Wrong number of arguments supplied"
  else
    data_dir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium/fastq/$1
    visium_dir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium

    samples=$(cut -f2 -d, $data_dir/outs/input_samplesheet.csv | tail -n $2)

    refpath=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/reference_seq/refdata-gex-GRCh38-2020-A

	for x in {5..20}
	do
	cd /omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium/aligned/Xenium/Probe_Binding/GB_panel
		if [ -d ${x}_BP_Binding/ ]
                then
                        echo "Folder already exists"
                else
			mkdir ${x}_BP_Binding
		fi
		
		probepath=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/reference_seq/costumized_probes/Xenium_CytAssist_probes/GB_panel/Revised_Probe_set_${x}_BP_Binding.csv
		outdir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium/aligned/Xenium/Probe_Binding/GB_panel/${x}_BP_Binding
		cd $outdir

    		for samp in $samples
        	do
     
    			if [ -d $2/ ]
      	  		then
          			echo "Error: Sample Directory already exists"
          			exit
        		else
          			echo $samp

				bsub -q "long" -n 32 -R "rusage[mem=200G]" "spaceranger count --id=$samp \
                                		                                              --fastqs=$data_dir/outs/fastq_path/$1 \
                                                     			                      --sample=$samp \
                                                  		                              --transcriptome=$refpath \
                                                                    			      --probe-set=$probepath \
                                                                               	              --image=$visium_dir/orig_images/${samp}.tif \
							            			      --cytaimage=$visium_dir/cytassist_images/${samp}.tif \
											      --no-bam \
                                                                         		      --unknown-slide=visium-2-large \
									 		      --localcores=32 \
											      --loupe-alignment=$visium_dir/manual_alignment/${samp}.json \
                                                                             	              --localmem=200"
			fi
		done
	done
fi


cd /omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium/aligned/Xenium/Probe_Binding/GB_panel
if [ -d 0_BP_Binding/ ]
then
	echo "Folder already exists"
else
        mkdir 0_BP_Binding
fi

probepath=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/reference_seq/costumized_probes/Xenium_CytAssist_probes/GB_Panel/Revised_Probe_set_0_BP_Binding.csv
outdir=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/visium/aligned/Xenium/Probe_Binding/GB_panel/0_BP_Binding
cd $outdir

samples=$(cut -f2 -d, $data_dir/outs/input_samplesheet.csv | tail -n $2)

for samp in $samples
do

	if [ -d $2/ ]
        then
		echo "Error: Sample Directory already exists"
                exit
        else
           	echo $samp
		bsub -q "long" -n 32 -R "rusage[mem=200G]" "spaceranger count --id=$samp \
                                                                              --fastqs=$data_dir/outs/fastq_path/$1 \     
                                                                              --sample=$samp \
                                                                              --transcriptome=$refpath \
                                                                              --probe-set=$probepath \
                                                                              --image=$visium_dir/orig_images/${samp}.tif \
                                                                              --cytaimage=$visium_dir/cytassist_images/${samp}.tif \
                                                                              --no-bam \
                                                                              --unknown-slide=visium-2-large \
                                                                              --localcores=32 \
									      --loupe-alignment=$visium_dir/manual_alignment/${samp}.json \
                                                                              --localmem=200"
	fi
done
