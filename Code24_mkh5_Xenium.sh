#Upload Xenium data to the cluster and rename main folder
# start script and use folder name and project name as argument
# $1 = project name
# $2 = folder name

if [ $# != 2 ]
  then
    echo "Wrong number of arguments supplied"
fi

module load R/4.3.0

RAW=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Xenium/Raw_data/$1/$2
cd $RAW

#Extract fodler and sample ID
samples=$(ls | tail -12 | cat)

#Go to the output directory
PROCESSED=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Xenium/Processed_data
cd $PROCESSED

for samp in $samples
        do
	ID=$(echo $samp | awk -F '__' '{ print $3 }')
        echo $ID
        cd $PROCESSED/$1/$2/$ID

        echo "Convert MEX to HDF5"
        Rscript /omics/groups/OE0146/internal/Micha/Unlimited_Power/R/functions/MEX_to_HDF5.R $1 $2 $samp $ID

done

