#Upload Xenium data to the cluster and rename main folder
# start script and use folder name and project name as argument
# $1 = project name
# $2 = folder name

if [ $# != 2 ]
  then
    echo "Wrong number of arguments supplied"
fi

#Activate conda environment
#conda activate xenium
RAW=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Xenium/Raw_data/$1/$2
cd $RAW

#Extract fodler and sample ID
samples=$(ls | tail -12 | cat)

#Go to the output directory
PROCESSED=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/Xenium/Processed_data
cd $PROCESSED

#Check if Run folder already exists
if [ -d $PROCESSED/$1/$2/ ]
	then
		echo "Directory alread exists"
	else
		mkdir $1/$2
	fi

for samp in $samples
	do
	ID=$(echo $samp | awk -F '__' '{ print $3 }')
	echo $ID
#Still needs to be worked on; should ask user, if an direcrtory, already existing for this sample should be overwritten
#	if [ -d ../Processed_data/$1/$samp/ ]
#        then
#          read -p "Project directory already exists, overwrite? (y/n) " yn
#		case $yn in 
#			y ) echo overwrite data;;
#			n ) echo exiting...;
#				exit;;
#			* ) echo invalid response;
#			exit 1;;
#		esac
#		
#	else

#Create sample directrory

	if [ -d $PROCESSED/$1/$2/$ID ]
        then
		echo "Directory alread exists"
	else
                mkdir $PROCESSED/$1/$2/$ID
        fi
	
	cd $PROCESSED/$1/$2/$ID

#Old Xenium Analyzer output
#	echo "unzip transcripts.csv"
#	gunzip $RAW/$samp/transcripts.csv.gz

#New Xenium Analyzer output
	#python $ANALYSIS/Unlimited_Power/Python/parquet_to_csv.py  -file $RAW/$samp/transcripts.parquet -out $RAW/$samp/transcripts.csv
	gunzip $RAW/$samp/cells.csv.gz

	echo "Filter Transcripts"
	python $ANALYSIS/Unlimited_Power/Python/filter_transcripts.py -transcript $RAW/$samp/transcripts.parquet
	echo "Map Transcripts"

	python $ANALYSIS/Unlimited_Power/Python/map_transcripts.py  -out nuclei_mtx

	echo "Convert MEX to HDF5"
	Rscript /omics/groups/OE0146/internal/Micha/Unlimited_Power/R/functions/MEX_to_HDF5.R $1 $2 $samp $ID

done
