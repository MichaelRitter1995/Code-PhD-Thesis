#Extract Visium and Xenium Sequences
#$1 is the pasta file of the Xenium probes

#Xen_query_seq=$(awk '0 == (NR) % 2' xenium_human_brain_gene_expression_panel_v1_probe_sequences.fasta | cat)
Xen_query_seq=$(awk '/^[^>]/' $1 | cat)
Vis_query_seq=$(awk -F ',' '{if (NR > 1) {print $2}}' Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv | cat)

#Get probes completely overlapping Visium Probes

for x in {0..5}
do
	#Count matching probe number
	no_match=0
	match=0
	
	#Prepare file for saving
	Binding=$((15+$x))

	echo "Visium alignment" $Binding "BP overlap"

        head -6 Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv > Revised_Probe_set_${Binding}_BP_Binding.csv
	
	#Circle through sequences
        for raw_seq in $Vis_query_seq
        do
		#remove first basepairs and get a 40 BP sequence
		a=$((10-$x))
		#get rid of the first space before each sequence
		b=$((1+$x))
                seq=${raw_seq:$a:40}
		#Extract sequence from Xenium probes
                Vis_seq=$(echo $Xen_query_seq | grep $seq)
		
		#Append probe file
                if [ -z "$Vis_seq" ]
                then
                        no_match=$(($no_match + 1))
                else
			echo $seq
                        match=$(($match + 1))
			Vis_seq=$(grep $seq Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv)
                        echo $Vis_seq >> Revised_Probe_set_${Binding}_BP_Binding.csv
			echo $Vis_seq >> 3_prime_end_seq.txt
        	fi

		#Repeat for 3' end -> remove last basepairs
                if [ $x != 5 ]
                then
                        seq=${raw_seq:$b:40}
                        Vis_seq=$(echo $Xen_query_seq | grep $seq)

                        if [ -z "$Vis_seq" ]
                        then
                                no_match=$(($no_match + 1))
                        else
				echo "3'"
				echo $seq
                                match=$(($match + 1))
				Vis_seq=$(grep $seq Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv)
                        	echo $Vis_seq >> Revised_Probe_set_${Binding}_BP_Binding.csv
				echo $Vis_seq >> 5_prime_end_seq.txt
                        fi
                fi
	done
        echo $Binding >> Summary.txt
        echo "Number of unique Sequences:" >> Summary.txt
        echo $no_match >> Summary.txt

        echo "Number of overlapping sequences:" >>Summary.txt
        echo $match >>Summary.txt
	echo " " >>Summary.txt

done


#Analyse Overlap of Xenium probes, this is done in the same way like it was done for the Visium probes

for x in {1..10}
do
	no_match=0
        match=0
	Binding=$((15-$x))
	echo $Binding
	head -6 Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv > Revised_Probe_set_${Binding}_BP_Binding.csv
	
	for raw_seq in $Xen_query_seq
	do
		
		z=$((40-$x))
		seq=${raw_seq:0:$z}
		Vis_seq=$(grep $seq Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv)

		if [ -z "$Vis_seq" ]
		then
     	 		no_match=$(($no_match + 1))
		else
      			match=$(($match + 1)) 
			echo $Vis_seq >> Revised_Probe_set_${Binding}_BP_Binding.csv
			echo $Vis_seq >> 3_prime_end_seq.txt
	fi

		if [ $x != 0 ]
		then
			seq=${raw_seq:$x:40}
      			Vis_seq=$(grep $seq Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv)
      	 		echo $seq >> .txt
			if [ -z "$Vis_seq" ]
      	 		then
        	       	 	no_match=$(($no_match + 1))
     			else
                		match=$(($match + 1))
                		echo $Vis_seq >> Revised_Probe_set_${Binding}_BP_Binding.csv
				echo $Vis_seq >> 5_prime_end_seq.txt
        		fi
		fi
	done

	#Remove duplicate lines
	awk '!x[$0]++'  Revised_Probe_set_${Binding}_BP_Binding.csv > temp.csv
	mv temp.csv Revised_Probe_set_${Binding}_BP_Binding.csv

	echo $x >> Summary.txt
	echo "Number of unique Sequences:" >> Summary.txt
	echo $no_match >> Summary.txt

	echo "Number of overlapping sequences:" >>Summary.txt
	echo $match >>Summary.txt
	echo " " >>Summary.txt

done


#Remove genes from previous probe set

for x in {5..19}

do
	#Remove lines present in other probe sets => each set is a unique set of probes
        mv Revised_Probe_set_${x}_BP_Binding.csv  temp.csv
        if      [ $x != 20 ]
        then
                y=$((x+1))
                head -6 Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv > Revised_Probe_set_${x}_BP_Binding.csv
                awk 'FNR==NR {a[$0]++; next} !($0 in a)' Revised_Probe_set_${y}_BP_Binding.csv temp.csv >> Revised_Probe_set_${x}_BP_Binding.csv
        fi
done

#Clean-up the probe set that is at the shift from Visium to Xenium probes for alignment

mv Revised_Probe_set_14_BP_Binding.csv temp.csv
for x in {15..20}

do
        cat Revised_Probe_set_${x}_BP_Binding.csv >> all.csv
done

head -6 Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv > Revised_Probe_set_14_BP_Binding.csv
awk 'FNR==NR {a[$0]++; next} !($0 in a)' all.csv temp.csv >> Revised_Probe_set_14_BP_Binding.csv

#Remove doublet rows from 5' and 3' files

awk '!x[$0]++'   5_prime_end_seq.txt > temp.csv
mv temp.csv  5_prime_end_seq.txt

awk '!x[$0]++'   3_prime_end_seq.txt > temp.csv
mv temp.csv  3_prime_end_seq.txt  

#Create a file containing all overlapping probes

for x in {5..20}

do
        cat Revised_Probe_set_${x}_BP_Binding.csv >> temp.csv
done

awk 'FNR==NR {a[$0]++; next} !($0 in a)' Ctrl_probes.csv temp.csv > All_Overlapping_probes.csv
awk '!x[$0]++'   temp.csv > All_probes.csv
rm temp.csv

#Create a final summary file

probes=0

for x in {5..20}
do
	count=$(cat Revised_Probe_set_${x}_BP_Binding.csv | wc -l)
	probes=$(($probes + $count-6))
done

echo "Total probes overallping between Visium and Xenium" >> Summary.txt
echo $probes >> Summary.txt


#Add control probes, needed to successfully run alignment

for x in {5..20}
do
cat Ctrl_probes.csv >> Revised_Probe_set_${x}_BP_Binding.csv
done

