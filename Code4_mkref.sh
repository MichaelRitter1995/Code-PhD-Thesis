#!/bin/bash
# 
#


path=/omics/odcf/analysis/OE0146_projects/singlecell-and-spatialseq/reference_seq/costumized_probes/GB_Probes 

bsub -q "long" -n 32 -R "rusage [mem=200G]" spaceranger mkref --genome=GRCh38_probes \
                                                              --fasta=$path/fasta/genome_probes.fa \
                                                              --genes=$path/genes/genes_probes.gtf 
