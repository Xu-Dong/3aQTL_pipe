#!/bin/bash
# use bedtools genomecov to convert bam files into bedgraph (*.wig)
# before running this script make a text file (e.g. sample_list.txt) containing all samples (each line contains one sample id)
# Alternative: to speed up, split samples into many small sets and processed in parallel.

echo "Start ..."
date
# change bamDir to user defined location (where the downloaded bam files exist)
# sample_list.txt contains all sample IDs (each line put one ID). The 445 sample IDs in Geuvadis RNA-seq project can be found in data/ of this repository
bamDir=${HOME}/Project_XXX/input/E-GEUV-1-RNA
selected_samples=${HOME}/Project_XXX/input/sample_list.txt

for sample in `cat $selected_samples|cut -f1`
do
	echo $sample
	bedtools genomecov -ibam ${bamDir}/${sample}.*.bam -bga -split -trackline > ${bamDir}/${sample}.wig
done

echo "Done!"
date
