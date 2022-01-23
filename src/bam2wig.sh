#!/bin/bash
# use bedtools genomecov to convert bam files into bedgraph (*.wig)
# samples listed in a text file (each line contains one sample id) were processed one by one
# Alternative: to speed up, split samples into many small sets and processed in parallel.

echo "Start ..."
date
bamDir=${HOME}/Project/input/E-GEUV-1-RNA
selected_samples=${HOME}/Project/input/sample_list.txt

for sample in `cat $selected_samples|cut -f1`
do
	echo $sample
	bedtools genomecov -ibam ${bamDir}/${sample}.*.bam -bga -split -trackline > ${bamDir}/${sample}.wig
done

echo "Done!"
date
