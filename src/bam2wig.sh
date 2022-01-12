#!/bin/bash

echo "Start ..."
date
bamDir=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/input/E-GEUV-1-RNA/download
selected_samples=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/remained_21

for sample in `cat $selected_samples|cut -f1`
do
	echo $sample
	bedtools genomecov -ibam ${bamDir}/${sample}.*.bam -bga -split -trackline > ${bamDir}/${sample}.wig
done

echo "Done!"
date
