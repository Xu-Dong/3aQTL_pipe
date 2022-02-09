#!/bin/bash

sample_list=sample_list.txt
indir=${HOME}/Project_XXX/input
outdir=${HOME}/Project_XXX/output
# load python2
module load basic/anaconda2

for sample in `cat ${indir}/$sample_list`
do
	echo $sample
	samtools flagstat -@ 8 ${indir}/E-GEUV-1-RNA/${sample}.*.bam > ${indir}/E-GEUV-1-RNA/${sample}.flagstat &
	wait
done

sleep 30

python extract_read_depth.py --sample_list ${indir}/$sample_list --path_wig ${indir}/E-GEUV-1-RNA --output ${outdir}/wigFile_and_readDepth.txt &
wait
echo "Done!"
