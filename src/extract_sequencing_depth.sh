#!/bin/bash

# sample_list.txt contains all sample IDs (each line put one ID). The 445 sample IDs in Geuvadis RNA-seq project can be found in data/ of this repository
sample_list="sample_list.txt"
indir=${HOME}/Project_XXX/input
outdir=${HOME}/Project_XXX/output
src=${HOME}/Project_XXX/src
# load python2
module load basic/anaconda2

# step 1: apply samtools flagstat to obtain the statistic info of a bam
for sample in `cat ${indir}/$sample_list`
do
	echo $sample
	samtools flagstat -@ 8 ${indir}/E-GEUV-1-RNA/${sample}.*.bam > ${indir}/E-GEUV-1-RNA/${sample}.flagstat &
	wait
done

sleep 30

# step 2: extract read depth from step 1 output
python ${src}/extract_read_depth.py --sample_list ${indir}/$sample_list --path_wig ${indir}/E-GEUV-1-RNA --output ${outdir}/wigFile_and_readDepth.txt &
wait
echo "Done!"
