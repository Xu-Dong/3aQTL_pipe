#!/bin/bash
dir=${HOME}/Project/input/E-GEUV-1-Genotype
chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr18 chr19 chr20 chr21 chr22)
# ---- extract genotype from vcf file ----
for (( i = 0 ; i < ${#chrs[@]} ; i++ ))
do
        echo "Start vcftools on $CHR"
        CHR=${chrs[$i]}
        echo "extract genotype..."
        vcftools --gzvcf ${dir}/GEUVADIS.${CHR}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz --out ${dir}/GEUVADIS.${CHR}.gt.filtering --remove-filtered-all --maf 0.05 --max-missing-count 10 --extract-FORMAT-info GT &
        wait
        echo "extract freq info..."
        vcftools --gzvcf ${dir}/GEUVADIS.${CHR}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz --out ${dir}/GEUVADIS.${CHR}.gt.filtering --remove-filtered-all --maf 0.05 --max-missing-count 10 --freq &
        wait
        echo "$CHR Done!"

# ---- encode genotype as 012 ----
for (( i = 0 ; i < ${#chrs[@]} ; i++ ))
do
        python recode_with_012.py --frq GEUVADIS.${chrs[$i]}.gt.filtering.frq --GT GEUVADIS.${chrs[$i]}.gt.filtering.GT.FORMAT --output GEUVADIS.${chrs[$i]}.gt012.bed &
        wait
done
sleep 30
# ---- merge all chromosomes into one file ----
cat GEUVADIS.*.gt012.bed  > GEUVADIS.all_chrs.gt012.bed
