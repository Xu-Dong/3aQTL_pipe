#!/bin/bash


# -- global settings
plink_path=/home/username/src
in_dir=/home/username/Project_XXX/input
out_dir=${in_dir}/E-GEUV-1-Genotype/plink_pca
curr_dir=`pwd`

# -- change work dir to the output dir
echo "Change to output directory..."
cd $out_dir
pwd

# -- transform vcf to plink format
# plink v1.9
vcf_path=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/input/E-GEUV-1-Genotype

date
echo "Convert vcf to plink format..."
for N in $(seq 1 22)
do
	${plink_path}/plink --vcf ${vcf_path}/GEUVADIS.chr${N}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz --const-fid --out chr${N}.plink &
	wait
# above command generate 3 files with ".bed", ".bim", and ".fam" suffix, code below write them to a file for further merging
	echo "chr${N}.plink.bed chr${N}.plink.bim chr${N}.plink.fam" >> merge_list.txt
done
sleep 10
# remove chr1 from merge_list.txt
cat merge_list.txt|tail -n+2 > tmp.txt
cat tmp.txt > merge_list.txt
rm tmp.txt

# -- merge all chromosomes
echo "merging genotype of all chromosomes into one file: merged_plink.* ..."
${plink_path}/plink --bfile chr1.plink --merge-list merge_list.txt --out merged_plink --allow-extra-chr &
wait

# -- split samples by sub-population: CEU, FIN, GBR, TSI, YRI
# get sample list in each sub-population
# format: one sample in a row with white-space separated two columns (FID IID), we use 0 for FID
cat ${in_dir}/subject_id2population.txt|tail -n+2|awk '$2=="CEU" {print '0',$1}' > CEU.sample
cat ${in_dir}/subject_id2population.txt|tail -n+2|awk '$2=="FIN" {print '0',$1}' > FIN.sample
cat ${in_dir}/subject_id2population.txt|tail -n+2|awk '$2=="GBR" {print '0',$1}' > GBR.sample
cat ${in_dir}/subject_id2population.txt|tail -n+2|awk '$2=="TSI" {print '0',$1}' > TSI.sample
cat ${in_dir}/subject_id2population.txt|tail -n+2|awk '$2=="YRI" {print '0',$1}' > YRI.sample

# extract population-specific samples and filtering
echo "split samples into sub-populations ..."
for i in *.sample
do
	${plink_path}/plink --bfile merged_plink --keep $i --geno 0.02 --hwe 0.000001 --maf 0.01 --make-bed --out ${i}_QC &
	wait
done
sleep 20

# -- pca analysis
echo "PCA analysis by plink..."
POPs=(CEU FIN GBR TSI YRI)
for (( i = 0 ; i < ${#POPs[@]} ; i++ ))
do
	echo "processing ${POPs[$i]}"
	${plink_path}/plink --bfile ${POPs[$i]}.sample_QC --indep-pairwise 50 5 0.2 --out ${POPs[$i]}.sample_QC &
	wait
	${plink_path}/plink --bfile ${POPs[$i]}.sample_QC --extract ${POPs[$i]}.sample_QC.prune.in --pca 30 --out ${POPs[$i]} &
	wait
done

echo "*.eigenvec files will be used as covariates"
echo "Done!"
date
