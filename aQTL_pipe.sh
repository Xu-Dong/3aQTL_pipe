#!/bin/bash

#Desc: Mapping 3'aQTL on GEUV-1 population dataset, 445 individuals in total, which belong to 5 different sub-populations, RNA-seq data in LCL cell line.
# APA level was estimated by Dapars2, Association test was conducted on each sub-population separately.
#
# -- prior data prepared -- beyond the code below, the following data should be obtained prioring:
# RNA-seq alignment files in bam format (445 bam files, aligned to hg19), downloaded from Array-Express (https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/)
# Genotype data of the 445 individuals, VCF4.1 files, downloaded from Array-Express (https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/)
# Refseq gene annotation file (hg19_refseq_whole_gene.bed) and refseq transcript id map to gene symbol (hg19_refseq_id_to_symbol.txt) both files can be downloaded from UCSC table browser.

main(){
	# step 2
	run_bam2wig sample_list.txt
	run_generate_3utr_annotation
	run_flagstat_bam sample_list.txt
	run_generate_configure_dapars2

	# step 3
	run_dapars2 Dapars2_GEUV_all_samples.joint_configure.txt
	run_merge_dapars2_res
	
	# step 4
	run_genotype_pca
	run_peer
	run_prepare_input4matrixEQTL
	run_prepare_gene_and_snp_location

	# step 5
	run_matrix_eqtl

	# step 6
	run_prepare_data_for_susieR_01
	run_prepare_data_for_susieR_02
	run_susieR CEU
	run_susieR FIN
	run_susieR GBR
	run_susieR TSI
	run_susieR YRI
}


# -- functions -- list order: bottom to up ---

# run susieR
function run_susieR(){
	src=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/src
	susie_DIR=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/output/susieR_analysis
	curr_dir=`pwd`
	POP=$1
	module load languages/R-3.6.3
	date
	cd ${susie_DIR}/$POP
	for gene in `cat ${susie_DIR}/input/${POP}.aGenes.loc_1Mb.txt|cut -f1`
	do
		if [ -d $gene ]
		then
			cd $gene
			if [ -f "3aQTL.vcf" -a -f "expr.phen" ]
			then
				Rscript ${src}/SuSiE_GTEx.r
			else
				continue
			fi
		else
			continue
		fi
		cd ${susie_DIR}/$POP
	done

	date

}

# prepare step 2: prepare 3aQTL.vcf and expr.phen for each selected aGene
function run_prepare_data_for_susieR_02(){
	src=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/src
	PDUI_DIR=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/output/matrix-eqtl/input
	susie_DIR=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/output/susieR_analysis
	curr_dir=`pwd`
	POPs=(CEU FIN GBR TSI YRI)

	# mkdir
#	cd $susie_DIR
#	mkdir -p CEU
#	mkdir -p FIN
#	mkdir -p GBR
#      	mkdir -p TSI
#	mkdir -p YRI
#
#	for (( i = 0 ; i < ${#POPs[@]} ; i++ ))
#	do
#		cd ${susie_DIR}/${POPs[$i]}
#		for gene in `cat ${susie_DIR}/input/${POPs[$i]}.aGenes.loc_1Mb.txt|cut -f1`
#		do
#			mkdir -p $gene
#			cd $gene
#			cat ${PDUI_DIR}/PDUI_mat.${POPs[$i]}.txt | awk -v aGENE=$gene -F"\t" 'BEGIN{OFS="\t"} {if (NR==1) { for (i=2;i<NF;++i) SAMPLES[i]=$i} if ($1==aGENE){ for(i=2;i<NF;++i) print SAMPLES[i],SAMPLES[i],$i}}' > expr.phen &
#			wait
#			cd ${susie_DIR}/${POPs[$i]}
#		done
#	done

	for (( i = 0 ; i < ${#POPs[@]} ; i++ ))
	do
		cd ${susie_DIR}/${POPs[$i]}
		while read line
		do
			gene=`echo $line|awk '{print $1}'`
			loc=`echo $line|awk '{print $2}'`
			cd $gene
			CHR=${loc%:*}
			COORD=${loc#*:}
			S=${COORD%-*}
			E=${COORD#*-}
			echo -e "$CHR\t$S\t$E" > gene_loc.bed
			cat ${susie_DIR}/input/Header.${POPs[$i]} > 3aQTL.vcf 
			bedtools intersect -a ${susie_DIR}/input/Genotype_mat.${POPs[$i]}.bed -b gene_loc.bed -wa |cut -f4- >> 3aQTL.vcf &
			wait
			rm gene_loc.bed
			cd ${susie_DIR}/${POPs[$i]}

		done < ${susie_DIR}/input/${POPs[$i]}.aGenes.loc_1Mb.txt

	done
}

# prepare step 1: prepare aGenes, locations, genotype formats
function run_prepare_data_for_susieR_01(){
	src=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/src
	GT=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/output/matrix-eqtl/input
	susie_IN=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/output/susieR_analysis/input

	module load basic/anaconda2
	# prepare the extended region of 3'UTR of each selected aGene,five txt files will be generated (*.aGenes.loc_1Mb.txt, * could be one of CEU/FIN/GBR/TSI/YRI)
	bash ${src}/submit_prepare_susieR_uniqGene_loc.sh &
	wait
	# transform the genotype 012 format data into bed file
	python ${src}/genotype_2_bed.py ${GT}/Genotype_mat.CEU.txt ${susie_IN}/Genotype_mat.CEU.bed ${susie_IN}/Header.CEU &
	wait
	python ${src}/genotype_2_bed.py ${GT}/Genotype_mat.FIN.txt ${susie_IN}/Genotype_mat.FIN.bed ${susie_IN}/Header.FIN &
	wait
	python ${src}/genotype_2_bed.py ${GT}/Genotype_mat.GBR.txt ${susie_IN}/Genotype_mat.GBR.bed ${susie_IN}/Header.GBR &
	wait
	python ${src}/genotype_2_bed.py ${GT}/Genotype_mat.TSI.txt ${susie_IN}/Genotype_mat.TSI.bed ${susie_IN}/Header.TSI &
	wait
	python ${src}/genotype_2_bed.py ${GT}/Genotype_mat.YRI.txt ${susie_IN}/Genotype_mat.YRI.bed ${susie_IN}/Header.YRI &
	wait

}

function run_matrix_eqtl(){
	POPs=(CEU FIN GBR TSI YRI)
	curr_dir=`pwd`
	for (( i = 0 ; i < ${#POPs[@]} ; i++ ))
	do
		POP=${POPs[$i]}
		echo "
#!/bin/bash
#PBS -N run_matrixEQTL_$POP
#PBS -q cu-1
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -S /bin/bash" > ${curr_dir}/submit_matrixEQTL_${POP}.sh
		echo '
cd $PBS_O_WORKDIR' >> ${curr_dir}/submit_matrixEQTL_${POP}.sh
		echo "
module load languages/R-3.6.3
SRC=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/src
population=$POP" >> ${curr_dir}/submit_matrixEQTL_${POP}.sh
		echo '
Rscript ${SRC}/run_Matrix_eQTL.R $population' >> ${curr_dir}/submit_matrixEQTL_${POP}.sh &
		wait
		qsub ${curr_dir}/submit_matrixEQTL_${POP}.sh &
		wait
	done
}

# 
function run_prepare_gene_and_snp_location(){
	src=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/src
	module load basic/anaconda2
	basedir=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project
	indir={basedir}/output
	outdir=${basedir}/output/matrix-eqtl/input

	echo "Extracting 3' UTR location ..."
	python extract_3UTR_location.py --dapars_res ${indir}/Dapars2_geuv_res.all_chromosomes.txt --output ${outdir}/3utr_location.txt &
	wait
	echo "Extracting SNP location ..."
	python extract_SNP_location.py --genotype_bed ${basedir}/input/E-GEUV-1-Genotype/GEUVADIS.all_chrs.gt012.bed --output ${outdir}/snp_location.txt &
	wait
	echo "Done!"

}
# prepare pdui matrix, covariates table, genotype matrix,all files with consistent sample order
function run_prepare_input4matrixEQTL(){
	src=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/src
	module load languates/R-3.6.3
	Rscript ${src}/prepare_input_files_matrixeqtl.R
}

# run impute and peer
function run_peer(){
	src=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/src
	module load languages/R-3.6.3
	echo "Start ..."
	date
	Rscript ${src}/3UTR_impute_peer.R &
	wait
	echo "Done!"
	date
}
# plink 1.9 is required
function run_genotype_pca(){
	src=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/src
	bash ${src}/merge_vcf_and_pca_by_plink.sh
}
# this is optional,only required if we ran dapars2 on separated chromosomes
# change the "dir" option if required
# change the input of "subject_id2population.txt" which only provides the sample list in the same order with
# the header of dapars2 result table.
# this will generated a file namely "Dapars2_geuv_res.all_chromosomes.txt" in current directory
function run_merge_dapars2_res(){
	module load language/R-4.1.0
	Rscript merge_dapars2_res_by_chr.R
}

function run_dapars2(){
	CHR_LIST=(chrX chrY chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)
	CONFIG=$1
	curr_dir=`pwd`
	for (( i = 0 ; i < ${#CHR_LIST[@]} ; i++ ))
	do
		CHR=${CHR_LIST[$i]}
		echo "
#!/bin/bash
#PBS -N run_dapars2_$CHR
#PBS -q cu-1
#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -S /bin/bash" > ${curr_dir}/submit_dapars2_${CHR}.sh
		echo '
cd $PBS_O_WORKDIR' >> ${curr_dir}/submit_dapars2_${CHR}.sh
		echo "
module load languages/R-3.6.3
module load basic/anaconda2
dapars2_path=/home/zouxudong/src/DaPars2-master/src
configure_file=$CONFIG
chromosome=$CHR" >> ${curr_dir}/submit_dapars2_${CHR}.sh
		echo '
python ${dapars2_path}/Dapars2_Multi_Sample.py ${configure_file} ${chromosome}' >> ${curr_dir}/submit_dapars2_${CHR}.sh &
		wait
		qsub ${curr_dir}/submit_dapars2_${CHR}.sh &
		wait
		
	done
}

# generate a configure file required by running dapars2
# wig files were generated by "run_bam2wig"
# the summary file contains all wig files and read depth were generated by "run_flagstat_bam"
# the reference of 3'UTR region were generated by "run_generate_3utr_annotation"
function run_generate_configure_dapars2(){
	src=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/src
	indir=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/output

	module load basic/anaconda2
	# use relative path for "--out_dir"
	python ${src}/generate_configure_for_Dapars2.py --annotation_3utr ${indir}/hg19_refseq_3utr_annotation.bed \
		--wigFile_depth ${indir}/wigFile_and_readDepth.txt \
		--coverage_threshold 15 \
		--threads 8 \
		--out_dir output/Dapars2_Output \
		--out_prefix Dapars2_geuv_res \
		--out_config_name Dapars2_GEUV_all_samples.joint_configure.txt
}

function run_flagstat_bam(){
	sample_list=$1
	indir=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/input
	outdir=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/output
	module load basic/anaconda2

	for sample in `cat ${indir}/$sample_list`
	do
		echo $sample
		samtools flagstat -@ 8 ${indir}/E-GEUV-1-RNA/download/${sample}.*.bam > ${indir}/E-GEUV-1-RNA/download/${sample}.flagstat &
		wait
	done
	sleep 30


	python extract_read_depth.py --sample_list ${indir}/$sample_list --path_wig ${indir}/E-GEUV-1-RNA/download --output ${outdir}/wigFile_and_readDepth.txt &
	wait
	echo "Done!"
}

function run_generate_3utr_annotation(){
	whole_gene_bed="/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/input/hg19_refseq_whole_gene.bed"
	refid2symbol="/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/input/hg19_refseq_transcript2geneSymbol.txt"
	ref_3utr_bed="/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/output/hg19_refseq_3utr_annotation.bed"
	path_dapars2_source="/home/zouxudong/src/DaPars2-master/src"

	python ${path_dapars2_source}/DaPars_Extract_Anno.py -b $whole_gene_bed -s $refid2symbol -o $ref_3utr_bed
}

function run_bam2wig(){
	sample_list=$1
	SRC=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/src
	
	bash ${SRC}/bam2wig.sh &
	wait

}


# -- main


main
