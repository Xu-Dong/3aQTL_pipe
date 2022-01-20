#!/bin/bash
# Split Genotype data and Phenotype data into individual gene level
# Usage: bash prepare_susie_input_02.sh CEU

# Global settings:
src=${HOME}/Projects/3aQTL_Pipe/src
PDUI_DIR=${HOME}/Projects/3aQTL_Pipe/output/matrix-eqtl/input
susie_DIR=${HOME}/Projects/3aQTL_Pipe/output/susieR_analysis

POP=$1

cd $susie_DIR
mkdir -p $POP


# --- preparing phenotype data of individual gene: format {sample_id sample_id PDUI}
cd ${susie_DIR}/${POP}
for gene in `cat ${susie_DIR}/input/${POP}.aGenes.loc_1Mb.txt|cut -f1`
do
        mkdir -p $gene
        cd $gene
        cat ${PDUI_DIR}/PDUI_mat.${POP}.txt | awk -v aGENE=$gene -F"\t" 'BEGIN{OFS="\t"} {if (NR==1) { for (i=2;i<NF;++i) SAMPLES[i]=$i} if ($1==aGENE){ for(i=2;i<NF;++i) print SAMPLES[i],SAMPLES[i],$i}}' > expr.phen &
        wait
        cd ${susie_DIR}/${POP}
done

# --- preparing genotype of SNPs within 1Mb region of selected gene

cd ${susie_DIR}/${POP}
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
        cat ${susie_DIR}/input/Header.${POP} > 3aQTL.vcf
        bedtools intersect -a ${susie_DIR}/input/Genotype_mat.${POP}.bed -b gene_loc.bed -wa |cut -f4- >> 3aQTL.vcf &
        wait
        rm gene_loc.bed
        cd ${susie_DIR}/${POP}

done < ${susie_DIR}/input/${POP}.aGenes.loc_1Mb.txt

