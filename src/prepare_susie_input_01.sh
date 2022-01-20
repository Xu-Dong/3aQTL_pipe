#!/bin/bash
# Preparing input data for running fine-mapping with susieR
# Two in-house scripts: prepare_susieR_uniqGene_location.py; genotype_2_bed.py

# global settings 1: Path of in-house scripts, path of matrix-eqtl results, path for susieR input
src=${HOME}/Projects/3aQTL_Pipe/src
base_in=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/output/matrix-eqtl
base_out=/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/output/susieR_analysis

# global settings 2: path of genotype data, path of data for susieR input       
GT=${HOME}/Projects/3aQTL_Pipe/output/matrix-eqtl/input
susie_IN=${HOME}/Projects/3aQTL_Pipe/output/susieR_analysis/input

module load basic/anaconda2
POPs=(CEU FIN GBR TSI YRI)

# --- Preparing all aGenes (APA genes with significant genetic associations) and their 3'UTR location (extend 1000000bp at both sides)
echo "Preparing aGenes..."
for (( i = 0; i < ${#POPs[@]}; i++ ))
do
        P=${POPs[$i]}
        echo ${P}

        python ${src}prepare_susieR_uniqGene_location.py --utr_loc_file ${base_in}/input/3utr_location.txt \
                --aQTL_map ${base_in}/output/${P}.cis_aQTL_all_control_gene_exprs.txt \
                --extend_size 1000000 \
                --outdir ${base_out}/input \
                --output ${P}.aGenes.loc_1Mb.txt &
done
wait



# --- Preprocessing genotype data: transform the genotype (012 format) data into bed file which will be further used by bedtools
echo "Preparing genotype..."
for (( i = 0 ; i < ${#POPs[@]} ; i++ ))
do
        P=${POPs[$i]}
        python ${src}/genotype_2_bed.py ${GT}/Genotype_mat.${P}.txt ${susie_IN}/Genotype_mat.${P}.bed ${susie_IN}/Header.$P &
        wait
done

echo "Done!"