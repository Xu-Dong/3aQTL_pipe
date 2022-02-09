#!/bin/bash

#prepare_susieR_uniqGene_location.py


date
base_in=/home/username/Project_XXX/matrix-eqtl
base_out=/home/username/Project_XXX/susieR_analysis
POPs=(CEU FIN GBR TSI YRI)
for (( i = 0; i < ${#POPs[@]}; i++ ))
do
	P=${POPs[$i]}
	echo ${P}

	python prepare_susieR_uniqGene_location.py --utr_loc_file ${base_in}/input/3utr_location.txt \
		--aQTL_map ${base_in}/output/${P}.cis_aQTL_all_control_gene_exprs.txt \
		--extend_size 1000000 \
		--outdir ${base_out}/input \
		--output ${P}.aGenes.loc_1Mb.txt &
done
wait
echo "Done!"
date
