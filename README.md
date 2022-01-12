# aQTL_pipe
aQTL pipeline applied on population scale genotype and transcriptomic data



* Desc: Mapping 3'aQTL on GEUV-1 population dataset, 445 individuals in total, which belong to 5 different sub-populations, RNA-seq data in LCL cell line.
* APA level was estimated by Dapars2, Association test was conducted on each sub-population separately.

## -- prior data prepared -- beyond the code below, the following data should be obtained prioring:
* RNA-seq alignment files in bam format (445 bam files, aligned to hg19), downloaded from Array-Express (https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/)
* Genotype data of the 445 individuals, VCF4.1 files, downloaded from Array-Express (https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/)
* Refseq gene annotation file (hg19_refseq_whole_gene.bed) and refseq transcript id map to gene symbol (hg19_refseq_id_to_symbol.txt) both files can be downloaded from UCSC table browser.

main(){
	> step 2
	run_bam2wig sample_list.txt
	run_generate_3utr_annotation
	run_flagstat_bam sample_list.txt
	run_generate_configure_dapars2

	> step 3
	run_dapars2 Dapars2_GEUV_all_samples.joint_configure.txt
	run_merge_dapars2_res
	
	> step 4
	run_genotype_pca
	run_peer
	run_prepare_input4matrixEQTL
	run_prepare_gene_and_snp_location

	> step 5
	run_matrix_eqtl

	> step 6
	run_susieR
}
