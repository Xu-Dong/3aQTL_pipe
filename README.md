# Mapping 3′UTR alternative polyadenylation quantitative trait loci through population-scale transcriptomic and genomic data

**Abbreviation** 
* APA: alternative polyadenylation
* 3'aQTL: 3′UTR alternative polyadenylation quantitative trait loci

This pipeline describes the detailed steps for analyzing dynamics alternative polyadenylation events across lymphoblastoid cell lines (LCL) samples from 445 unrelated individuals and performing association analysis between common genetic variants and APA events to obtain a map of genetic regulation of APA. The whole pipeline includes APA quantitative analysis across samples, association test between common genetic variants and APA usage (mapping 3'aQTL), fine-mapping 3'aQTLs, and other steps for preparing phenotype, genotype data and processing output of above analyses. The final outputs of this pipeline including the matrix of APA usage profile across samples, the table of association between common genetic variants and APA usage (3'aQTLs), the table of fine-mapped 3'aQTLs.

The scripts in this repository were used to analyze Geuvadis RNA-seq Project dataset which contains RNA sequencing data of 445 unrelated individuals and corresponding genotype data from 1000 Genome Project. The scripts with a few edits can be reused to perform APA analysis and 3'aQTL mapping on other data source.

For conditions to reuse of these scripts please refer to LICENSE file.

## Using this pipeline
Details on how to prepare environment and use the script can be found on [GitHub wiki](https://github.com/Xu-Dong/3aQTL_pipe/wiki) pages for this repository.

