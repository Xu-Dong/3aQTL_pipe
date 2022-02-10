# Mapping 3′UTR alternative polyadenylation quantitative trait loci through population-scale transcriptomic and genomic data

aQTL pipeline applied on population scale genotype and transcriptomic data



* Desc: Mapping 3'aQTL on GEUV-1 population dataset, 445 individuals in total, which belong to 5 different sub-populations, RNA-seq data in LCL cell line.
* APA level was estimated by Dapars2, Association test was conducted on each sub-population separately.


## Before starting the pipeline the following data should be obtained:
* RNA-seq alignment files in bam format (445 bam files, aligned to hg19), downloaded from Array-Express (https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/)
* Genotype data of the 445 individuals, VCF4.1 files, downloaded from Array-Express (https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/)
* Refseq gene annotation file (hg19_refseq_whole_gene.bed) and refseq transcript id map to gene symbol (hg19_refseq_id_to_symbol.txt) both files can be downloaded from UCSC table browser.

## Make a workspace (root directory) for this project
***NOTE:*** We use pseudo-path in all the codes in src/. Users apply the codes should modify the pathes to their own. We recommend making a root directory for
testing the whole project.
* Use the command below to build a workspace
```
> cd $HOME # change to home directory
> mkdir Project_XXX # change XXX to user defined project name
> cd Project_XXX # change path to the built root directory
> mkdir input output src matrix-eqtl susieR_analysis
```
