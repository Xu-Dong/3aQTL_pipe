#!/opt/app/languages/R-3.6.3/bin/Rscript

library(MatrixEQTL)

usage <- function(){
	cat(
	    'usage: run_Matrix_eQTL.R <sub_population>')
}

args <- commandArgs(trailingOnly=TRUE)
POP <- args[1] #CEU/FIN/GBR/TSI/YRI

# -- check input args
if(is.na(POP)) {
	usage()
	quit(save='no',status=1)
}


# -- settings
DIR <- "/home/username/Project_XXX/matrix-eqtl/"
CIS_DISTANCE <- 1e6

# - Use linear model
useModel = modelLINEAR # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# - Genotype file name
SNP_file_name = paste0(DIR,"input/","Genotype_mat.",POP,".txt")
snps_location_file_name = paste0(DIR,"input/","snp_location.txt")

# - APA expression file name
expression_file_name = paste0(DIR,"input/","PDUI_mat.",POP,"imputed_qnorm.txt")
gene_location_file_name = paste0(DIR,"input/","3utr_location.txt")

# - Covariates file name
covariates_file_name = paste0(DIR,"input/","Covariates.",POP,".txt")

# - output file name
output_file_name_cis = paste0(DIR,"output/",POP,".cis_aQTL_all_control_gene_exprs.txt")
output_file_name_tra = paste0(DIR,"output/",POP,".trans_aQTL_all_control_gene_exprs.txt")
output_figure_name_cis = paste0(DIR,"output/",POP,".cis_aQTL_genotype_info_control_gene_exprs.pdf")
pdf(output_figure_name_cis)

# - threshold
pvOutputThreshold_cis = 1e-2;
pvOutputThreshold_tra = 1e-5;

# - Error covariance matrix
# set to numeric() for identity
errorCovariance = numeric();

# - Distance for local gene-SNP pairs
cisDist = CIS_DISTANCE;

# -- load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";
snps$fileOmitCharacters = "NA";
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

# -- load apa expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

# -- load covariates data
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
	cvrt$LoadFile(covariates_file_name);
}


## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
		      snps = snps,
		      gene = gene,
		      cvrt = cvrt,
		      output_file_name     = output_file_name_tra,
		      pvOutputThreshold     = pvOutputThreshold_tra,
		      useModel = useModel,
		      errorCovariance = errorCovariance,
		      verbose = TRUE,
		      output_file_name.cis = output_file_name_cis,
		      pvOutputThreshold.cis = pvOutputThreshold_cis,
		      snpspos = snpspos,
		      genepos = genepos,
		      cisDist = cisDist,
		      pvalue.hist = "qqplot",
		      min.pv.by.genesnp = TRUE,
		      noFDRsaveMemory = FALSE);

gz1 <- paste0(paste0(DIR,"output/",POP,".data.RData"),"w")
save.image(gz1)
#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

# -- Results
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
cat('Detected local aQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant aQTLs:', '\n');
show(me$trans$eqtls)

write.table(me$cis$min.pv.gene, paste0(DIR,"output/",POP, ".cis.min.pv.gene.txt"))
save(gene,snps,file=paste0(DIR,"output/",POP,".RData"))
save(snps,genepos,file=paste0(DIR,"output/",POP,".permutation.RData"))

# -- plot the Q-Q plot of local and distant p-values
plot(me)
dev.off()
