#!/opt/app/languages/R-3.6.3/bin/Rscript

args <- commandArgs(trailingOnly=TRUE)
POP <- args[1] #CEU/FIN/GBR/TSI/YRI
basedir <- "/media/Rome/zouxd/Projects/2021-12-05-3aQTL-STAR-Protocol-Project/output/susieR_analysis"

setwd(basedir)

gene_list <- read.table(paste0("./input/",POP,".aGenes.loc_1Mb.txt"),header=F,sep="\t")
gene_list <- gene_list[,1]
independent_snp_count <- c()

susie_df <- data.frame(locus_id=c(),variant_id=c(),pip=c(),cs=c(),cs_size=c(),cs_purity=c())
for(idx in 1:length(gene_list)){
	file_name <- paste0("./",POP,"/",gene_list[idx],"/3aQTL.SuSiE.txt")
	cat(file_name)
	if (file.exists(file_name)){
		df <- read.table(file_name,header=T,sep=" ")
		independent_snp_count[idx] <- dim(df)[1]

		if (dim(df)[1]>0){
			df$locus_id <- gene_list[idx]
			susie_df <- rbind(susie_df,df)
		}
	}else{
		independent_snp_count[idx] <- NA
	}
}

summary_susie <- data.frame(Gene=gene_list,Count=independent_snp_count)

write.table(susie_df,file=paste0("susieR_res.all_genes.",POP,".txt"),quote=F,row.names=F,sep="\t")
write.table(summary_susie,file=paste0("susieR_res.stat.",POP,".txt"),quote=F,row.names=F,sep="\t")
