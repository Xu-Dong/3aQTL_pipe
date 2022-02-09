#!/opt/app/languages/R-4.1.0/bin/Rscript

# merge Dapars2 output by chromosome
# -- global variable
dir <- "/home/username/Project_XXX" # change this path to user defined root path of the project
# -- functions

load_dapars2_res <- function(chromosome,new_header){
	input_file <- paste0(dir,"output/Dapars2_geuv_res_",chromosome,"/Dapars2_GEUV_result_temp.",chromosome,".txt")
	dap_res <- read.table(input_file,header=T, sep="\t")
	names(dap_res) <- new_header

	return(dap_res)
}


# -- main

# load samples
dat <- read.table(paste0(dir,"input/subject_id2population.txt"),header=T,sep="\t")
sample_list <- dat$subject_id
col_names <- c("Gene","fit_value","Predicted_Proximal_APA","Loci",sample_list)
chrs <- paste0("chr",c(seq(22),"X","Y"))

res.df <- data.frame()

for(chr in chrs){
	temp.df <- load_dapars2_res(chr,col_names)
	print(paste(chr,dim(temp.df)[1],sep=":"))
	res.df <- rbind(res.df,temp.df)
}

dim(res.df)
write.table(res.df,file="Dapars2_geuv_res.all_chromosomes.txt",quote=F,sep="\t",row.names=F)


