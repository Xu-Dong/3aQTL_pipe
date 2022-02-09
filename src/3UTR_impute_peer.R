#!/opt/app/languages/R-3.6.3/bin/Rscript
# 2022-01-08

# print usage
usage <- function() {
  cat(
    'usage: Rscript 3UTR_impute_peer.R
')
}

library(dplyr)
library(magrittr)
library(peer)
library(impute)

setwd("/home/username/Project_XXX/output") # change this path to user defined project workspace

# --------------- load known covariates -----------------## --------------------------------
# genotype pca
gt_pca.ceu <- read.table("./PCA_genotype/CEU.eigenvec",header=F,sep=" ",stringsAsFactors=F)
gt_pca.fin <- read.table("./PCA_genotype/FIN.eigenvec",header=F,sep=" ",stringsAsFactors=F)
gt_pca.gbr <- read.table("./PCA_genotype/GBR.eigenvec",header=F,sep=" ",stringsAsFactors=F)
gt_pca.tsi <- read.table("./PCA_genotype/TSI.eigenvec",header=F,sep=" ",stringsAsFactors=F)
gt_pca.yri <- read.table("./PCA_genotype/YRI.eigenvec",header=F,sep=" ",stringsAsFactors=F)

gt_pca.ceu$V1 <- NULL;gt_pca.ceu <- gt_pca.ceu[,1:6]
gt_pca.fin$V1 <- NULL;gt_pca.fin <- gt_pca.fin[,1:6]
gt_pca.gbr$V1 <- NULL;gt_pca.gbr <- gt_pca.gbr[,1:6]
gt_pca.tsi$V1 <- NULL;gt_pca.tsi <- gt_pca.tsi[,1:6]
gt_pca.yri$V1 <- NULL;gt_pca.yri <- gt_pca.yri[,1:6]

names(gt_pca.ceu) <- c("subject_id","PC_1","PC_2","PC_3","PC_4","PC_5")
names(gt_pca.fin) <- c("subject_id","PC_1","PC_2","PC_3","PC_4","PC_5")
names(gt_pca.gbr) <- c("subject_id","PC_1","PC_2","PC_3","PC_4","PC_5")
names(gt_pca.tsi) <- c("subject_id","PC_1","PC_2","PC_3","PC_4","PC_5")
names(gt_pca.yri) <- c("subject_id","PC_1","PC_2","PC_3","PC_4","PC_5")

# subject population info, can be download in data/ of the current repository
# current format is:
# subject_id    population(CEU/FIN/GBR/TSI/YRI) sex(female/male)
sub2pops <- read.table("/home/username/Project_XXX/input/subject_population_sex.445.txt",header=T, sep="\t",stringsAsFactors=F)
sub2pops$sex[sub2pops$sex=="male"] <- 1
sub2pops$sex[sub2pops$sex=="female"] <- 0
known_cov.ceu <- sub2pops %>% dplyr::filter(population == "CEU") %>% select(subject_id,sex)
known_cov.fin <- sub2pops %>% dplyr::filter(population == "FIN") %>% select(subject_id,sex)
known_cov.gbr <- sub2pops %>% dplyr::filter(population == "GBR") %>% select(subject_id,sex)
known_cov.tsi <- sub2pops %>% dplyr::filter(population == "TSI") %>% select(subject_id,sex)
known_cov.yri <- sub2pops %>% dplyr::filter(population == "YRI") %>% select(subject_id,sex)

# merge sex with top 5 pca
known_cov.ceu <- merge(known_cov.ceu,gt_pca.ceu,by="subject_id")
known_cov.fin <- merge(known_cov.fin,gt_pca.fin,by="subject_id")
known_cov.gbr <- merge(known_cov.gbr,gt_pca.gbr,by="subject_id")
known_cov.tsi <- merge(known_cov.tsi,gt_pca.tsi,by="subject_id")
known_cov.yri <- merge(known_cov.yri,gt_pca.yri,by="subject_id")

known_cov.ceu$sex <- as.numeric(known_cov.ceu$sex)
known_cov.fin$sex <- as.numeric(known_cov.fin$sex)
known_cov.gbr$sex <- as.numeric(known_cov.gbr$sex)
known_cov.tsi$sex <- as.numeric(known_cov.tsi$sex)
known_cov.yri$sex <- as.numeric(known_cov.yri$sex)

rm(gt_pca.ceu,gt_pca.fin,gt_pca.gbr,gt_pca.tsi,gt_pca.yri)
# --------------- load PDUI matrix ------ # ------------------
pdui_mat <- read.table("Dapars2_geuv_res.all_chromosomes.txt", stringsAsFactors=FALSE, header=TRUE,sep="\t")
pdui_mat <- pdui_mat[,-c(2,3,4)]

pdui_mat.ceu <- pdui_mat %>% dplyr::select(all_of(known_cov.ceu$subject_id))
pdui_mat.fin <- pdui_mat %>% dplyr::select(all_of(known_cov.fin$subject_id))
pdui_mat.gbr <- pdui_mat %>% dplyr::select(all_of(known_cov.gbr$subject_id))
pdui_mat.tsi <- pdui_mat %>% dplyr::select(all_of(known_cov.tsi$subject_id))
pdui_mat.yri <- pdui_mat %>% dplyr::select(all_of(known_cov.yri$subject_id))

pdui_mat.ceu <- as.matrix(pdui_mat.ceu)
pdui_mat.fin <- as.matrix(pdui_mat.fin)
pdui_mat.gbr <- as.matrix(pdui_mat.gbr)
pdui_mat.tsi <- as.matrix(pdui_mat.tsi)
pdui_mat.yri <- as.matrix(pdui_mat.yri)
rownames(pdui_mat.ceu) <- pdui_mat[,1]
rownames(pdui_mat.fin) <- pdui_mat[,1]
rownames(pdui_mat.gbr) <- pdui_mat[,1]
rownames(pdui_mat.tsi) <- pdui_mat[,1]
rownames(pdui_mat.yri) <- pdui_mat[,1]
#remove genes with more than 50% entries missing and individuals with more than 80% missing data
pdui_mat.ceu <- pdui_mat.ceu[, colMeans(is.na(pdui_mat.ceu)) <= 0.8];pdui_mat.ceu <-  pdui_mat.ceu[rowMeans(is.na(pdui_mat.ceu)) < 0.5,]
pdui_mat.fin <- pdui_mat.fin[, colMeans(is.na(pdui_mat.fin)) <= 0.8];pdui_mat.fin <-  pdui_mat.fin[rowMeans(is.na(pdui_mat.fin)) < 0.5,]
pdui_mat.gbr <- pdui_mat.gbr[, colMeans(is.na(pdui_mat.gbr)) <= 0.8];pdui_mat.gbr <-  pdui_mat.gbr[rowMeans(is.na(pdui_mat.gbr)) < 0.5,]
pdui_mat.tsi <- pdui_mat.tsi[, colMeans(is.na(pdui_mat.tsi)) <= 0.8];pdui_mat.tsi <-  pdui_mat.tsi[rowMeans(is.na(pdui_mat.tsi)) < 0.5,]
pdui_mat.yri <- pdui_mat.yri[, colMeans(is.na(pdui_mat.yri)) <= 0.8];pdui_mat.yri <-  pdui_mat.yri[rowMeans(is.na(pdui_mat.yri)) < 0.5,]


class(pdui_mat.ceu) <- 'numeric'
class(pdui_mat.fin) <- 'numeric'
class(pdui_mat.gbr) <- 'numeric'
class(pdui_mat.tsi) <- 'numeric'
class(pdui_mat.yri) <- 'numeric'

# convert data.frame to matrix 
rownames(known_cov.ceu) <- known_cov.ceu$subject_id;known_cov.ceu <- as.matrix(known_cov.ceu[,2:7])
rownames(known_cov.fin) <- known_cov.fin$subject_id;known_cov.fin <- as.matrix(known_cov.fin[,2:7])
rownames(known_cov.gbr) <- known_cov.gbr$subject_id;known_cov.gbr <- as.matrix(known_cov.gbr[,2:7])
rownames(known_cov.tsi) <- known_cov.tsi$subject_id;known_cov.tsi <- as.matrix(known_cov.tsi[,2:7])
rownames(known_cov.yri) <- known_cov.yri$subject_id;known_cov.yri <- as.matrix(known_cov.yri[,2:7])

covariate.list <- list(CEU=known_cov.ceu,FIN=known_cov.fin,GBR=known_cov.gbr,TSI=known_cov.tsi, YRI=known_cov.yri)

pdui_mat.list <- list(CEU=pdui_mat.ceu,FIN=pdui_mat.fin, GBR=pdui_mat.gbr, TSI=pdui_mat.tsi, YRI=pdui_mat.yri)
# ----------------- run peer -------------------- # --------------------------------------------------

save.image(file="run_peer_impute.RData")
loop.pop <- c("ceu","fin","gbr","tsi","yri")
for(i in 1:5){
	print(loop.pop[i])
	model <- PEER()
	covs_se <- covariate.list[[i]]

	PEER_setCovariates(model, covs_se)
	dim(PEER_getCovariates(model))
	# impute missing values in 3'UTR expression
	mat.ds <- pdui_mat.list[[i]]
	mat_impute <- impute.knn(mat.ds)
	# save imputed PDUI matrix without quantile normalization
	saveRDS(as.matrix(mat_impute$data),file=paste0(loop.pop[i],".pdui_mat.imputed.RDS"))
	
	# save imputed PDUI matrix after quantile normalization
	df_w <- as.data.frame(mat_impute$data)
	for(gene iin 1:nrow(df_w)){
		mat = df_w[gene,]
		mat = apply(mat,1,rank,ties.method = "average")
		mat = qnorm(mat / (ncol(df_w)+1))
		df_w[gene,] = mat
	}
	saveRDS(cbind(rownames(mat_impute$data),df_w),file=paste0(loop.pop[i],".pdui_mat.imputed_qnorm.RDS"))
	write.table(cbind(rownames(mat_impute$data),df_w),file=paste0(loop.pop[i],".pdui_mat.imputed_qnorm.txt",row.names=F,col.names=T,quote=F,sep="\t"))
	
	PEER_setPhenoMean(model, t(as.matrix(mat_impute$data)))

	dim(PEER_getPhenoMean(model))


	# set number of peer factors
	## N < 150, use 15  PEERs, 150<=N<250, use 30 PEERs, N >=250 use 35 PEERs
	if (ncol(mat.ds) < 150) {
		numcov <- 15
	} else if (ncol(mat.ds) < 250) {
		numcov <- 30
	} else if (ncol(mat.ds) >= 250) {
		numcov <- 35
	}

	PEER_setNk(model, numcov)
	PEER_getNk(model)

	PEER_update(model)

	# diag
	pdf(paste0(loop.pop[i], '.peer.diag.pdf'), width=6, height=8)
	PEER_plotModel(model)
	dev.off()


	factors = t(PEER_getX(model))
	weights = PEER_getW(model)
	precision = PEER_getAlpha(model)

	residuals = t(PEER_getResiduals(model))
	rownames(residuals) <- rownames(mat.ds)
	colnames(residuals) <- colnames(mat.ds)

	rownames(factors) <- c("sex", "PC1", "PC2", "PC3", "PC4", "PC5", paste0("PEER_",1:numcov))
	colnames(factors) <- colnames(mat.ds)

	residuals.ds <- residuals

	#png(paste0(loop.pop[i], '.expr.peer.clust.png'), width=8, height=8, res=150, units='in')
	#heatmap.2(as.matrix(residuals.ds), distfun=function(x) dist(x,method='euclidian'), hclustfun=function(x) hclust(x,method='ward.D2'),
        #  trace='none', dendrogram='both', Rowv=TRUE, Colv=TRUE, breaks=pairs.breaks, col=colorRampPalette(myCols), scale='none', symkey=T, na.color='grey', density.info='histogram', cexRow=0.2, cexCol=0.5, main=paste0(TISSUE, '\nexpr clustering'))
	#dev.off()

	covariate_file <- paste0(loop.pop[i], ".pdui.peer.covariates.txt")
	write.table(cbind(rownames(factors), factors), file=covariate_file, row.names=FALSE, col.names=c("id",colnames(factors)), quote=FALSE, sep='\t')

	gz1 <- paste0(loop.pop[i], ".pdui.peer.residuals.txt")
	write.table(cbind(rownames(residuals), residuals), file=gz1, row.names=FALSE, col.names=c("id",colnames(residuals)), quote=FALSE, sep='\t')
	rm(model,cov_se,mat.ds,mat_impute,factors,weights,precision,residuals,covariate_file,gz1)
}
