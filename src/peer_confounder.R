#!/opt/app/languages/R-3.6.3/bin/Rscript

# perform covariates analysis with peer

# -- load libraries
library(peer)
library(dplyr)
library(magrittr)

# -- Main
setwd("/home/username/Project_XXX/output")
# --- load input data
# imputed PUDI matrix, sex info, PCA of genotype
pdui_mat <- read.table("Dapars2_geuv_PDUIs.imputed.txt",header=T,sep="\t",stringsAsFactors=F)

# genotyep pca
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

# subject population info
# current format is:
# subject_id	population(CEU/FIN/GBR/TSI/YRI)	sex(female/male)
sub2pops <- read.table("../input/subject_population_sex.445.txt",header=T, sep="\t",stringsAsFactors=F)
sub2pops$sex[sub2pops$sex=="male"] <- 1
sub2pops$sex[sub2pops$sex=="female"] <- 0
pop_sex.ceu <- sub2pops %>% dplyr::filter(population == "CEU") %>% select(subject_id,sex)
pop_sex.fin <- sub2pops %>% dplyr::filter(population == "FIN") %>% select(subject_id,sex)
pop_sex.gbr <- sub2pops %>% dplyr::filter(population == "GBR") %>% select(subject_id,sex)
pop_sex.tsi <- sub2pops %>% dplyr::filter(population == "TSI") %>% select(subject_id,sex)
pop_sex.yri <- sub2pops %>% dplyr::filter(population == "YRI") %>% select(subject_id,sex)

# merge sex with top 5 pca
pop_sex.ceu <- merge(pop_sex.ceu,gt_pca.ceu,by="subject_id")
pop_sex.fin <- merge(pop_sex.fin,gt_pca.fin,by="subject_id")
pop_sex.gbr <- merge(pop_sex.gbr,gt_pca.gbr,by="subject_id")
pop_sex.tsi <- merge(pop_sex.tsi,gt_pca.tsi,by="subject_id")
pop_sex.yri <- merge(pop_sex.yri,gt_pca.yri,by="subject_id")

pop_sex.ceu$sex <- as.numeric(pop_sex.ceu$sex)
pop_sex.fin$sex <- as.numeric(pop_sex.fin$sex)
pop_sex.gbr$sex <- as.numeric(pop_sex.gbr$sex)
pop_sex.tsi$sex <- as.numeric(pop_sex.tsi$sex)
pop_sex.yri$sex <- as.numeric(pop_sex.yri$sex)

# --- reformating PDUI matrix for peer
pdui_mat.ceu <- pdui_mat %>% dplyr::select(all_of(pop_sex.ceu$subject_id))
pdui_mat.fin <- pdui_mat %>% dplyr::select(all_of(pop_sex.fin$subject_id))
pdui_mat.gbr <- pdui_mat %>% dplyr::select(all_of(pop_sex.gbr$subject_id))
pdui_mat.tsi <- pdui_mat %>% dplyr::select(all_of(pop_sex.tsi$subject_id))
pdui_mat.yri <- pdui_mat %>% dplyr::select(all_of(pop_sex.yri$subject_id))


rownames(pop_sex.ceu) <- pop_sex.ceu$subject_id;pop_sex.ceu <- as.matrix(pop_sex.ceu[,2:7])
rownames(pop_sex.fin) <- pop_sex.fin$subject_id;pop_sex.fin <- as.matrix(pop_sex.fin[,2:7])
rownames(pop_sex.gbr) <- pop_sex.gbr$subject_id;pop_sex.gbr <- as.matrix(pop_sex.gbr[,2:7])
rownames(pop_sex.tsi) <- pop_sex.tsi$subject_id;pop_sex.tsi <- as.matrix(pop_sex.tsi[,2:7])
rownames(pop_sex.yri) <- pop_sex.yri$subject_id;pop_sex.yri <- as.matrix(pop_sex.yri[,2:7])


# add to list
covariate.list <- list(CEU=t(pop_sex.ceu),FIN=t(pop_sex.fin),GBR=t(pop_sex.gbr),TSI=t(pop_sex.tsi), YRI=t(pop_sex.yri))

# transpose PDUI matrix
pdui_mat.ceu <- t(pdui_mat.ceu)
pdui_mat.fin <- t(pdui_mat.fin)
pdui_mat.gbr <- t(pdui_mat.gbr)
pdui_mat.tsi <- t(pdui_mat.tsi)
pdui_mat.yri <- t(pdui_mat.yri)

# add to a running list
pdui_mat.list <- list(CEU=pdui_mat.ceu,FIN=pdui_mat.fin, GBR=pdui_mat.gbr, TSI=pdui_mat.tsi, YRI=pdui_mat.yri)

# --- run peer
print("Start peer...")
loop.pop <- c("ceu","fin","gbr","tsi","yri")
for(i in 1:5){
	print(paste("Start peer",loop.pop[i]))
	# 15 factors for sample number < 150
	model = PEER() # 1.create the model object

	PEER_setPhenoMean(model,pdui_mat.list[[i]]) # 2.set the observed data
	PEER_setNk(model,15)

	#PEER_setCovariates(model,covariate.list[[i]]) # 3. add known covariates

	PEER_update(model) # 4.perform the inference
	# --- extract output
	# the result is the model object
	# we can get the posterior mean of the inferred confounders (NxK matrix)
	# their weights (GxK matrix), precision (inverse variance) of the weights, and the residual dataset (NxG)

	factors <- PEER_getX(model)
	rownames(factors) <- rownames(pdui_mat.list[[i]])
	colnames(factors) <- paste0("PEER_",1:15)
	factors.t <- t(factors)
	known_factors <- covariate.list[[i]]
	known_factors[,match(colnames(factors.t),colnames(known_factors))]
	factors.combined <- rbind(factors.t,known_factors)
	saveRDS(factors.t,file=paste0("factors.",loop.pop[i],".RDS"))
	write.table(factors.combined,file=paste0("Covariates.",loop.pop[i],".txt"),quote=F,sep="\t",row.names=T,col.names=T)

	residuals <- PEER_getResiduals(model)
	rownames(residuals) <- rownames(pdui_mat.list[[i]])
	colnames(residuals) <- colnames(pdui_mat.list[[i]])
	residuals.t <- t(residuals)
	saveRDS(residuals.t,file=paste0("residuals_pduis.",loop.pop[i],".RDS"))
	write.table(residuals.t,file=paste0("residuals_pduis.",loop.pop[i],".txt"),quote=F,sep="\t",row.names=T,col.names=T)
	rm(model,factors,factors.t,residuals,residuals.t,known_factors,factors.combined)
}

print("Done")
