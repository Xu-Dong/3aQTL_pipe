#!/opt/app/languages/R-3.6.3/bin/Rscript

# prepare input files for matrix-eqtl

# -- load libraries
library(dplyr)
library(magrittr)
library(data.table)

# -- Main
setwd("/home/username/Project_XXX/output")
# --- load input data
# imputed PUDI matrix, sex info, PCA of genotype
pdui_mat.ceu <- readRDS("CEU.pdui_mat.imputed.RDS")
pdui_mat.fin <- readRDS("FIN.pdui_mat.imputed.RDS")
pdui_mat.gbr <- readRDS("GBR.pdui_mat.imputed.RDS")
pdui_mat.tsi <- readRDS("TSI.pdui_mat.imputed.RDS")
pdui_mat.yri <- readRDS("YRI.pdui_mat.imputed.RDS")

pdui_mat.ceu <- as.data.frame(pdui_mat.ceu);pdui_mat.ceu$id <- rownames(pdui_mat.ceu)
pdui_mat.fin <- as.data.frame(pdui_mat.fin);pdui_mat.fin$id <- rownames(pdui_mat.fin)
pdui_mat.gbr <- as.data.frame(pdui_mat.gbr);pdui_mat.gbr$id <- rownames(pdui_mat.gbr)
pdui_mat.tsi <- as.data.frame(pdui_mat.tsi);pdui_mat.tsi$id <- rownames(pdui_mat.tsi)
pdui_mat.yri <- as.data.frame(pdui_mat.yri);pdui_mat.yri$id <- rownames(pdui_mat.yri)

# current format is:
# subject_id	population(CEU/FIN/GBR/TSI/YRI)	sex(female/male)
sub2pops <- read.table("../input/subject_population_sex.445.txt",header=T, sep="\t",stringsAsFactors=F)
pop_sex.ceu <- sub2pops %>% dplyr::filter(population == "CEU") %>% select(subject_id)
pop_sex.fin <- sub2pops %>% dplyr::filter(population == "FIN") %>% select(subject_id)
pop_sex.gbr <- sub2pops %>% dplyr::filter(population == "GBR") %>% select(subject_id)
pop_sex.tsi <- sub2pops %>% dplyr::filter(population == "TSI") %>% select(subject_id)
pop_sex.yri <- sub2pops %>% dplyr::filter(population == "YRI") %>% select(subject_id)

# --- reformating PDUI matrix for peer
pdui_mat.ceu %<>% dplyr::select("id",all_of(pop_sex.ceu$subject_id))
pdui_mat.fin %<>% dplyr::select("id",all_of(pop_sex.fin$subject_id))
pdui_mat.gbr %<>% dplyr::select("id",all_of(pop_sex.gbr$subject_id))
pdui_mat.tsi %<>% dplyr::select("id",all_of(pop_sex.tsi$subject_id))
pdui_mat.yri %<>% dplyr::select("id",all_of(pop_sex.yri$subject_id))

write.table(pdui_mat.ceu,file="./matrix-eqtl/input/PDUI_mat.CEU.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(pdui_mat.fin,file="./matrix-eqtl/input/PDUI_mat.FIN.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(pdui_mat.gbr,file="./matrix-eqtl/input/PDUI_mat.GBR.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(pdui_mat.tsi,file="./matrix-eqtl/input/PDUI_mat.TSI.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(pdui_mat.yri,file="./matrix-eqtl/input/PDUI_mat.YRI.txt",quote=F,sep="\t",row.names=F,col.names=T)

# --- load covariates
covariates.ceu <- read.table("ceu.pdui.peer.covariates.txt",header=T,sep="\t",stringsAsFactors=F)
covariates.fin <- read.table("fin.pdui.peer.covariates.txt",header=T,sep="\t",stringsAsFactors=F)
covariates.gbr <- read.table("gbr.pdui.peer.covariates.txt",header=T,sep="\t",stringsAsFactors=F)
covariates.tsi <- read.table("tsi.pdui.peer.covariates.txt",header=T,sep="\t",stringsAsFactors=F)
covariates.yri <- read.table("yri.pdui.peer.covariates.txt",header=T,sep="\t",stringsAsFactors=F)

covariates.ceu %<>% dplyr::select("id",all_of(pop_sex.ceu$subject_id))
covariates.fin %<>% dplyr::select("id",all_of(pop_sex.fin$subject_id))
covariates.gbr %<>% dplyr::select("id",all_of(pop_sex.gbr$subject_id))
covariates.tsi %<>% dplyr::select("id",all_of(pop_sex.tsi$subject_id))
covariates.yri %<>% dplyr::select("id",all_of(pop_sex.yri$subject_id))

write.table(covariates.ceu,file="./matrix-eqtl/input2/Covariates.CEU.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(covariates.fin,file="./matrix-eqtl/input2/Covariates.FIN.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(covariates.gbr,file="./matrix-eqtl/input2/Covariates.GBR.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(covariates.tsi,file="./matrix-eqtl/input2/Covariates.TSI.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(covariates.yri,file="./matrix-eqtl/input2/Covariates.YRI.txt",quote=F,sep="\t",row.names=F,col.names=T)


# --- load genotype matrix
gt_mat <- fread("../input/E-GEUV-1-Genotype/GEUVADIS.all_chrs.gt012.bed",header=T)
gt_mat <- as.data.frame(gt_mat)

gt_mat.ceu <- gt_mat %>% dplyr::select("id",all_of(pop_sex.ceu$subject_id))
gt_mat.fin <- gt_mat %>% dplyr::select("id",all_of(pop_sex.fin$subject_id))
gt_mat.gbr <- gt_mat %>% dplyr::select("id",all_of(pop_sex.gbr$subject_id))
gt_mat.tsi <- gt_mat %>% dplyr::select("id",all_of(pop_sex.tsi$subject_id))
gt_mat.yri <- gt_mat %>% dplyr::select("id",all_of(pop_sex.yri$subject_id))

write.table(gt_mat.ceu,file="./matrix-eqtl/input2/Genotype_mat.CEU.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(gt_mat.fin,file="./matrix-eqtl/input2/Genotype_mat.FIN.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(gt_mat.gbr,file="./matrix-eqtl/input2/Genotype_mat.GBR.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(gt_mat.tsi,file="./matrix-eqtl/input2/Genotype_mat.TSI.txt",quote=F,sep="\t",row.names=F,col.names=T)
write.table(gt_mat.yri,file="./matrix-eqtl/input2/Genotype_mat.YRI.txt",quote=F,sep="\t",row.names=F,col.names=T)

print("Done!")
