library(dplyr)
library(ranger)
library(rtracklayer)
library(matrixStats)
library(parallel)

load("./score_matrix/score/filter.rda")
load("./score_matrix/score/score_A.rda")
load("./score_matrix/score/score_T.rda")
load("./score_matrix/score/score_C.rda")
load("./score_matrix/score/score_G.rda")
#load("./score_matrix/score/ref_nucs.rda")

load("./score_matrix/score/known_genotype_blood.rda")
load("./rda/subjid_train.rda")



###################################################################
## build genotype call training model on all geuvadis samples
genotype_call_model <- function(idx, range_idx){
  result <- list()
  ## build GRange objects for snp with snp index idx = sample_snps_idx[i]: ranges of each sample snp centers at the snp position, and with range diameter 2000
  single_snp_gr <- GRanges(seqnames=chrlist[idx],
                           ranges=IRanges(poslist[idx] - 250 * range_idx, poslist[idx] + 250 * range_idx),
                           strand="*")

  distance_cluster_snps_idx <- findOverlaps(single_snp_gr,filter) %>% subjectHits()

  result[[1]] <- distance_cluster_snps_idx


  subjects <- colnames(known_genotype_blood)
  ini_mat <- matrix(rep(0,length(distance_cluster_snps_idx) * length(subjects)), nrow = length(distance_cluster_snps_idx))
  ini_mat <- apply(ini_mat,c(1,2),as.integer
  snp_subset_A <-  snp_subset_C <- snp_subset_G <- snp_subset_T <- ini_mat
  rm(ini_mat)
  gc()

  colnames(snp_subset_A) <- colnames(snp_subset_C) <- colnames(snp_subset_G) <- colnames(snp_subset_T) <- subjects

  for (i in 1:length(subjects)){
    runs <- sra_gtex_blood[which(sra_gtex_blood$SUBJID == subjects[i]), "Run"]
    col_id <- match(runs, colnames(score_bw))
    snp_subset_A[,i] <- rowSums(score_A[,col_id,drop = FALSE])
    snp_subset_C[,i] <- rowSums(score_C[,col_id,drop = FALSE])
    snp_subset_G[,i] <- rowSums(score_G[,col_id,drop = FALSE])
    snp_subset_T[,i] <- rowSums(score_T[,col_id,drop = FALSE])
  }

  snp_total_counts <- snp_subset_A + snp_subset_C + snp_subset_G + snp_subset_T
  ## use all snps' ACGT in each distence cluster to predict one snp's genotype
  ## t(snp_subset_A): each row is one snp's A counts. number of rows equal number of snps in this distance cluster
  if (length(distance_cluster_snps_idx)==1) {
    snpcombine <- data.frame(genotype = as.factor(known_genotype_blood[idx,]),
                             nucA = snp_subset_A / snp_total_counts,# A counts for all snps in this distance cluster
                             nucC = snp_subset_C / snp_total_counts,# C counts for all snps in this distance cluster
                             nucG = snp_subset_G / snp_total_counts,# G counts for all snps in this distance cluster
                             nucT = snp_subset_T / snp_total_counts)
    snpcombine[is.na(snpcombine)] <- 0

  }else{
    snpcombine <- data.frame(genotype = as.factor(known_genotype_blood[idx,]),
                             nucA = t(snp_subset_A / snp_total_counts),# A counts for all snps in this distance cluster
                             nucC = t(snp_subset_C / snp_total_counts),# C counts for all snps in this distance cluster
                             nucG = t(snp_subset_G / snp_total_counts),# G counts for all snps in this distance cluster
                             nucT = t(snp_subset_T / snp_total_counts))# T counts for all snps in this distance cluster
    snpcombine[is.na(snpcombine)] <- 0
  }

  model_snpcombine <- ranger(genotype ~ ., data = snpcombine, num.trees = 100, write.forest = TRUE)

  result[[2]] <- model_snpcombine

  return(result)
}




chrlist = data.frame(chr = seqnames(filter))
poslist = pos(GPos(filter))
training_model <- mclapply(c(1:length(poslist)),
                                    genotype_call_model,
                                    range_idx = 2,
                                    mc.cores = 10)
save(training_model, file = "./rda/training_model_ratio.rda", compress = TRUE)

q(save = "no")
# shell script: echo R CMD BATCH genotype_model_call_train_ratio.R | qsub -cwd -pe local 10 -R y -l mf=20G,h_vmem=20G
