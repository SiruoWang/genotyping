library(dplyr)
library(ranger)
library(rtracklayer)
library(matrixStats)
library(parallel)

load("./rda/chrlist.rda")
load("./rda/poslist.rda")
load("./rda/score_A.rda")
load("./rda/score_T.rda")
load("./rda/score_C.rda")
load("./rda/score_G.rda")
load("./rda/known_genotype_blood.rda")
load("./rda/snp_major_minor.rda")
load("./rda/sra_gtex_blood.rda")
load("./rda/subjid_train.rda")

## build genotype call training model on all geuvadis samples
genotype_call_model <- function(use_seed=12345){

  set.seed(use_seed)

  result <- list()

  subjects <- colnames(known_genotype_blood)
  snp_subset_A <-  snp_subset_C <- snp_subset_G <- snp_subset_T <- matrix(rep(0,length(nrow(score_bw)) * length(subjects)), nrow = length(nrow(score_bw)))
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
  rm(score_A)
  rm(score_C)
  rm(score_G)
  rm(score_T)
  gc()

  train_idx <- sample(length(subjects), 2/3 * length(subjects))

  genotype_train <- c(known_genotype_blood[,train_idx])
  major_train <- rep(snp_major_minor[,1], length(subjects[train_idx]))
  minor_train <- rep(snp_major_minor[,2], length(subjects[train_idx]))
  nucA_train <- c(snp_subset_A[,train_idx])
  nucC_train <- c(snp_subset_C[,train_idx])
  nucG_train <- c(snp_subset_G[,train_idx])
  nucT_train <- c(snp_subset_T[,train_idx])

  genotype_test <- c(known_genotype_blood[,-train_idx])
  major_test <- rep(snp_major_minor[,1], length(subjects[-train_idx]))
  minor_test <- rep(snp_major_minor[,2], length(subjects[-train_idx]))
  nucA_test <- c(snp_subset_A[,-train_idx])
  nucC_test <- c(snp_subset_C[,-train_idx])
  nucG_test <- c(snp_subset_G[,-train_idx])
  nucT_test <- c(snp_subset_T[,-train_idx])

  train_df <- data.frame(genotype = as.factor(genotype_train),
                         nucA = nucA_train, nucC = nucC_train, nucG = nucG_train, nucT = nucT_train,
                         major = major_train, minor = minor_train)
  model_allpos <- ranger(genotype ~ nucA + nucC + nucG + nucT + major + minor, data = train_df, num.trees = 100, write.forest = TRUE)
  result[[1]] <- model_allpos

  
  pred.snp <- predict(rg.snp, dat = snp.test)



  return(result)
}


subjects_blood <- colnames(known_genotype_blood)



## geuvadis_training_model.rda: a list of length 50531, each element of which is also a list with the first element recording snps indexes that used to train the snp specific model, the second element recording the snp specific model
training_model <- mclapply(colnames(known_genotype_blood),
                                    genotype_call_model,
                                    mc.cores = 10)
save(training_model, file = "./rda/training_model_ratio.rda", compress = TRUE)

q(save = "no")
# shell script: echo R CMD BATCH genotype_model_call_train_ratio.R | qsub -cwd -pe local 10 -R y -l mf=20G,h_vmem=20G
