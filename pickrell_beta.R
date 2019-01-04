library(rtracklayer)
library(ranger)
library(caret)
library(parallel)

pickrell_data <- read.table("./rda/pickrell_SraRunTable.txt", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
pickrell_data <- pickrell_data[,c("Run","cell_line")]
pickrell_data_uni <- pickrell_data[!duplicated(pickrell_data[,"cell_line"]),]

## only subset pickrell samples with known genotypes
load("./rda/conv_geno.rda")
overlap_cell_lines <- pickrell_data_uni[which(pickrell_data_uni[,"cell_line"] %in% colnames(conv_geno)), "cell_line"]
pickrell_subset <- pickrell_data[which(pickrell_data[,"cell_line"] %in% overlap_cell_lines),]

pickrell_prediction_prepare <- function(pickrell_runs){

  sra_path_list <- list()
  sra_A <- read.table("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/sra_A.txt",stringsAsFactors = FALSE)
  sra_C <- read.table("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/sra_C.txt",stringsAsFactors = FALSE)
  sra_G <- read.table("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/sra_G.txt",stringsAsFactors = FALSE)
  sra_T <- read.table("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/sra_T.txt",stringsAsFactors = FALSE)

  output_A_bw <- sapply(sra_A[,1], function(x) strsplit(x,"/")[[1]][9])
  output_C_bw <- sapply(sra_C[,1], function(x) strsplit(x,"/")[[1]][9])
  output_G_bw <- sapply(sra_G[,1], function(x) strsplit(x,"/")[[1]][9])
  output_T_bw <- sapply(sra_T[,1], function(x) strsplit(x,"/")[[1]][9])

  output_A <- sapply(output_A_bw, function(x) sub(".A.bw","",x))
  output_C <- sapply(output_C_bw, function(x) sub(".C.bw","",x))
  output_G <- sapply(output_G_bw, function(x) sub(".G.bw","",x))
  output_T <- sapply(output_T_bw, function(x) sub(".T.bw","",x))

  sharenames_AC <- intersect(output_A, output_C)
  sharenames_ACT <- intersect(sharenames_AC, output_T)
  sharenames_ACTG <- intersect(sharenames_ACT,output_G)

  pickrell_runs <- pickrell_runs[pickrell_runs %in% sharenames_ACTG]

  sra_path_A <- sra_A[match(pickrell_runs,output_A),]
  sra_path_C <- sapply(sra_path_A, function(x) sub(".A.bw",".C.bw",x))
  sra_path_G <- sapply(sra_path_A, function(x) sub(".A.bw",".G.bw",x))
  sra_path_T <- sapply(sra_path_A, function(x) sub(".A.bw",".T.bw",x))
  sra_path_bw <- sapply(sra_path_A, function(x) sub(".A.bw",".bw",x))

  sra_path_list[[1]] <- sra_path_A
  sra_path_list[[2]] <- sra_path_C
  sra_path_list[[3]] <- sra_path_G
  sra_path_list[[4]] <- sra_path_T
  sra_path_list[[5]] <- sra_path_bw
  sra_path_list[[6]] <- pickrell_runs

  return(sra_path_list)
}

load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/snp_ref_nucs.rda")
load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/chrlist.rda")
load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/poslist.rda")
load("./rda/beta_model.rda")

ratio_est <- function(counts, total, alpha, beta){
  return((counts + alpha) / (total + alpha+ beta))
}

genotype_prediction <- function(x, path){

  snp_pos_gr <- GRanges(seqnames=chrlist,ranges=IRanges(poslist,poslist),strand="*")

  scoreA_list <- scoreC_list <- scoreG_list <- scoreT_list <- scorebw_list <- genotype_list <- c()

  tempA <- import(path[[1]][x],format="bigwig") # path[[1]][x]: nucleotide A counts bigwig file path
  overlap_loci <- findOverlaps(snp_pos_gr,tempA)
  scoreA_list[queryHits(overlap_loci)] <- tempA$score[subjectHits(overlap_loci)]


  tempC <- import(path[[2]][x],format="bigwig") # path[[2]][x]: nucleotide C counts bigwig file path
  overlap_loci <- findOverlaps(snp_pos_gr,tempC)
  scoreC_list[queryHits(overlap_loci)] <- tempC$score[subjectHits(overlap_loci)]


  tempG <- import(path[[3]][x],format="bigwig") # path[[3]][x]: nucleotide G counts bigwig file path
  overlap_loci <- findOverlaps(snp_pos_gr,tempG)
  scoreG_list[queryHits(overlap_loci)] <- tempG$score[subjectHits(overlap_loci)]


  tempT <- import(path[[4]][x],format="bigwig") # path[[4]][x]: nucleotide T counts bigwig file path
  overlap_loci <- findOverlaps(snp_pos_gr,tempT)
  scoreT_list[queryHits(overlap_loci)] <- tempT$score[subjectHits(overlap_loci)]


  tempbw <- import(path[[5]][x],format="bigwig") # path[[5]][x]: nucleotide bw counts bigwig file path
  overlap_loci <- findOverlaps(snp_pos_gr,tempbw)
  scorebw_list[queryHits(overlap_loci)] <- tempbw$score[subjectHits(overlap_loci)]


  rm(tempA)
  rm(tempC)
  rm(tempG)
  rm(tempT)
  rm(tempbw)
  gc()

  ref_score_list <- scorebw_list - scoreA_list - scoreC_list - scoreG_list - scoreT_list

  for (i in c(1:length(snp_ref_nucs))){
    if (snp_ref_nucs[i]=="A") {scoreA_list[i] <- ref_score_list[i]}
    else if (snp_ref_nucs[i]=="C") {scoreC_list[i] <- ref_score_list[i]}
    else if (snp_ref_nucs[i]=="G") {scoreG_list[i] <- ref_score_list[i]}
    else if (snp_ref_nucs[i]=="T") {scoreT_list[i] <- ref_score_list[i]}
  } ## end for (i in c(1:length(snp_ref_nucs))){


  #### make predictions using null_model
  total = scoreA_list + scoreC_list + scoreG_list + scoreT_list

  ## ratio A
  alpha.prior <- 0.04189064
  beta.prior <- 0.14943819
  A_est_valset <- ratio_est(scoreA_list, total, alpha.prior, beta.prior)
  gc()

  ## ratio C
  alpha.prior <- 0.05040109
  beta.prior <- 0.12654410
  C_est_valset <- ratio_est(scoreC_list, total, alpha.prior, beta.prior)
  gc()

  ## ratio G
  alpha.prior <- 0.04767691
  beta.prior <- 0.12520538
  G_est_valset <- ratio_est(scoreG_list, total, alpha.prior, beta.prior)
  gc()

  ## ratio T

  alpha.prior <- 0.04291391
  beta.prior <- 0.15213301
  T_est_valset <- ratio_est(scoreT_list, total, alpha.prior, beta.prior)
  gc()

  ratioA = A_est_valset
  ratioC = C_est_valset
  ratioG = G_est_valset
  ratioT = T_est_valset

  non_na_id = which(!is.na(ratioA))
  pred_geno = rep(NA, length(ratioA))

  rf_pred_df <- data.frame(ratioA = ratioA[non_na_id], ratioC = ratioC[non_na_id], ratioG = ratioG[non_na_id], ratioT = ratioT[non_na_id])
  prediction <- predict(beta_model, dat = rf_pred_df, num.threads = 20)

  pred_geno[non_na_id] <- as.character(prediction$predictions)

  return(pred_geno)

}

pickrell_path_list <- pickrell_prediction_prepare(pickrell_subset[,"Run"])
pred_geno <- mclapply(c(1:length(pickrell_path_list[[1]])),
                                    genotype_prediction,
                                    path = pickrell_path_list,
                                    mc.cores = 10)

save(pred_geno, file = "./rda/pickrell_beta_pred_geno.rda", compress = TRUE)

q(save = "no")

## shell script: echo R CMD BATCH pickrell_beta.R |  qsub -N pic_beta -cwd -pe local 10 -R y -l mf=20G,h_vmem=20G
