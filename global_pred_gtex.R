library(rtracklayer)
library(ranger)
library(parallel)



load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/chrlist.rda")
load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/poslist.rda")
load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/snp_ref_nucs.rda")
load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/snp_major_minor.rda")
load("/users/swang1/2017_recount_genotype/model1/global/rda/null_model.rda")

load("/users/swang1/2017_recount_genotype/gtex/rda/sra_gtex_blood.rda")
load("/users/swang1/2017_recount_genotype/gtex/rda/subjid_val.rda")



prediction_prepare <- function(){
  sra_path_list <- list()

  output_A <- read.table(file="/users/swang1/2017_recount_genotype/gtex/path/path_A_blood.txt",stringsAsFactors = FALSE)
  output_C <- read.table(file="/users/swang1/2017_recount_genotype/gtex/path/path_G_blood.txt",stringsAsFactors = FALSE)
  output_G <- read.table(file="/users/swang1/2017_recount_genotype/gtex/path/path_T_blood.txt",stringsAsFactors = FALSE)
  output_T <- read.table(file="/users/swang1/2017_recount_genotype/gtex/path/path_C_blood.txt",stringsAsFactors = FALSE)
  output_bw <- read.table(file="/users/swang1/2017_recount_genotype/gtex/path/path_bw_blood.txt",stringsAsFactors = FALSE)

  #run_val_id <- which(sra_gtex_blood$SUBJID %in% subjid_val[1:10])
  run_val_id <- match(subjid_val[1:10],sra_gtex_blood$SUBJID)

  sra_path_list[[1]] <- output_A[run_val_id,1]
  sra_path_list[[2]] <- output_C[run_val_id,1]
  sra_path_list[[3]] <- output_G[run_val_id,1]
  sra_path_list[[4]] <- output_T[run_val_id,1]
  sra_path_list[[5]] <- output_bw[run_val_id,1]
  sra_path_list[[6]] <- sra_gtex_blood[run_val_id,"Run"]

  return(sra_path_list)
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
  }

  #genotype_per_sample <- rep(NA, length(poslist)) # length(poslist)

  ratioA = scoreA_list / scorebw_list
  ratioC = scoreC_list / scorebw_list
  ratioG = scoreG_list / scorebw_list
  ratioT = scoreT_list / scorebw_list

  ratioA[is.na(ratioA)] = 0
  ratioC[is.na(ratioC)] = 0
  ratioG[is.na(ratioG)] = 0
  ratioT[is.na(ratioT)] = 0

  rf_pred_df <- data.frame(ratioA = ratioA, ratioC = ratioC, ratioG = ratioG, ratioT = ratioT)
  prediction <- predict(null_model, dat = rf_pred_df, num.threads = 20)

  pred_geno <- prediction$predictions

  return(pred_geno)

}

sra_path_list <- prediction_prepare()
#run_names <- sra_path_list[[6]]
#save(run_names, file = "./rda/run_names.rda", compress = TRUE)

predict_genotype <- lapply(c(1:10), function(x) genotype_prediction(x,sra_path_list))
output <- t(matrix(unlist(predict_genotype), nrow = 10, byrow = TRUE))
colnames(output) = subjid_val[1:10]

save(output, file = "output.rda")
q(save = "no")

## shell script: echo R CMD BATCH geuv_pred_gtex_blood.R | qsub -cwd -pe local 10 -R y -l mf=20G,h_vmem=20G
