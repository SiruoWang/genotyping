library(rtracklayer)
library(ranger)
library(parallel)



load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/chrlist.rda")
load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/poslist.rda")
load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/snp_ref_nucs.rda")
load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/snp_major_minor.rda")
load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/geuvadis_training_model_ratio.rda")

load("./rda/sra_gtex_blood.rda")
load("./rda/subjid_val.rda")



prediction_prepare <- function(){
  sra_path_list <- list()

  output_A <- read.table(file="./path/path_A_blood.txt",stringsAsFactors = FALSE)
  output_C <- read.table(file="./path/path_G_blood.txt",stringsAsFactors = FALSE)
  output_G <- read.table(file="./path/path_T_blood.txt",stringsAsFactors = FALSE)
  output_T <- read.table(file="./path/path_C_blood.txt",stringsAsFactors = FALSE)
  output_bw <- read.table(file="./path/path_bw_blood.txt",stringsAsFactors = FALSE)

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

  genotype_per_sample <- rep(NA, length(poslist)) # length(poslist)

  for (i in c(1:length(poslist))){
    distance_cluster <- unlist(geuvadis_training_model[[i]][1])
    snp_subset_A <- scoreA_list[distance_cluster]
    snp_subset_C <- scoreC_list[distance_cluster]
    snp_subset_G <- scoreG_list[distance_cluster]
    snp_subset_T <- scoreT_list[distance_cluster]

    total <- sum(c(snp_subset_A, snp_subset_C, snp_subset_G, snp_subset_T))

    if (total == 0){
      genotype <- NA
    }else{
      snpcombine <- data.frame(nucA = t(snp_subset_A),# A counts for all snps in this distance cluster
                               nucC = t(snp_subset_C),# C counts for all snps in this distance cluster
                               nucG = t(snp_subset_G),# G counts for all snps in this distance cluster
                               nucT = t(snp_subset_T))# T counts for all snps in this distance cluster

      predict_genotype <- predict(geuvadis_training_model[[i]][2][[1]], dat = snpcombine)
      genotype <- predict_genotype$predictions

      # if (genotype_num == 0) {genotype <- c(snp_major_minor[i,"major"], snp_major_minor[i,"major"])
      # }else if (genotype_num == 1) {genotype <- c(snp_major_minor[i,"major"], snp_major_minor[i,"minor"])
      # }else if (genotype_num == 2) {genotype <- c(snp_major_minor[i,"minor"], snp_major_minor[i,"minor"])}

    } ## close else(total >0)

    genotype_per_sample[i] <- as.character(genotype)

  } # close for loop to go over 50531 snps
  return(genotype_per_sample)
}

sra_path_list <- prediction_prepare()
#run_names <- sra_path_list[[6]]
#save(run_names, file = "./rda/run_names.rda", compress = TRUE)


predict_genotype <- mclapply(c(1:length(sra_path_list[[1]])),
                                    genotype_prediction,
                                    path = sra_path_list,
                                    mc.cores = 10)
save(predict_genotype, file = "geuv_pred_gtex_genotype.rda", compress = TRUE)
q(save = "no")


## shell script: echo R CMD BATCH geuv_pred_gtex_blood.R | qsub -cwd -pe local 10 -R y -l mf=20G,h_vmem=20G
