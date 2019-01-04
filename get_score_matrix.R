args <- commandArgs(trailingOnly=FALSE)
chr <- args[length(args)]
chr <- sub("-","",chr)
print(chr)


#source("gene_annotation_hg38.R")
load("./rda/gene_info_hg38_chr_gr.rda") # gene_info_hg38_chr_gr.rda saved in gene_annotation_hg38.R
library(rtracklayer)
library(matrixStats)


# get_bp(tempA) separates ranges of Granges objects with width greater than 1 into Granges objects with width = 1 using GPos()
get_bp <- function(temp,chr_in){
  if (length(temp) == 0) {
    return(temp)
  }else{
    temp_pos <- GPos(temp)
    # create GRange object gr with the three keys
    gr <- GRanges(seqnames=chr_in,ranges=IRanges(pos(temp_pos),pos(temp_pos)),
                  strand=strand(temp_pos))
    overlap_idx <- findOverlaps(gr,temp)
    # add score to the new gr object
    gr$score[queryHits(overlap_idx)] <- temp$score[subjectHits(overlap_idx)]
    return(gr)
  }
}

output_A <- read.table(file="./path/path_A_blood.txt",stringsAsFactors = FALSE)
output_T <- read.table(file="./path/path_C_blood.txt",stringsAsFactors = FALSE)
output_C <- read.table(file="./path/path_G_blood.txt",stringsAsFactors = FALSE)
output_G <- read.table(file="./path/path_T_blood.txt",stringsAsFactors = FALSE)
output_bw <- read.table(file="./path/path_bw_blood.txt",stringsAsFactors = FALSE)

load("./rda/sra_gtex.rda")
temp_Run <- sapply(output_bw[,1], function(x) strsplit(x,"/")[[1]][8])
Run <- sapply(temp_Run, function(x) strsplit(x,"_")[[1]][1])
sra_gtex_blood <- sra_gtex[match(Run, sra_gtex[,"Run"]),] ## runs in sra_gtex_blood (row) is the same order as path_bw_blood.txt (row)
#save(sra_gtex_blood, file = "./rda/sra_gtex_blood.rda", compress = TRUE)


create_testIndex <- function(num, break_num){
  Index_list <- list()
  group_idx <- sample(cut(seq(1,num),breaks=break_num,label=FALSE))
  for (i in 1:break_num){
    Index_list[[i]] <- which(group_idx == i)
  }
  return(Index_list)
}

set.seed(2018)
subjid <- unique(sra_gtex_blood$SUBJID)
random_index <- create_testIndex(length(subjid),3)
subjid_train <- subjid[unlist(random_index[1:2])]
subjid_val <- subjid[unlist(random_index[3])]
#save(subjid_train, file = "./rda/subjid_train.rda", compress = TRUE)
#save(subjid_val, file = "./rda/subjid_val.rda", compress = TRUE)

run_train <- which(sra_gtex_blood$SUBJID %in% subjid_train)
run_names <- sra_gtex_blood[run_train,"Run"]

chr_size_table <- read.table("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes")
chr_size <- chr_size_table[which(chr_size_table[,1] == chr),2]
#chunks <- parallel::splitIndices(chr_size, ncl = 100)
print(chr_size)

score_data_matrix_bw <- c()
score_data_matrix_A <- c()
score_data_matrix_C <- c()
score_data_matrix_G <- c()
score_data_matrix_T <- c()
filter_gr <- c()

for (i in 1:20){
  print(i)
  which <- GRanges(seqnames=chr,ranges=IRanges((i-1)*ceiling(chr_size/20)+1,i*ceiling(chr_size/20)),strand="*")
  print(which)

  togetherA <- lapply(output_A[run_train,1],function(x) {
    tempA <- import(x,format="bigwig", which = which)
    tempA <- tempA[tempA$score != 0]
    tempA <- get_bp(tempA, chr)
  })


  togetherT <- lapply(output_T[run_train,1],function(x) {
    tempT <- import(x,format="bigwig", which = which)
    tempT <- tempT[tempT$score != 0]
    tempT <- get_bp(tempT, chr)
  })


  togetherC <- lapply(output_C[run_train,1],function(x) {
    tempC <- import(x,format="bigwig", which = which)
    tempC <- tempC[tempC$score != 0]
    tempC <- get_bp(tempC, chr)
  })


  togetherG <- lapply(output_G[run_train,1],function(x) {
    tempG <- import(x,format="bigwig", which = which)
    tempG <- tempG[tempG$score != 0]
    tempG <- get_bp(tempG, chr)
  })

  togetherbw <- lapply(output_bw[run_train,1],function(x) {
    tempbw <- import(x,format="bigwig", which = which)
    tempbw <- tempbw[tempbw$score != 0]
    tempbw <- get_bp(tempbw, chr)
  })

  together_ATCG <- c(togetherA, togetherT, togetherC, togetherG, togetherbw)
  # apply do.call(c, files) to multiple GRange files to collapse into one GRange files
  new_together_ATCG <- do.call(c, together_ATCG)
  # reduce returns an object of the same type as x containing reduced ranges for each distinct (seqname, strand) pairing.
  # set reduce range width to be 1, so ranges of GRange objects (width == 1) do not collapse into width > 1
  reduce_together_ATCG <- reduce(new_together_ATCG,min.gapwidth=0L)
  # match reduce_together_ATCG ranges with gene positions, and only subset positions that contain genes
  overlap_gene_index <- findOverlaps(gene_info_hg38_chr_gr,reduce_together_ATCG)
  # genes overlap, so we can only take out the unique subjectHits of reduce_together_ATCG position
  unique_overlap_gene_index <- unique(subjectHits(overlap_gene_index))
  reduce_together_ATCG <- reduce_together_ATCG[unique_overlap_gene_index]

  #save(reduce_together_ATCG, file = paste0("reduce_together_ATCG_",chr,".rda"), compress = TRUE)
  rm(together_ATCG)
  rm(new_together_ATCG)
  rm(overlap_gene_index)
  rm(unique_overlap_gene_index)
  gc()

  if (length(reduce_together_ATCG) != 0){
    # score_data_matrix_ATCG has the same locus potion order as in the reduce_together_ATCG (row order matches)
    score_bw <- matrix(data = 0, nrow = length(reduce_together_ATCG),ncol = length(run_train))
    score_bw <- apply(score_bw,c(1,2),as.integer)

    for(i in 1:length(togetherbw)){
      overlap_loci <- findOverlaps(reduce_together_ATCG,togetherbw[[i]])
      score_bw[queryHits(overlap_loci),i] <- togetherbw[[i]]$score[subjectHits(overlap_loci)]
    }

    ## overlap reduce_together_ATCG with known genotype position first (if we have trustful genotypes already), and then filter by row median

    total_reads_median <- rowMedians(score_bw)
    filter_row_idx <- which(total_reads_median > 5)
    filter <- reduce_together_ATCG[filter_row_idx]
    filter_gr <- c(filter_gr, filter)

    rm(reduce_together_ATCG)
    gc()

    score_bw <- score_bw[filter_row_idx,]
    score_data_matrix_bw <- rbind(score_data_matrix_bw,score_bw)
    print("rbind score_data_matrix_bw done")
    rm(score_bw)
    rm(togetherbw)
    gc()


    score_A <- matrix(data = 0, nrow = length(filter),ncol = length(run_train))
    score_A <- apply(score_A,c(1,2),as.integer)

    for(i in 1:length(togetherA)){
      overlap_loci <- findOverlaps(filter,togetherA[[i]])
      score_A[queryHits(overlap_loci),i] <- togetherA[[i]]$score[subjectHits(overlap_loci)]
    }

    score_data_matrix_A <- rbind(score_data_matrix_A,score_A)
    print("rbind score_data_matrix_A done")
    rm(score_A)
    rm(togetherA)
    gc()

    score_T <- matrix(data = 0, nrow = length(filter),ncol = length(run_train))
    score_T <- apply(score_T,c(1,2),as.integer)

    for(i in 1:length(togetherT)){
      overlap_loci <- findOverlaps(filter,togetherT[[i]])
      score_T[queryHits(overlap_loci),i] <- togetherT[[i]]$score[subjectHits(overlap_loci)]
    }

    score_data_matrix_T <- rbind(score_data_matrix_T,score_T)
    print("rbind score_data_matrix_T done")
    rm(score_T)
    rm(togetherT)
    gc()


    score_C <- matrix(data = 0, nrow = length(filter),ncol = length(run_train))
    score_C <- apply(score_C,c(1,2),as.integer)

    for(i in 1:length(togetherC)){
      overlap_loci <- findOverlaps(filter,togetherC[[i]])
      score_C[queryHits(overlap_loci),i] <- togetherC[[i]]$score[subjectHits(overlap_loci)]
    }

    score_data_matrix_C <- rbind(score_data_matrix_C,score_C)
    print("rbind score_data_matrix_C done")
    rm(score_C)
    rm(togetherC)
    gc()

    score_G <- matrix(data = 0, nrow = length(filter),ncol = length(run_train))
    score_G <- apply(score_G,c(1,2),as.integer)

    for(i in 1:length(togetherG)){
      overlap_loci <- findOverlaps(filter,togetherG[[i]])
      score_G[queryHits(overlap_loci),i] <- togetherG[[i]]$score[subjectHits(overlap_loci)]
    }

    score_data_matrix_G <- rbind(score_data_matrix_G,score_G)
    print("rbind score_data_matrix_G done")
    rm(score_G)
    rm(togetherG)
    gc()
 }

}

colnames(score_data_matrix_bw) <- colnames(score_data_matrix_A) <- colnames(score_data_matrix_T) <- colnames(score_data_matrix_C) <- colnames(score_data_matrix_G) <- run_names

save(filter_gr, file = paste0("./score_matrix/filter_gr_",chr,".rda"), compress = TRUE)
save(score_data_matrix_bw, file = paste0("./score_matrix/score_data_matrix_bw_",chr,".rda"),compress = TRUE)
save(score_data_matrix_A, file = paste0("./score_matrix/score_data_matrix_A_",chr,".rda"),compress = TRUE)
save(score_data_matrix_T, file = paste0("./score_matrix/score_data_matrix_T_",chr,".rda"),compress = TRUE)
save(score_data_matrix_C, file = paste0("./score_matrix/score_data_matrix_C_",chr,".rda"),compress = TRUE)
save(score_data_matrix_G, file = paste0("./score_matrix/score_data_matrix_G_",chr,".rda"),compress = TRUE)

q(save="no")
