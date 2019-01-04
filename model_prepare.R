library(Biostrings)
library(stringr)
library(dplyr)
library(data.table)
library(rtracklayer)
#library(bigmemory)
#library(bigalgebra)


bim <- data.frame(fread("/dcl01/leek/data/sra/sara_genotypes/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6.bim"))
colnames(bim) <- c("chr","variant_id","dummy","pos","minor","major")

## system(liftover.R;liftover.sh) ## get gtex_hg38.bed
gtex_hg38 <- data.frame(fread("gtex_hg38.bed"))
colnames(gtex_hg38) <- c("chr","start","end","id")
gtex_hg38_gr <- GRanges(seqnames= gtex_hg38$chr,
                        ranges=IRanges(gtex_hg38$start,gtex_hg38$start),
                        strand="*",
                        id = gtex_hg38$id)
rm(gtex_hg38)


model_prepare <- function(chr_idx){
  score_A <- c()
  score_C <- c()
  score_T <- c()
  score_G <- c()
  ref_nucs <- c()
  filter <- c()

  ## for loop going over chromosomes 1~22
  for (j in chr_idx){

    chr = paste0("chr",j)

    load(paste0("./score_matrix/filter_gr_",chr,".rda")) # load: filter_gr
    filter_gr <- do.call(c,filter_gr)
    overlap_loci <- findOverlaps(filter_gr,gtex_hg38_gr)
    filter_gr <- filter_gr[queryHits(overlap_loci)]

    load(paste0("./score_matrix/score_data_matrix_A_",chr,".rda")) # load: score_data_matrix_A
    score_data_matrix_A <- score_data_matrix_A[queryHits(overlap_loci),]

    load(paste0("./score_matrix/score_data_matrix_T_",chr,".rda")) # load: score_data_matrix_T
    score_data_matrix_T <- score_data_matrix_T[queryHits(overlap_loci),]

    load(paste0("./score_matrix/score_data_matrix_C_",chr,".rda")) # load: score_data_matrix_C
    score_data_matrix_C <- score_data_matrix_C[queryHits(overlap_loci),]

    load(paste0("./score_matrix/score_data_matrix_G_",chr,".rda")) # load: score_data_matrix_G
    score_data_matrix_G <- score_data_matrix_G[queryHits(overlap_loci),]

    load(paste0("./score_matrix/score_data_matrix_bw_",chr,".rda")) # load: score_data_matrix_bw
    score_data_matrix_bw <- score_data_matrix_bw[queryHits(overlap_loci),]


    snp_position_hg38 <- filter_gr %>% GPos() %>% pos()
    sequence_chr_hg38 <- as.character(readDNAStringSet(paste0("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/",chr,".fa.gz")))
    ref_nucs_chr <- str_sub(sequence_chr_hg38, start = snp_position_hg38, end=snp_position_hg38)
    ref_score_matrix <- score_data_matrix_bw - score_data_matrix_A - score_data_matrix_C - score_data_matrix_G - score_data_matrix_T


    for (i in c(1:length(ref_nucs_chr))){
      if (ref_nucs_chr[i]=="A") {score_data_matrix_A[i,] <- ref_score_matrix[i,]}
      else if (ref_nucs_chr[i]=="T") {score_data_matrix_T[i,] <- ref_score_matrix[i,]}
      else if (ref_nucs_chr[i]=="C") {score_data_matrix_C[i,] <- ref_score_matrix[i,]}
      else if (ref_nucs_chr[i]=="G") {score_data_matrix_G[i,] <- ref_score_matrix[i,]}

    } ## end for loop to go through ref_nucs_chr

    rm(ref_score_matrix)
    gc()

    filter <- c(filter, filter_gr)
    score_A <- rbind(score_A, score_data_matrix_A)
    score_C <- rbind(score_C, score_data_matrix_C)
    score_G <- rbind(score_G, score_data_matrix_G)
    score_T <- rbind(score_T, score_data_matrix_T)
    ref_nucs <- c(ref_nucs, ref_nucs_chr)
  } ## end inner for loop to go through 1:22 chromosomes

  save(filter, file = "./score_matrix/score/filter.rda", compress = TRUE)
  save(score_A, file = "./score_matrix/score/score_A.rda", compress = TRUE)
  save(score_T, file = "./score_matrix/score/score_T.rda", compress = TRUE)
  save(score_C, file = "./score_matrix/score/score_C.rda", compress = TRUE)
  save(score_G, file = "./score_matrix/score/score_G.rda", compress = TRUE)
  save(ref_nucs, file = "./score_matrix/score/ref_nucs.rda", compress = TRUE)


} # end function model_prepare(chr_idx)

model_prepare(c(1:22))


load("./score_matrix/score/filter.rda")
load("./score_matrix/score/score_A.rda")
load("./score_matrix/score/score_T.rda")
load("./score_matrix/score/score_C.rda")
load("./score_matrix/score/score_G.rda")
load("./score_matrix/score/ref_nucs.rda")


### find known genotypes using bim file

filter <- do.call(c,filter)
overlap_loci <- findOverlaps(filter,gtex_hg38_gr)
filter <- filter[queryHits(overlap_loci)]
lftover_id <- gtex_hg38_gr[subjectHits(overlap_loci)]$id
snp_major_minor <- bim[lftover_id, c("major","minor")]

select_snp_name <- bim[lftover_id, "variant_id"]
write.table(select_snp_name,file="select_snp.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

## use select_snp.txt to split ped file and only get thoese snps out in convert_bed2ped.sh
system("bash plink --noweb --bfile /dcl01/leek/data/sra/sara_genotypes/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6 --extract select_snp.txt --recode --transpose --tab --out /dcl01/leek/data/sra/sara_genotypes/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6_select_blood
")

tped <- data.frame(fread("/dcl01/leek/data/sra/sara_genotypes/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6_select_blood.tped")) ## tped: row is positions; column is samples
fam <-  data.frame(fread("/dcl01/leek/data/sra/sara_genotypes/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6.fam"))


load("./rda/sra_gtex_blood.rda")
temp_subject<- sapply(fam[,1], function(x) strsplit(x,"-")[[1]][2])
fam_subject <- paste0("GTEX-",temp_subject) ## has some special tissues not starting at GTEX. eg. NA12878_C
blood_genotype_df <- tped[,-c(1:4)]
matched_subject_id <- match(unique(sra_gtex_blood$SUBJID),fam_subject)
blood_genotype_df <- blood_genotype_df[, matched_subject_id]
colnames(blood_genotype_df) <- unique(sra_gtex_blood$SUBJID)

known_genotype_blood <- matrix(rep(NA,nrow(blood_genotype_df)*ncol(blood_genotype_df)), nrow = nrow(blood_genotype_df))

for (i in 1:nrow(blood_genotype_df)){
  for (j in 1:ncol(blood_genotype_df)){
    if (blood_genotype_df[i,j] == paste0(snp_major_minor[i,]$major, " ", snp_major_minor[i,]$minor) | blood_genotype_df[i,j] == paste0(snp_major_minor[i,]$minor, " ", snp_major_minor[i,]$major)){
      known_genotype_blood[i,j] <- 1
    }else if (blood_genotype_df[i,j] == paste0(snp_major_minor[i,]$major, " ", snp_major_minor[i,]$major)){
      known_genotype_blood[i,j] <- 0
    }else if(blood_genotype_df[i,j] == paste0(snp_major_minor[i,]$minor, " ", snp_major_minor[i,]$minor)){
      known_genotype_blood[i,j] <- 2
    }else{
      known_genotype_blood[i,j] <- NA
    }
  }
}
colnames(known_genotype_blood) <- colnames(blood_genotype_df)
save(known_genotype_blood, file = "./rda/known_genotype_blood.rda", compress = TRUE)

## only subset rows with complete genotype values for all 456 geuvadis samples
no_na_snp_idx <- complete.cases(known_genotype_blood) # complete.cases(): check which rows contain no missing genotype values
chrlist <- chrlist[no_na_snp_idx]
known_genotype_blood <- known_genotype_blood[no_na_snp_idx,]
poslist <- poslist[no_na_snp_idx]
score_A <- score_A[no_na_snp_idx,]
score_T <- score_T[no_na_snp_idx,]
score_C <- score_C[no_na_snp_idx,]
score_G <- score_G[no_na_snp_idx,]
snp_major_minor <- snp_major_minor[no_na_snp_idx,]
snp_ref_nucs <- snp_ref_nucs[no_na_snp_idx]

save(chrlist,file = "./rda/chrlist.rda",compress = TRUE)
save(poslist,file = "./rda/poslist.rda", compress = TRUE)
save(score_A, file = "./rda/score_A.rda", compress = TRUE)
save(score_T, file = "./rda/score_T.rda", compress = TRUE)
save(score_C, file = "./rda/score_C.rda", compress = TRUE)
save(score_G, file = "./rda/score_G.rda", compress = TRUE)
save(known_genotype_blood, file = "./rda/known_genotype_blood.rda", compress = TRUE)
save(snp_major_minor, file = "./rda/snp_major_minor.rda", compress = TRUE)
save(snp_ref_nucs, file = "./rda/snp_ref_nucs.rda", compress = TRUE)

q(save = "no")
# shell script: echo R CMD BATCH model_prepare.R | qsub -N data_pre -cwd -l mf=5G,h_vmem=5G
