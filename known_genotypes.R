
################ get genotypes for gtex sample_individuals (191 individuals)
library(data.table)
library(rtracklayer)

load("./rda/chrlist.rda")
load("./rda/poslist.rda")
load("./rda/snp_ref_nucs.rda")
load("./rda/snp_major_minor.rda")
model_gr <- GRanges(seqnames=chrlist,ranges=IRanges(poslist,poslist),strand="*")
bim <- data.frame(fread("/dcl01/leek/data/sra/sara_genotypes/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6.bim"))
colnames(bim) <- c("chr","variant_id","dummy","pos","minor","major")

#gtex_hg38 <- read.table("gtex_hg38.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
gtex_hg38 <- data.frame(fread("gtex_hg38.bed"))
colnames(gtex_hg38) <- c("chr","start","end","id")
gtex_hg38_gr <- GRanges(seqnames= gtex_hg38$chr,
                        ranges=IRanges(gtex_hg38$start,gtex_hg38$start),
                        strand="*",
                        id = gtex_hg38$id)
rm(gtex_hg38)
overlap_loci <- findOverlaps(model_gr,gtex_hg38_gr)
save(overlap_loci, file = "./rda/overlap_loci.rda", compress = TRUE)

model_mm <- snp_major_minor[queryHits(overlap_loci),] ## 1063 string mismatches to gtex_hg38_mm

id <- gtex_hg38_gr[subjectHits(overlap_loci)]$id

gtex_hg38_mm <- bim[id,c("major","minor")] ## 1063 string mismatches to model_mm
select_snp_name <- bim[id,"variant_id"]
# write.table(select_snp_name,file="select_snp.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
## use select_snp.txt to split ped file and only get thoese snps out in convert_bed2ped.sh
# system("convert_bed2ped.sh")

tped <- data.frame(fread("/dcl01/leek/data/sra/sara_genotypes/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6_select.tped"))
fam <-  data.frame(fread("/dcl01/leek/data/sra/sara_genotypes/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6.fam"))


temp_subject<- sapply(fam[,1], function(x) strsplit(x,"-")[[1]][2])
fam_subject <- paste0("GTEX-",temp_subject) ## has some special tissues not starting at GTEX

known_genotype_df <- tped[,-c(1:4)]


load("./rda/sra_gtex_df.rda")
matched_subject_id <- match(unique(sra_gtex_df$SUBJID),fam_subject)
known_genotype_df <- known_genotype_df[, matched_subject_id]

known_genotype_matrix <- matrix(rep(NA,nrow(known_genotype_df)*ncol(known_genotype_df)), nrow = nrow(known_genotype_df))

for (i in 1:nrow(known_genotype_df)){
  for (j in 1:ncol(known_genotype_df)){
    if (known_genotype_df[i,j] == paste0(gtex_hg38_mm[i,]$major, " ", gtex_hg38_mm[i,]$minor) | known_genotype_df[i,j] == paste0(gtex_hg38_mm[i,]$minor, " ", gtex_hg38_mm[i,]$major)){
      known_genotype_matrix[i,j] <- 1
    }else if (known_genotype_df[i,j] == paste0(gtex_hg38_mm[i,]$major, " ", gtex_hg38_mm[i,]$major)){
      known_genotype_matrix[i,j] <- 0
    }else if(known_genotype_df[i,j] == paste0(gtex_hg38_mm[i,]$minor, " ", gtex_hg38_mm[i,]$minor)){
      known_genotype_matrix[i,j] <- 2
    }else{
      known_genotype_matrix[i,j] <- NA
    }
  }
}

save(known_genotype_matrix, file = "./rda/known_genotype_matrix.rda", compress = TRUE)
q(save = "no")

# script: echo R CMD BATCH known_genotypes.R | qsub -N known_geno -cwd -l mf=8G,h_vmem=8G
