setwd("~/Documents/2017_recount_genotype/gtex")

library(data.table)
library(dplyr)
library(rtracklayer)
library(parallel)
#ped <- fread("GTEX_5M_Genotype_Data_matched.chr22.ped")
#bim <- fread("GTEX_5M_Genotype_Data_matched.chr22.bim")
fam <- fread("/dcl01/leek/data/sra/sara_genotypes/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6.fam")

load("sample_individuals.Rdata")
load("gtexmetadata.rda")
gtex_meta = gtexmetadata
gtex_meta = cbind(gtex_meta,usegtex)
rm(gtexmetadata,usegtex)
usegtex = gtex_meta$usegtex
pheno = gtex_meta
pheno = pheno[usegtex,]

pheno_subset <- pheno[,c("Run", "SUBJID","Histological_Type", "Body_Site")]

tempid_1<- sapply(as.matrix(fam[,1]), function(x) strsplit(x,"-")[[1]][2])
sra_subjid <- paste0("GTEX-",tempid_1)
sra_gtex <- pheno_subset[which(pheno_subset$SUBJID %in% sra_subjid),]
save(sra_gtex, file = "./rda/sra_gtex.rda", compress = TRUE)


## full paths for all gtex samples
path_A <- paste0("find /dcl01/leek/data/gtex_v2/batch_*/coverage_bigwigs/", sra_gtex[,"Run"], "*?[0-9].A.bw >> path_A.txt")
write.table(path_A, file = "path_A.sh", row.names = FALSE, col.names = FALSE, quote=FALSE)

path_C <- paste0("find /dcl01/leek/data/gtex_v2/batch_*/coverage_bigwigs/", sra_gtex[,"Run"], "*?[0-9].C.bw >> path_C.txt")
write.table(path_C, file = "path_C.sh", row.names = FALSE, col.names = FALSE, quote=FALSE)

path_G <- paste0("find /dcl01/leek/data/gtex_v2/batch_*/coverage_bigwigs/", sra_gtex[,"Run"], "*?[0-9].G.bw >> path_G.txt")
write.table(path_G, file = "path_G.sh", row.names = FALSE, col.names = FALSE, quote=FALSE)

path_T <- paste0("find /dcl01/leek/data/gtex_v2/batch_*/coverage_bigwigs/", sra_gtex[,"Run"], "*?[0-9].T.bw >> path_T.txt")
write.table(path_T, file = "path_T.sh", row.names = FALSE, col.names = FALSE, quote=FALSE)

path_bw <- paste0("find /dcl01/leek/data/gtex_v2/batch_*/coverage_bigwigs/", sra_gtex[,"Run"], "*?[0-9].bw >> path_bw.txt")
write.table(path_bw, file = "path_bw.sh", row.names = FALSE, col.names = FALSE, quote=FALSE


tissue_paths <- function(tissue){
  ## full paths for subset gtex samples of certain tissue
  sub <- sra_gtex[which(sra_gtex$Histological_Type == tissue),]

  path_A_skin <- paste0("find /dcl01/leek/data/gtex_v2/batch_*/coverage_bigwigs/", sub[,"Run"], "*?[0-9].A.bw >> path_A_",tissue,".txt")
  write.table(path_A_skin, file = "path_A_skin.sh", row.names = FALSE, col.names = FALSE, quote=FALSE)

  path_C_skin <- paste0("find /dcl01/leek/data/gtex_v2/batch_*/coverage_bigwigs/", sub[,"Run"], "*?[0-9].C.bw >> path_C_",tissue,".txt")
  write.table(path_C_skin, file = "path_C_skin.sh", row.names = FALSE, col.names = FALSE, quote=FALSE)

  path_G_skin <- paste0("find /dcl01/leek/data/gtex_v2/batch_*/coverage_bigwigs/", sub[,"Run"], "*?[0-9].G.bw >> path_G_",tissue,".txt")
  write.table(path_G_skin, file = "path_G_skin.sh", row.names = FALSE, col.names = FALSE, quote=FALSE)

  path_T_skin <- paste0("find /dcl01/leek/data/gtex_v2/batch_*/coverage_bigwigs/", sub[,"Run"], "*?[0-9].T.bw >> path_T_",tissue,".txt")
  write.table(path_T_skin, file = "path_T_skin.sh", row.names = FALSE, col.names = FALSE, quote=FALSE)

  path_bw_skin <- paste0("find /dcl01/leek/data/gtex_v2/batch_*/coverage_bigwigs/", sub[,"Run"], "*?[0-9].bw >> path_bw_",tissue,".txt")
  write.table(path_bw_skin, file = "path_bw_skin.sh", row.names = FALSE, col.names = FALSE, quote=FALSE)
}

tissue_paths("skin")
tissue_paths("blood")
tissue_paths("brain")

q(save = "no")
