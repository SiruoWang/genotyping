library(data.table)
library(dplyr)

bim <- fread("/dcl01/leek/data/sra/GTEX_Genotypes/Merged/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6.bim") %>% data.frame()
colnames(bim) <- c("chr","variant_id","dummy","pos","minor","major")

# for(k in 1:22){
#
#   bim_sub <- bim[which(bim[,1]==k),]
#   hg19_output <- data.frame(chromosome=paste0("chr",bim_sub$chr),start=bim_sub$pos,end=bim_sub$pos+1,id=c(1:nrow(bim_sub)))
#
#   options(scipen=999)
#   write.table(hg19_output,file=paste0("gtex_hg19_chr",k,".bed"),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
#
#   system(paste0("bash liftover.sh gtex_hg19_chr",k,".bed ","gtex_hg38_chr",k,".bed"))
# }


hg19_output <- data.frame(chromosome=paste0("chr",bim$chr),start=bim_sub$pos,end=bim_sub$pos+1,id=c(1:nrow(bim)))

options(scipen=999)
write.table(hg19_output,file="gtex_hg19.bed"),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)

system("bash liftover.sh gtex_hg19.bed gtex_hg38.bed")

#q(save="no")
# qsub -cwd -l mem_free=5G,h_vmem=5G run_liftover.sh
# script: echo R CMD BATCH liftover.R | qsub -N liftover -cwd -l mf=8G,h_vmem=8G
