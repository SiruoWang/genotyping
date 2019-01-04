
load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/snp_major_minor.rda")
load("/users/swang1/2017_recount_genotype/geuvadis/genotype_call/rda/snp_genotype.rda")


conv_geno <- matrix(rep(NA, nrow(snp_genotype) * ncol(snp_genotype)), nrow = nrow(snp_genotype))

for (i in 1:nrow(snp_genotype)){
  for (j in 1:ncol(snp_genotype)){
    if (snp_genotype[i,j] == 1){
      conv_geno[i,j] <- paste0(snp_major_minor[i,1], snp_major_minor[i,2])
    }else if (snp_genotype[i,j] == 0){
      conv_geno[i,j] <- paste0(snp_major_minor[i,1], snp_major_minor[i,1])
    }else if(snp_genotype[i,j] == 2){
      conv_geno[i,j] <- paste0(snp_major_minor[i,2], snp_major_minor[i,2])
    }else{
      conv_geno[i,j] <- NA
    }
  }
}
colnames(conv_geno) = colnames(snp_genotype)

save(conv_geno, file = "conv_geno.rda", compress = TRUE)

# shell script: echo R CMD BATCH convert_known_geno.R | qsub -N conv_geno -cwd -l mf=3G,h_vmem=3G
