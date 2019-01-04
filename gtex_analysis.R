library(caret)
library(rtracklayer)
load("./rda/sra_gtex_df.rda")
load("./rda/gtex_genotype_prediction.rda")
load("./rda/known_genotype_matrix.rda")
load("./rda/overlap_loci.rda")

pred_accuracy <- sapply(unique(sra_gtex_df$Histological_Type), function(x) {

  tissue_id <- which(sra_gtex_df$Histological_Type == x)
  tissue_df <- sra_gtex_df[tissue_id,]

  ## get loci that have both known genotypes and predicted genotypes using overlap_loci (saved from known_genotypes.R)
  predict_genotype_tissue <- lapply(gtex_genotype_prediction[tissue_id], function(x) {x[queryHits(overlap_loci)]})

  ## Run in nrow(sra_gtex_df) is the same order as in the gtex_genotype_prediction
  ## column in known_genotype_matrix is matched with known_geno_subject (== unique(sra_gtex_df$SUBJID))
  known_geno_subject <- unique(sra_gtex_df$SUBJID) # unique 184 subjects in gtex_genotype_prediction
  uniq_tissue_subject <- unique(tissue_df$SUBJID) # unique 41 subjects from brain tissues

  subj_tissue_accuracy <- sapply(uniq_tissue_subject, function(i) {
    ## find known genotypes for subject i
    known_genotype_i <- known_genotype_matrix[, which(known_geno_subject == i)]

    ## find predicted genotypes for multiple sra runs of subject i
    ## list order of predict_genotype_tissue matches subjects order (row order) of tissue_df
    predicted_genotype_i <- predict_genotype_tissue[which(tissue_df$SUBJID == i)]
    compare_runs <- lapply(predicted_genotype_i, function(x) {
      confusionMatrix(factor(x, levels = 0:2), factor(known_genotype_i, levels = 0:2))$overall['Accuracy']
      })
    mean(unlist(compare_runs))
  })

  mean(subj_tissue_accuracy)

})

tissue_accuracy_df <- data.frame(tissue = unique(sra_gtex_df$Histological_Type), pred_accuracy)
save(tissue_accuracy_df, file = "./rda/tissue_accuracy_df.rda", compress = TRUE)

# script: echo R CMD BATCH gtex_analysis.R | qsub -N analysis -cwd -l mf=4G,h_vmem=4G
