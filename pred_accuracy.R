library(caret)
library(rtracklayer)

load("./rda/predict_genotype.rda")
load("./rda/known_genotype_matrix.rda")
load("./rda/sra_gtex_blood.rda")
load("./rda/subjid_val.rda")
load("./rda/run_names.rda")


pred_accuracy <- sapply(subjid_val, function(x) {

  runs <- sra_gtex_blood[which(sra_gtex_blood$SUBJID == x), "Run"]
  predict_run_id <- match(runs, run_names)
  known_col <- which(colnames(known_genotype_matrix) == x)

  compare_runs <- lapply(predict_genotype[predict_run_id], function(pred_geno){
    confusionMatrix(factor(pred_geno, levels = 0:2), factor(known_genotype_matrix[, known_col], levels = 0:2))$overall['Accuracy']
    })
  accuracy <- mean(unlist(compare_runs))

})

subject_accuracy_df <- data.frame(subject = subjid_val, accuracy = pred_accuracy)
save(subject_accuracy_df, file = "./rda/subject_accuracy_df.rda", compress = TRUE)

# script: echo R CMD BATCH gtex_analysis.R | qsub -N analysis -cwd -l mf=4G,h_vmem=4G
