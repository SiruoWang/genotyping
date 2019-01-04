
library(ranger)
load("./rda/score_A.rda")
load("./rda/score_C.rda")
load("./rda/score_G.rda")
load("./rda/score_T.rda")
load("./rda/conv_geno.rda") ## 16 unique genotypes: "TT" "AA" "CC" "AG" "GG" "CT" "GC" "GA" "TA" "TC" "CG" "TG" "CA" "AC" "GT" "AT"

### split data into prior calculation set, model fitting set, and validation set
set.seed(2018)
prior_id = c(1:ceiling(ncol(conv_geno)/3))
model_id = c((prior_id[length(prior_id)] + 1) : (prior_id[length(prior_id)] +ceiling(ncol(conv_geno)/3) ))
val_id = c((model_id[length(model_id)] + 1) : ncol(conv_geno))

#row_id = sample(c(1:nrow(conv_geno)),30000)

### fit random forest model using ranger package

total = score_A + score_C + score_G + score_T

geno = c(conv_geno[,model_id])
ratioA = c(score_A[,model_id] / total[,model_id])
ratioC = c(score_C[,model_id] / total[,model_id])
ratioG = c(score_G[,model_id] / total[,model_id])
ratioT = c(score_T[,model_id] / total[,model_id])

ratioA[is.na(ratioA)] = 0
ratioC[is.na(ratioC)] = 0
ratioG[is.na(ratioG)] = 0
ratioT[is.na(ratioT)] = 0

rf_model_df = data.frame(geno = geno, ratioA = ratioA, ratioC = ratioC, ratioG = ratioG, ratioT = ratioT)
null_model = ranger(geno ~ ratioA + ratioC + ratioG + ratioT, data = rf_model_df, num.trees = 100, write.forest = TRUE, num.threads = 20)


save(null_model, file = "./rda/null_model.rda", compress = TRUE)

rm(rf_model_df)
gc()

#### make predictions using null_model on the validation set
ratioA = c(score_A[,val_id] / total[,val_id])
ratioC = c(score_C[,val_id] / total[,val_id])
ratioG = c(score_G[,val_id] / total[,val_id])
ratioT = c(score_T[,val_id] / total[,val_id])

ratioA[is.na(ratioA)] = 0
ratioC[is.na(ratioC)] = 0
ratioG[is.na(ratioG)] = 0
ratioT[is.na(ratioT)] = 0

rf_pred_df <- data.frame(ratioA = ratioA, ratioC = ratioC, ratioG = ratioG, ratioT = ratioT)
prediction <- predict(null_model, dat = rf_pred_df, num.threads = 20)

pred_geno <- prediction$predictions
pred_geno <- matrix(pred_geno, nrow = nrow(score_A))

save(pred_geno,file = "./rda/null_pred_geno.rda", compress = TRUE)
q(save = "no")

#shell script: echo R CMD BATCH null_model.R | qsub -N null_model -cwd -l mf=50G,h_vmem=50G
