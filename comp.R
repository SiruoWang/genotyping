library(ranger)

load("./rda/score_A.rda")
load("./rda/score_C.rda")
load("./rda/score_G.rda")
load("./rda/score_T.rda")
load("./rda/conv_geno.rda")

prior_id = c(1:ceiling(ncol(conv_geno)/3))
model_id = c((prior_id[length(prior_id)] + 1) : (prior_id[length(prior_id)] +ceiling(ncol(conv_geno)/3) ))
val_id = c((model_id[length(model_id)] + 1) : ncol(conv_geno))

### calculate the posterior ratio estimates using alpha.prior and beta.prior
ratio_est <- function(counts, total, alpha, beta){
  return((counts + alpha) / (total + alpha+ beta))
}

total = score_A + score_C + score_G + score_T

## ratio A

alpha.prior <- 0.04189064
beta.prior <- 0.14943819
A_est_modelset <- ratio_est(c(score_A[,model_id]), c(total[,model_id]), alpha.prior, beta.prior)
A_est_valset <- ratio_est(c(score_A[,val_id]), c(total[,val_id]), alpha.prior, beta.prior)
gc()

## ratio C

alpha.prior <- 0.05040109
beta.prior <- 0.12654410
C_est_modelset <- ratio_est(c(score_C[,model_id]), c(total[,model_id]), alpha.prior, beta.prior)
C_est_valset <- ratio_est(c(score_C[,val_id]), c(total[,val_id]), alpha.prior, beta.prior)
gc()

## ratio G

alpha.prior <- 0.04767691
beta.prior <- 0.12520538
G_est_modelset <- ratio_est(c(score_G[,model_id]), c(total[,model_id]), alpha.prior, beta.prior)
G_est_valset <- ratio_est(c(score_G[,val_id]), c(total[,val_id]), alpha.prior, beta.prior)
gc()

## ratio T

alpha.prior <- 0.04291391
beta.prior <- 0.15213301
T_est_modelset <- ratio_est(c(score_T[,model_id]), c(total[,model_id]), alpha.prior, beta.prior)
T_est_valset <- ratio_est(c(score_T[,val_id]), c(total[,val_id]), alpha.prior, beta.prior)
gc()
############################ build model in the model set and make predictions in the val set

### fit random forest model using ranger package

geno = c(conv_geno[,model_id])
ratioA = A_est_modelset
ratioC = C_est_modelset
ratioG = G_est_modelset
ratioT = T_est_modelset

ratioA[is.na(ratioA)] = 0
ratioC[is.na(ratioC)] = 0
ratioG[is.na(ratioG)] = 0
ratioT[is.na(ratioT)] = 0

rf_model_df = data.frame(geno = geno, ratioA = ratioA, ratioC = ratioC, ratioG = ratioG, ratioT = ratioT)
beta_model = ranger(geno ~ ratioA + ratioC + ratioG + ratioT, data = rf_model_df, num.trees = 100, write.forest = TRUE, num.threads = 20)

save(beta_model, file = "./rda/beta_model.rda", compress = TRUE)

rm(rf_model_df)
gc()

#### make predictions using znb_model on the validation set

ratioA = A_est_valset
ratioC = C_est_valset
ratioG = G_est_valset
ratioT = T_est_valset

ratioA[is.na(ratioA)] = 0
ratioC[is.na(ratioC)] = 0
ratioG[is.na(ratioG)] = 0
ratioT[is.na(ratioT)] = 0

rf_pred_df <- data.frame(ratioA = ratioA, ratioC = ratioC, ratioG = ratioG, ratioT = ratioT)
prediction <- predict(beta_model, dat = rf_pred_df, num.threads = 20)

pred_geno <- prediction$predictions
pred_geno <- matrix(pred_geno, nrow = nrow(conv_geno))

save(pred_geno,file = "./rda/beta_pred_geno.rda", compress = TRUE)
q(save = "no")

#shell script: echo R CMD BATCH comp.R | qsub -N comp_beta -cwd -l mf=40G,h_vmem=40G
