library(ranger)
library(VGAM)

load("./rda/score_A.rda")
load("./rda/score_C.rda")
load("./rda/score_G.rda")
load("./rda/score_T.rda")
load("./rda/conv_geno.rda")

prior_id = c(1:ceiling(ncol(conv_geno)/3))
model_id = c((prior_id[length(prior_id)] + 1) : (prior_id[length(prior_id)] +ceiling(ncol(conv_geno)/3) ))
val_id = c((model_id[length(model_id)] + 1) : ncol(conv_geno))

####### fit beta distribution to data, and calculate priors by maximizing the log likelihood of beta distribution

beta_prior <- function(counts, total, alpha0, beta0){

  ## sum of negative log likelihood
  ll <- function(alpha, beta) {
    -sum(VGAM::dbetabinom.ab(counts, total, alpha, beta, log = TRUE))
  }

  ## calculate minimum of log likelihood
  m <- mle(ll, start = list(alpha = alpha0, beta = beta0), method = "L-BFGS-B", lower = c(0.00001, 0.01))

  ## get prior of beta distribution
  ab <- coef(m)

  return(ab)
}

### calculate the posterior ratio estimates using alpha.prior and beta.prior
ratio_est <- function(counts, total, alpha, beta){
  return((counts + alpha) / (total + alpha+ beta))
}

row_id = sample(c(1:nrow(conv_geno)),10000)
total = score_A + score_C + score_G + score_T

## ratio A
A_ab = beta_prior(c(score_A[row_id,prior_id]),c(total[row_id,prior_id]),2,2)
print(A_ab)
alpha.prior <- A_ab[1]
beta.prior <- A_ab[2]
A_est_modelset <- ratio_est(c(score_A[,model_id]), c(total[,model_id]), alpha.prior, beta.prior)
A_est_valset <- ratio_est(c(score_A[,val_id]), c(total[,val_id]), alpha.prior, beta.prior)
gc()

## ratio C
C_ab = beta_prior(c(score_C[row_id,prior_id]),c(total[row_id,prior_id]),2,2)
print(C_ab)
alpha.prior <- C_ab[1]
beta.prior <- C_ab[2]
C_est_modelset <- ratio_est(c(score_C[,model_id]), c(total[,model_id]), alpha.prior, beta.prior)
C_est_valset <- ratio_est(c(score_C[,val_id]), c(total[,val_id]), alpha.prior, beta.prior)
gc()

## ratio G
G_ab = beta_prior(c(score_G[row_id,prior_id]),c(total[row_id,prior_id]),2,2)
print(G_ab)
alpha.prior <- G_ab[1]
beta.prior <- G_ab[2]
G_est_modelset <- ratio_est(c(score_G[,model_id]), c(total[,model_id]), alpha.prior, beta.prior)
G_est_valset <- ratio_est(c(score_G[,val_id]), c(total[,val_id]), alpha.prior, beta.prior)
gc()

## ratio T
T_ab = beta_prior(c(score_T[row_id,prior_id]),c(total[row_id,prior_id]),2,2)
print(T_ab)
alpha.prior <- T_ab[1]
beta.prior <- T_ab[2]
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

#shell script: echo R CMD BATCH beta_model.R | qsub -N beta_model -cwd -l mf=60G,h_vmem=60G
