## before correcting
load("~/Documents/2017_recount_genotype/geuvadis/genotype_call/rda/score_A.rda")
load("~/Documents/2017_recount_genotype/geuvadis/genotype_call/rda/score_C.rda")
load("~/Documents/2017_recount_genotype/geuvadis/genotype_call/rda/score_G.rda")
load("~/Documents/2017_recount_genotype/geuvadis/genotype_call/rda/score_T.rda")

total = score_A[,val_id] + score_C[,val_id] + score_G[,val_id] + score_T[,val_id]
ratio_A = score_A[,val_id] / total
ratio_A[is.na(ratio_A)] = 0
plot(density(c(ratio_A[,1:10])))


## after correcting
load("./rda/A_est_valset.rda")
load("./rda/C_est_valset.rda")
load("./rda/G_est_valset.rda")
load("./rda/T_est_valset.rda")

total_est = A_est_valset + C_est_valset + G_est_valset + T_est_valset
total_est = matrix(total_est, nrow = nrow(score_A))
score_Aest = matrix(A_est_valset, nrow = nrow(score_A))
ratio_Aest = score_Aest / total_est
ratio_Aest[is.na(ratio_Aest)] = 0
lines(density(c(ratio_Aest[,1:10])), col = "red")


###################

# ratioA = score_A[,val_id[1:5]]/total[,val_id[1:5]]
# ratioA[is.na(ratioA)] = 0
# plot(density(c(ratioA)))
#
# alpha.prior
# beta.prior
#
# #lines(density(rbeta(3000,alpha.prior,beta.prior )), col = "red")
#
# est_ratioA = (c(score_A[,val_id[1:5]]) + alpha.prior) / (c(total[,val_id[1:5]]) + alpha.prior + beta.prior)
# lines(density(est_ratioA),col = "red")
