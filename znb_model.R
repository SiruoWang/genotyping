library(ranger)
load("./rda/score_A.rda")
load("./rda/score_C.rda")
load("./rda/score_G.rda")
load("./rda/score_T.rda")
load("./rda/conv_geno.rda")

prior_id = c(1:ceiling(ncol(conv_geno)/3))
model_id = c((prior_id[length(prior_id)] + 1) : (prior_id[length(prior_id)] +ceiling(ncol(conv_geno)/3) ))
val_id = c((model_id[length(model_id)] + 1) : ncol(conv_geno))

####### fit zero inflated negative binomial distribution to 1/3 of data, and calculate priors using em


### znb function fits the zero inflated negative binomial distribution and calculate priors using em
znb <- function(nuc_counts){
  ## filter out outliers -- counts over 500 (1%)
  y <- nuc_counts[which(nuc_counts < 500)]


  ## set initial value for parameters -- mu, phi, p
  mu_int <- 20
  phi_int <- 5
  p_int <- 0.5

  ## sum of log likelihood function of zero inflated negative binomial
  ll <- function(param_est,y,I_0){
    mu <- param_est[1]
    phi <- param_est[2]
    p <- param_est[3]
    #f_nb <- (gamma(y + 1/phi) / (gamma(y+1) * gamma(1/phi))) * ((1/phi) / (1/phi + mu))^(1/phi) * (mu/(mu + 1/phi))^y
    return(sum(I_0*log(p) + (1-I_0)*(log(1-p)+lgamma(y+1/phi)-lgamma(y+1)-lgamma(1/phi)+(1/phi)*log((1/phi)/(1/phi + mu))+y*log(mu/(mu + 1/phi)) )))
  }

  ## E-step function in the EM algorithm. I(0) is the only random variable, so we need to calculate expection of Indicator function
  Expect_I <- function(y,mu,phi,p){
    f_nb <-  ((1/phi) / (1/phi + mu))^(1/phi)
    output <- rep(p / (p + (1-p)*f_nb),length(y))
    output[y>0] <- 0
    return(output)
  }

  ## In the K iteration of M-step, parameters denote as mu_k, phi_k, p_k
  mu_k <- mu_int
  phi_k <- phi_int
  p_k <- p_int

  ## keep track of value of log likehood, mu, phi, p through each iteration in the M-step
  ll_record = mu_record = phi_record = p_record = c()
  ll_record[1] <- 0
  mu_record[1] <- 0
  phi_record[1] <- 0
  p_record[1] <- 0

  ## indicator function I(0)
  ## E-step: expectation of I(0). Replace the value of I(0) in the kth iteration by expectation
  I_0_k <- Expect_I(y,mu_k,phi_k,p_k)

  ## M-step: maximize log likehood function
  ## set constrains for mu, phi, p while maximazing the log likelihood. mu>0, phi>0, 0<=p<=1
  opt_result <- constrOptim(c(mu_k,phi_k,p_k), f = ll, grad = NULL, ui=rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(0,0,-1)), ci= c(0,0,0,-1), control = list(fnscale = -1), y=y, I_0=I_0_k)

  ## keep track of values of parameters -- mu, phi, p, values of ll, in the first iteration of M step
  mu_record[2] <- mu_k <- opt_result$par[1]
  phi_record[2] <- phi_k <- opt_result$par[2]
  p_record[2] <- p_k <- opt_result$par[3]
  ll_record[2] <- opt_result$value

  i = 2

  ## while the distence between log likehihood functions is small enough (less than 0.01), the iteration stops, and the values of priors are maximized;
  ## otherwise, stay in the iteration loop of EM algorithm.
  ## Note: the value of log likehood should be increading through iterations and become the max when the iteration stops
  while(abs(ll_record[i] - ll_record[i-1]) >= 0.01){

    I_0_k <- Expect_I(y,mu_k,phi_k,p_k)

    opt_result <- constrOptim(c(mu_k,phi_k,p_k), f = ll, grad = NULL, ui=rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(0,0,-1)), ci= c(0,0,0,-1), control = list(fnscale = -1), y=y, I_0=I_0_k)
    #opt_result <- constrOptim(c(mu_k,phi_k,p_k), f = ll, grad = NULL, ui=rbind(c(0,0,1),c(0,0,-1)), ci= c(0,-1), control = list(fnscale = -1), y=y, I_0=I_0_k)
    #opt_result <- optim(c(mu_k,phi_k,p_k), fn = ll, gr = NULL, y, I_0_k, control = list(fnscale = -1))
    mu_record[i] <- mu_k <- opt_result$par[1]
    phi_record[i] <- phi_k <- opt_result$par[2]
    p_record[i] <- p_k <- opt_result$par[3]
    i <- i + 1
    ll_record[i] <- opt_result$value

  } ## end while loop

  ## get the estimated priors from EM
  lamda = mu_record[length(mu_record)]
  phi = phi_record[length(phi_record)]
  p = p_record[length(p_record)]

  ## return priors
  return(c(lamda,phi,p))
}

### calculate posterior expection for estimated counts using estimated priors from EM
count_est <- function(y,lamda,phi,p){
  f_nb <-  ((1/phi) / (1/phi + lamda))^(1/phi)
  prob_I_0 <- p / (p + (1-p)*f_nb)
  return((y+1/phi)/(1/(lamda*phi)+1) * (1 - prob_I_0))
}

## subset 20% snps from all available to find priors in znb function
row_id = sample(c(1:nrow(conv_geno)),20000)

## score A
## calculate values of priors
A_est = znb(c(score_A[row_id,prior_id]))
print(A_est)
## calculate posterior estimates of count values in the model set
A_est_modelset = count_est(c(score_A[,model_id]),A_est[1],A_est[2],A_est[3])
## calculate posterior estimates of count values in the validation set
A_est_valset = count_est(c(score_A[,val_id]),A_est[1],A_est[2],A_est[3])
rm(score_A)
gc()


## score C
C_est = znb(c(score_C[row_id,prior_id]))
print(C_est)
C_est_modelset = count_est(c(score_C[,model_id]),C_est[1],C_est[2],C_est[3])
C_est_valset = count_est(c(score_C[,val_id]),C_est[1],C_est[2],C_est[3])
rm(score_C)
gc()

## score G
G_est = znb(c(score_G[row_id,prior_id]))
print(G_est)
G_est_modelset = count_est(c(score_G[,model_id]),G_est[1],G_est[2],G_est[3])
G_est_valset = count_est(c(score_G[,val_id]),G_est[1],G_est[2],G_est[3])
rm(score_G)
gc()

## score T
T_est = znb(c(score_T[row_id,prior_id]))
print(T_est)
T_est_modelset = count_est(c(score_T[,model_id]),T_est[1],T_est[2],T_est[3])
T_est_valset = count_est(c(score_T[,val_id]),T_est[1],T_est[2],T_est[3])
rm(score_T)
gc()

############################ build model in the model set and make predictions in the val set

### fit random forest model using ranger package
total = A_est_modelset + C_est_modelset + G_est_modelset + T_est_modelset

geno = c(conv_geno[,model_id])
ratioA = A_est_modelset / total
ratioC = C_est_modelset / total
ratioG = G_est_modelset / total
ratioT = T_est_modelset / total

ratioA[is.na(ratioA)] = 0
ratioC[is.na(ratioC)] = 0
ratioG[is.na(ratioG)] = 0
ratioT[is.na(ratioT)] = 0

rf_model_df = data.frame(geno = geno, ratioA = ratioA, ratioC = ratioC, ratioG = ratioG, ratioT = ratioT)
znb_model = ranger(geno ~ ratioA + ratioC + ratioG + ratioT, data = rf_model_df, num.trees = 100, write.forest = TRUE, num.threads = 20)

save(znb_model, file = "./rda/znb_model.rda", compress = TRUE)

rm(rf_model_df)
gc()


#### make predictions using znb_model on the validation set
total = A_est_valset + C_est_valset + G_est_valset + T_est_valset

ratioA = A_est_valset / total
ratioC = C_est_valset / total
ratioG = G_est_valset / total
ratioT = T_est_valset / total

ratioA[is.na(ratioA)] = 0
ratioC[is.na(ratioC)] = 0
ratioG[is.na(ratioG)] = 0
ratioT[is.na(ratioT)] = 0

rf_pred_df <- data.frame(ratioA = ratioA, ratioC = ratioC, ratioG = ratioG, ratioT = ratioT)
prediction <- predict(znb_model, dat = rf_pred_df, num.threads = 20)

pred_geno <- prediction$predictions
pred_geno <- matrix(pred_geno, nrow = nrow(conv_geno))

save(pred_geno,file = "./rda/znb_pred_geno.rda", compress = TRUE)
q(save = "no")

#shell script: echo R CMD BATCH znb_model.R | qsub -N znb_model -cwd -l mf=60G,h_vmem=60G
