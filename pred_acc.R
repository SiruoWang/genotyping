#library(ranger)
library(caret)

load("./rda/conv_geno.rda")
load("./rda/null_pred_geno.rda")

prior_id = c(1:ceiling(ncol(conv_geno)/3))
model_id = c((prior_id[length(prior_id)] + 1) : (prior_id[length(prior_id)] +ceiling(ncol(conv_geno)/3) ))
val_id = c((model_id[length(model_id)] + 1) : ncol(conv_geno))

known_geno = conv_geno[,val_id]

known_geno[which(known_geno == "CA")] = "AC"
known_geno[which(known_geno == "GA")] = "AG"
known_geno[which(known_geno == "TA")] = "AT"
known_geno[which(known_geno == "GC")] = "CG"
known_geno[which(known_geno == "TC")] = "CT"
known_geno[which(known_geno == "GT")] = "TG"

pred_geno[which(pred_geno == "CA")] = "AC"
pred_geno[which(pred_geno == "GA")] = "AG"
pred_geno[which(pred_geno == "TA")] = "AT"
pred_geno[which(pred_geno == "GC")] = "CG"
pred_geno[which(pred_geno == "TC")] = "CT"
pred_geno[which(pred_geno == "GT")] = "TG"

geno_level = c("AA","AC","AG","AT","CC","CG","CT","TT","TG","GG")
## check overall prediction accuracy per sample
accuracy_list_null = c()
for (i in c(1:ncol(known_geno))){

  sample_accuracy = confusionMatrix(factor(pred_geno[,i], levels = geno_level),factor(known_geno[,i],levels = geno_level))$overall['Accuracy']
  accuracy_list_null = c(accuracy_list_null, sample_accuracy)
}

confusionMatrix(factor(c(pred_geno), levels = geno_level),factor(c(known_geno),levels = geno_level))

###################################################
load("./rda/znb_pred_geno.rda")


pred_geno[which(pred_geno == "CA")] = "AC"
pred_geno[which(pred_geno == "GA")] = "AG"
pred_geno[which(pred_geno == "TA")] = "AT"
pred_geno[which(pred_geno == "GC")] = "CG"
pred_geno[which(pred_geno == "TC")] = "CT"
pred_geno[which(pred_geno == "GT")] = "TG"

## check overall prediction accuracy per sample
accuracy_list_znb = c()
for (i in c(1:ncol(known_geno))){

  sample_accuracy = confusionMatrix(factor(pred_geno[,i], levels = geno_level),factor(known_geno[,i],levels = geno_level))$overall['Accuracy']
  accuracy_list_znb = c(accuracy_list_znb, sample_accuracy)
}

confusionMatrix(factor(c(pred_geno), levels = geno_level),factor(c(known_geno),levels = geno_level))

#######################################################
load("./rda/beta_pred_geno.rda")

pred_geno[which(pred_geno == "CA")] = "AC"
pred_geno[which(pred_geno == "GA")] = "AG"
pred_geno[which(pred_geno == "TA")] = "AT"
pred_geno[which(pred_geno == "GC")] = "CG"
pred_geno[which(pred_geno == "TC")] = "CT"
pred_geno[which(pred_geno == "GT")] = "TG"

## check overall prediction accuracy per sample
accuracy_list_beta = c()
for (i in c(1:ncol(known_geno))){

  sample_accuracy = confusionMatrix(factor(pred_geno[,i], levels = geno_level),factor(known_geno[,i],levels = geno_level))$overall['Accuracy']
  accuracy_list_beta = c(accuracy_list_beta, sample_accuracy)
}

confusionMatrix(factor(c(pred_geno), levels = geno_level),factor(c(known_geno),levels = geno_level))

boxplot(accuracy_list_null,accuracy_list_znb,accuracy_list_beta)

#### pickrell data
#############################################################
load("./rda/pickrell_geno.rda") ## known_geno
load("./rda/pickrell_null_pred_geno.rda") ## pred_geno

pred_geno <- t(matrix(unlist(pred_geno), ncol = length(pred_geno[[1]]), byrow = TRUE))

known_geno[which(known_geno == "CA")] = "AC"
known_geno[which(known_geno == "GA")] = "AG"
known_geno[which(known_geno == "TA")] = "AT"
known_geno[which(known_geno == "GC")] = "CG"
known_geno[which(known_geno == "TC")] = "CT"
known_geno[which(known_geno == "GT")] = "TG"

pred_geno[which(pred_geno == "CA")] = "AC"
pred_geno[which(pred_geno == "GA")] = "AG"
pred_geno[which(pred_geno == "TA")] = "AT"
pred_geno[which(pred_geno == "GC")] = "CG"
pred_geno[which(pred_geno == "TC")] = "CT"
pred_geno[which(pred_geno == "GT")] = "TG"

geno_level = c("AA","AC","AG","AT","CC","CG","CT","TT","TG","GG")

pic_accuracy_list_null = c()
for (i in c(1:ncol(known_geno))){

  sample_accuracy = confusionMatrix(factor(pred_geno[,i], levels = geno_level),factor(known_geno[,i],levels = geno_level))$overall['Accuracy']
  pic_accuracy_list_null = c(pic_accuracy_list_null, sample_accuracy)
}

confusionMatrix(factor(c(pred_geno), levels = geno_level),factor(c(known_geno),levels = geno_level))

pickrell_seq_depth <- c(10.0, 10.3, 10.3, 9.4, 8.5, 8.3, 9.1, 10.2, 10.5, 10.6,
                        9.7, 5.8, 7.7, 8.1, 8.5, 8.5, 8.4, 7.1, 9.8, 10.1,
                        10.6, 10.6, 10.7, 10.1, 7.5, 4.6, 8.4, 7.6, 4.6, 8.2,
                        7.6, 8.4, 4.2, 4.0, 3.8, 6.2, 9.4, 8.9, 7.0, 7.1,
                        8.7, 8.8, 8.8, 9.0, 8.9, 9.0, 9.0, 9.2, 8.7, 9.2,
                        9.0, 8.7, 8.3, 8.2, 8.2, 10.0, 10.0, 10.0, 9.7, 7.1,
                        8.3, 8.5, 8.7, 8.6, 7.9, 8.4, 9.2, 8.7, 6.4, 6.2,
                        5.9, 6.8, 6.5, 10.3, 10.1, 9.1, 9.5, 11.7, 11.7, 11.6,
                        8.5, 7.7, 7.7, 6.8, 5.2, 6.0, 6.4, 9.3, 9.1, 7.7,
                        8.3, 8.8, 6.1, 6.6, 8.0, 7.2, 7.0, 8.0, 4.0, 4.0,
                        2.3, 2.0, 2.1)

pickrell_pred_df <- data.frame(pred_accuracy = pic_accuracy_list_null, seq_depth = pickrell_seq_depth)

model.lo <- loess(pred_accuracy ~ seq_depth, pickrell_pred_df, span=1.3)
max_depth <- max(pickrell_pred_df$seq_depth)
min_depth <- min(pickrell_pred_df$seq_depth)
prediction_curve <- predict(model.lo, data.frame(seq_depth = seq(min_depth, max_depth, 0.005)), se = TRUE)
fit_df = data.frame(seq_depth = seq(min_depth, max_depth, 0.005), pred_accuracy = prediction_curve$fit)
p = ggplot(pickrell_pred_df, aes(seq_depth, pred_accuracy)) +
    geom_point() +
    geom_point(data = fit_df, color = "blue", size = 0.1) +
    labs(x = "sequencing depth", y = "prediction accuracy", title = "global model" )
p


###############################################################
load("./rda/pickrell_beta_pred_geno.rda") ## pred_geno

pred_geno <- t(matrix(unlist(pred_geno), ncol = length(pred_geno[[1]]), byrow = TRUE))

pred_geno[which(pred_geno == "CA")] = "AC"
pred_geno[which(pred_geno == "GA")] = "AG"
pred_geno[which(pred_geno == "TA")] = "AT"
pred_geno[which(pred_geno == "GC")] = "CG"
pred_geno[which(pred_geno == "TC")] = "CT"
pred_geno[which(pred_geno == "GT")] = "TG"

confusionMatrix(factor(c(pred_geno), levels = geno_level),factor(c(known_geno),levels = geno_level))

###############################################################
load("./rda/pickrell_znb_pred_geno.rda") ## pred_geno

pred_geno <- t(matrix(unlist(pred_geno), ncol = length(pred_geno[[1]]), byrow = TRUE))

pred_geno[which(pred_geno == "CA")] = "AC"
pred_geno[which(pred_geno == "GA")] = "AG"
pred_geno[which(pred_geno == "TA")] = "AT"
pred_geno[which(pred_geno == "GC")] = "CG"
pred_geno[which(pred_geno == "TC")] = "CT"
pred_geno[which(pred_geno == "GT")] = "TG"

confusionMatrix(factor(c(pred_geno), levels = geno_level),factor(c(known_geno),levels = geno_level))
