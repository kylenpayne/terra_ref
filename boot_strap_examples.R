obs_sub <- read.csv("obs_sub.csv")
ind <- sample(2, nrow(obs_sub), replace = TRUE, prob=c(0.8, 0.2))
obs_sub$genotype <- factor(obs_sub$genotype)
obs_sub$plotid <- factor(obs_sub$plotid)

RMSPE <- function(y, prediction){
  return(norm(as.matrix(y-prediction), type="F")/norm(as.matrix(y), type="F"))
}



MASE_MAS <- function(y, prediction, lcl, ucl, alpha){
  MASE <- mean(abs(y-prediction))
  
  divide <- 0
  for(k in 2:length(y)){
    divide <- divide+abs(y[k] - y[k-1])
  }
  divide <- (1/(length(y)-1))*divide
  
  MASE<-MASE*divide
  
  s <- ucl - lcl + (2/alpha)*(lcl - y)*(y < lcl) + (2/alpha)*(y-ucl)*(y > ucl)
  mean_s <- mean(s)
  
  d<-1/(length(s)-1)
  denom <- 0
  for(k in 2:length(s)){
    denom <- denom + abs(s[k] - s[k-1])
  }
  denom <- d*denom
  
  s_rel <- mean_s/denom
  
  return(MASE + s_rel)
}



library(mgcv)

train <- obs_sub[ind==1,]
train_sub <- train[sample(1:nrow(train), size=10000,replace = T),]
gam1 <- gam(Stem ~ plotid + precip + Leaf + Root + genotype,
            data=train_sub)
gam2 <- gam(Stem ~ Rhizome + LAI + NDVI + precip, 
            data=train_sub)

test <- obs_sub[ind==2,] 
p1 <- predict(gam1, newdata=test, type = "link", se.fit = TRUE)
p2 <- predict(gam2, newdata=test, type = "link", se.fit = TRUE)



### same model, two different prediction intervals. 
n_iter <- 1000
results <- matrix(NA, ncol=2, nrow=n_iter)

eps <- 0.01
for(k in 2:n_iter){
  ### normal approximation interval
  upr_norm <- p2$fit + (qnorm(.95) * (p2$se.fit))
  lwr_norm <- p2$fit - (qnorm(.95) * (p2$se.fit))
  
  boot_pred_means <- replicate(k, mean(sample(p2$fit, 100, replace=F)))
  boot_pred_means_std <- (boot_pred_means-mean(boot_pred_means))/(sd(boot_pred_means))
  
  piv_up <- quantile(boot_pred_means_std, 0.95)
  piv_down <- quantile(boot_pred_means_std, 0.05)
  
  
  upr_boot <- p2$fit + (piv_up*p2$se.fit)
  lwr_boot <- p2$fit + (piv_down*p2$se.fit) 
  
  
  p_norm <- MASE_MAS(y = test[,'Stem'], prediction = p2$fit, ## alpha = 0.05
                   lcl = lwr_norm, ucl = upr_norm, alpha = 0.1)
  
  p_boot<-MASE_MAS(y = test[,'Stem'], prediction = p2$fit, ## alpha = 0.05
                 lcl = lwr_boot, ucl = upr_boot, alpha = 0.1)
  
  results[k,] <- c(p_norm, p_boot)
}


results<-results[-1,]

colnames(results) <- c("norm", "boot")
results <- data.frame(results)
library(reshape2)

results_m <- melt(results)
results_m <- cbind(results_m, rep(2:n_iter, 2))



library(ggplot2)
colnames(results_m) <- c("int_type", "R", "boot_samps") 
ggplot(results_m,aes(x=boot_samps, y=R))+geom_line(aes(colour=int_type))



### ----- import dow data

variety_data<-read.csv('Variety_SC_K_Colon_Moisture_Yes_1A.csv')
variety_data$GenoID <- factor(variety_data$GenoID)

variety_sub <- variety_data[sample(1:nrow(variety_data), size=10000, replace = FALSE),]

ind <- sample(2, nrow(variety_sub), replace = TRUE, prob=c(0.8, 0.2))

variety_data_train <- variety_sub[ind==1,]
variety_data_test <- variety_sub[ind==2,]

library(dummies)

dummy_genoID <- variety_data_train$GenoID
genos <- dummy(dummy_genoID)

y <- as.matrix(variety_data_train$Yield_TPH)
library(glmnet)
library(mgcv)
gam_mod<-gam(Yield_TPH~Moisture, data=variety_data_train)