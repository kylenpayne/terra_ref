
library(randomForest)
observations <- read.csv("~/Desktop/observations.csv")
obs_sub <- observations[sample(nrow(observations), size = 100000),]
write.csv("obs_sub.csv", x=obs_sub)

obs_sub <- read.csv("obs_sub.csv")
ind <- sample(2, nrow(obs_sub), replace = TRUE, prob=c(0.8, 0.2))
obs_sub$genotype <- factor(obs_sub$genotype)
obs_sub$plotid <- factor(obs_sub$plotid)

RMSPE <- function(y, prediction){
  return(norm(as.matrix(y-prediction), type="F")/norm(as.matrix(y), type="F"))
}

obs_rf <- randomForest(Stem ~ plotid + precip, 
                       data=obs_sub[ind == 1,])

rmspe <- numeric(100)
for(k in 1:100){
  test <- obs_sub[ind==2,] 
  test['Stem'] <- test['Stem'] + rnorm(n=nrow(obs_sub[ind==2,]),
                                        mean = 0, sd=k/10)
  obs_pd <- predict(obs_rf, test)
  rmspe[k] <- RMSPE(test[,'Stem'], obs_pd)
}

library(ggplot2)
results <- data.frame(rmspe, noise_sd =(1:100)/10)
ggplot(results, aes(x=noise_sd, rmspe))+geom_line()+theme_bw()+
  xlab("Standard Deviation of Gaussian Noise") + 
  ylab("Relative Mean Squared Prediction Error") 
### penalized relative mean squared prediction error
PRMSPE <- function(y, prediction, lcl, ucl){
  return((norm(as.matrix(y-prediction), type="F")
          +sum(ucl - lcl))/norm(as.matrix(y), type="F"))
}

model_v_obs <- read.csv("~/Desktop/terra_ref/model_v_obs.csv")
library(reshape2)
model_v_obs_tmp <- split(model_v_obs, f= model_v_obs$type)
total <- merge(model_v_obs_tmp$model,model_v_obs_tmp$observed,by="genotype")


## penalized relative mean squared prediction error.

PRMSPE <- function(y, prediction, lcl, ucl, alpha){
  r <- RMSPE(y, prediction)
  s <- ucl - lcl + (2/alpha)*(lcl - y)*(y < lcl) + (2/alpha)*(y-ucl)*(y > ucl)
  s_rel <- mean(s)/norm(as.matrix(y), type="F")
  
  return(r + s_rel)
}


library(mgcv)
niter<-10
results <- matrix(NA, nrow=niter, ncol=4)
for(k in 1:niter){
  train <- obs_sub[ind==1,]
  train_sub <- train[sample(1:nrow(train), size=10000,replace = T),]
  gam1 <- gam(Stem ~ plotid + precip + Leaf + Root + genotype,
                         data=train_sub)

  gam2 <- gam(Stem ~ Rhizome + LAI + NDVI + precip, 
                         data=train_sub)
  
  test <- obs_sub[ind==2,] 
  p1 <- predict(gam1, newdata=test, type = "link", se.fit = TRUE)
  p2 <- predict(gam2, newdata=test, type = "link", se.fit = TRUE)
  
  upr1 <- p1$fit + (pnorm(.975) * p1$se.fit)
  lwr1 <- p1$fit - (pnorm(.975) * p1$se.fit)
  
  upr2 <- p2$fit + (pnorm(.975) * p2$se.fit)
  lwr2 <- p2$fit - (pnorm(.975) * p2$se.fit)
  
  upr3 <- p2$fit + (pnorm(.95) * p2$se.fit)
  lwr3 <- p2$fit - (pnorm(.95) * p2$se.fit)
  
  upr4 <- p2$fit + (pnorm(.995) * p2$se.fit)
  lwr4 <- p2$fit - (pnorm(.995) * p2$se.fit)
  
  results[k,] <- c(
                  PRMSPE(y = test[,'Stem'], prediction = p1$fit, ## alpha = 0.05
                           lcl = lwr1, ucl = upr1, alpha = 0.05),
                  PRMSPE(y = test[,'Stem'], prediction = p2$fit, 
                         lcl = lwr2, ucl = upr2, alpha = 0.05),
                  PRMSPE(y=test[,'Stem'], prediction = p2$fit, ## alpha = 0.10
                         lcl=lwr3, ucl=upr3, alpha=0.10),
                  PRMSPE(y=test[,'Stem'], prediction = p2$fit, ## alpha = 0.10
                         lcl=lwr4, ucl=upr4, alpha=0.01)
                  )
}
colnames(results) <- c("PRMSPE1", "PRMSPE2", "PRMSPE.10", "PRMSPE.01")
library(reshape2)
res_m <- melt(results)
colnames(res_m) <- c("iteration","measure", "value")
library(ggplot2)
res_m <- data.frame(res_m)
ggplot(res_m, aes(x=iteration, y=value)) + geom_line(aes(colour=measure))



### fitting a large lasso regression of genotype

library(glmnet)

X <- data.frame(genotype=obs_sub$genotype, plotid=obs_sub$plotid)
X <- dummy(X)
X <- data.matrix(X)
Y <- data.matrix(obs_sub$LAI)
lasso_mod <- glmnet(x=X, y=Y)

## set the seed to make your partition reproductible
set.seed(123)

## 75% of the sample size
smp_size <- floor(0.75 * nrow(obs_sub))

train_ind <- sample(seq_len(nrow(obs_sub)), size = smp_size)


train_x <- X[train_ind, ]
test_x <- X[-train_ind, ]
train_y <- Y[train_ind, ]
test_y <- Y[-train_ind, ]

X_0 <- train_x[,]
X_1 <- train_x[, c(157, 213, 92)]


test_x_0 <- test_x[,c(20,25,67)]
test_x_1 <- test_x[,c(157, 213, 92)]

modl_0 <- glmnet(x=X_0, y=train_y)
modl_1 <- glmnet(x=X_1, y=train_y)

preds_0 <- predict(modl_0, type="response", newx=test_x_0, s = 150)
preds_1 <- predict(modl_1, newx=test_x_1, s = 150)



PRMSPE <- function(y, prediction, lcl, ucl, alpha){
  r <- RMSPE(y, prediction)
  s <- ucl - lcl + (2/alpha)*(lcl - y)*(y < lcl) + (2/alpha)*(y-ucl)*(y > ucl)
  s_rel <- mean(s)/norm(as.matrix(y), type="F")
  
  return(r + s_rel)
}

resamps_0 <- lapply(1:20, function(i)
  sample(preds_0, replace = T))
resamps_1 <- lapply(1:20, function(i)
  sample(preds_1, replace = T))

r_sd_0 <- sapply(resamps_0, sd)
r_sd_1 <- sapply(resamps_1, sd)




library(mgcv)
niter<-10
results <- matrix(NA, nrow=niter, ncol=4)
for(k in 1:niter){
  train <- obs_sub[ind==1,]
  train_sub <- train[sample(1:nrow(train), size=10000,replace = T),]
  gam1 <- gam(Stem ~ plotid + precip + Leaf + Root + genotype,
              data=train_sub)
  
  gam2 <- gam(Stem ~ Rhizome + LAI + NDVI + precip, 
              data=train_sub)
  
  test <- obs_sub[ind==2,] 
  p1 <- predict(gam1, newdata=test, type = "link", se.fit = TRUE)
  p2 <- predict(gam2, newdata=test, type = "link", se.fit = TRUE)
  
  upr1 <- p1$fit + (pnorm(.975) * p1$se.fit)
  lwr1 <- p1$fit - (pnorm(.975) * p1$se.fit)
  
  upr2 <- p2$fit + (pnorm(.975) * p2$se.fit)
  lwr2 <- p2$fit - (pnorm(.975) * p2$se.fit)
  
  upr3 <- p2$fit + (pnorm(.95) * p2$se.fit)
  lwr3 <- p2$fit - (pnorm(.95) * p2$se.fit)
  
  upr4 <- p2$fit + (pnorm(.995) * p2$se.fit)
  lwr4 <- p2$fit - (pnorm(.995) * p2$se.fit)
  
  results[k,] <- c(
    PRMSPE(y = test[,'Stem'], prediction = p1$fit, ## alpha = 0.05
           lcl = lwr1, ucl = upr1, alpha = 0.05),
    PRMSPE(y = test[,'Stem'], prediction = p2$fit, 
           lcl = lwr2, ucl = upr2, alpha = 0.05),
    PRMSPE(y=test[,'Stem'], prediction = p2$fit, ## alpha = 0.10
           lcl=lwr3, ucl=upr3, alpha=0.10),
    PRMSPE(y=test[,'Stem'], prediction = p2$fit, ## alpha = 0.10
           lcl=lwr4, ucl=upr4, alpha=0.01)
  )
}
colnames(results) <- c("PRMSPE1", "PRMSPE2", "PRMSPE.10", "PRMSPE.01")
library(reshape2)
res_m <- melt(results)
colnames(res_m) <- c("iteration","measure", "value")
library(ggplot2)
res_m <- data.frame(res_m)
ggplot(res_m, aes(x=iteration, y=value)) + geom_line(aes(color=measure,size=2))




ggplot(res_m, aes(x=iteration, y=value)) + geom_line(aes(color=measure,size=2)) + 
  annotate("text", x = 2.5, y = 0.65, label = c("Model 2 with 99% Interval"), parse=FALSE) + 
  annotate("text", x = 2.5, y = 0.475, label = c("Model 2 with 90% Itner"), parse=FALSE) + 
  geom_segment(aes(x = 2.5, y = .67, xend = 2.5, yend = .71), arrow = arrow(length =unit(0.5, "cm"))) + 
  geom_segment(aes(x = 2.5, y = .45, xend = 2.5, yend = .375), arrow = arrow(length = unit(0.5, "cm"))) + 
  theme_bw()

n_alpha <- 50

alpha <- seq(0.001, 0.5, length.out = n_alpha)

results <- matrix(NA, nrow=n_alpha, ncol=2)
for(k in 1:n_alpha){
  train <- obs_sub[ind==1,]
  train_sub <- train[sample(1:nrow(train), size=10000,replace = T),]
  gam1 <- gam(Stem ~ plotid + precip + Leaf + Root + genotype,
              data=train_sub)
  
  gam2 <- gam(Stem ~ Rhizome + LAI + NDVI + precip, 
              data=train_sub)
  
  test <- obs_sub[ind==2,] 
  p1 <- predict(gam1, newdata=test, type = "link", se.fit = TRUE)
  p2 <- predict(gam2, newdata=test, type = "link", se.fit = TRUE)
  
  upr1 <- p1$fit + (pnorm(1-alpha[k]/2) * p1$se.fit)
  lwr1 <- p1$fit - (pnorm(1-alpha[k]/2) * p1$se.fit)
  
  upr2 <- p2$fit + (pnorm(1-alpha[k]/2) * p2$se.fit)
  lwr2 <- p2$fit - (pnorm(1-alpha[k]/2) * p2$se.fit)
  
  
  results[k,] <- c(
    PRMSPE(y = test[,'Stem'], prediction = p1$fit, ## alpha = 0.05
           lcl = lwr1, ucl = upr1, alpha = 1-alpha[k]/2),
    PRMSPE(y = test[,'Stem'], prediction = p2$fit, 
           lcl = lwr2, ucl = upr2, alpha = 1-alpha[k]/2)
  )
}



