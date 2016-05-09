## contains code for the Mean Absolute Percent Score

map_score <- function(true, prediction){
  true_1 <- true[-1]
  true_n <- true[-length(true)]
  
  return(mean(abs(true-prediction))/((1/(length(true)-1))*sum(abs(true_1-true_n))))
}

RMSPE <- function(y, prediction){
  return(norm(as.matrix(y-prediction), type="F")/norm(as.matrix(y), type="F"))
}

PRMSPE <- function(y, prediction, lcl, ucl, alpha){
  r <- RMSPE(y, prediction)
  s <- ucl - lcl + (2/alpha)*(lcl - y)*(y < lcl) + (2/alpha)*(y-ucl)*(y > ucl)
  s_rel <- mean(s)/norm(as.matrix(y), type="F")
  
  return(r + s_rel)
}

pmap_score <- function(y, prediction, lcl, ucl, alpha){
    r <- map_score(y, prediction)
    s <- ucl - lcl + (2/alpha)*(lcl - y)*(y < lcl) + (2/alpha)*(y-ucl)*(y > ucl)
    s_rel <- mean(s)/norm(as.matrix(y), type="F")
    
    return(r + s_rel)
}

obs_sub <- read.csv("obs_sub.csv")
ind <- sample(2, nrow(obs_sub), replace = TRUE, prob=c(0.8, 0.2))
obs_sub$genotype <- factor(obs_sub$genotype)
obs_sub$plotid <- factor(obs_sub$plotid)


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
  
  upr1 <- p1$fit + (2 * p1$se.fit)
  lwr1 <- p1$fit - (2 * p1$se.fit)
  
  upr2 <- p2$fit + (2 * p2$se.fit)
  lwr2 <- p2$fit - (2 * p2$se.fit)
  
  preds <- data.frame(p1$fit, lwr1, upr1, p2$fit, lwr2, upr2)
  
  results[k,1:2] <- c(map_score(test[,'Stem'], p1$fit), map_score(test[,'Stem'], p2$fit))
  results[k,3:4] <- c(pmap_score(y = test[,'Stem'], prediction = p1$fit, 
                             lcl = lwr1, ucl = upr1, alpha = 0.1),
                      pmap_score(y = test[,'Stem'], prediction = p2$fit, 
                             lcl = lwr2, ucl = upr2, alpha = 0.1))
}
colnames(results) <- c("map 1", "map 2", "p map 1", "p map 2")
library(reshape2)
res_m <- melt(results)
colnames(res_m) <- c("iteration","measure", "value")
library(ggplot2)
res_m <- data.frame(res_m)
ggplot(res_m, aes(x=iteration, y=value)) + geom_line(aes(colour=measure))



#### ----- bootstrap MAP score

obs_sub <- read.csv("obs_sub.csv")
ind <- sample(2, nrow(obs_sub), replace = TRUE, prob=c(0.8, 0.2))
obs_sub$genotype <- factor(obs_sub$genotype)
obs_sub$plotid <- factor(obs_sub$plotid)
