
library(randomForest)
obs_sub <- observations[sample(nrow(observations), size = 10000),]

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
ggplot(x=results$rmspe, y=results$noise_sd)