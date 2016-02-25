load("/samples.Rdata")

phenotypes <- cbind(genotype = 9914+1:500, 
                    ensemble.samples$Miscanthus_x_giganteus)

knitr::kable(summary(phenotypes))
