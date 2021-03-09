#Clear
#...variables
rm(list=ls())
#...console
cat("\014\n")
#...graphs
dev.off()
dev.new()

library(tidyverse)
library(edgeR)
library(MASS)

################################################################################
#Run the setup script
source("Scripts/setup1.R")

################################################################################

#Make factor variables
batch_FCT <- factor(c(1,1,1,2,2,2))

#...levels: baseline should go first
treat_FCT <- factor(
  c("treat","cont","mock","treat","cont","mock"),
  levels = c("cont","treat","mock")
)

tibble(sample = colnames(df_c), batch_FCT, treat_FCT)

#Make design matrix
design1 <- model.matrix(~batch_FCT + treat_FCT)
rownames(design) <- colnames(df_c)
design1

################################################################################
#Before attempting to loop through every gene, 
#...try the first two genes

#check offset

#glm:poisson and glm.nb require integer outcome values

#Poisson with batch as a covariate - check first 2 genes 
#test1b: also check if offset option can be written 
#...as part of the formula instead

test1 <- glm(
  formula = df_c[1,] ~ batch_FCT + treat_FCT,
  family = "poisson",
  offset = log(s_a*colSums(df_c))
)
test1
coef(summary(test1))

test1b <- glm(
  formula = df_c[1,] ~ batch_FCT + treat_FCT + offset(log(s_a*colSums(df_c))),
  family = "poisson"
)
test1b
coef(summary(test1b))


test2 <- glm(
  formula = df_c[2,] ~ batch_FCT + treat_FCT,
  family = "poisson",
  offset = log(s_a*colSums(df_c))
)
test2
coef(summary(test2))

#Poisson without batch as a covariate - check first 2 genes
test3 <- glm(
  formula = df_c[1,] ~ treat_FCT,
  family = "poisson",
  offset = log(s_a*colSums(df_c))
)
test3
coef(summary(test3))

test4 <- glm(
  formula = df_c[2,] ~ treat_FCT,
  family = "poisson",
  offset = log(s_a*colSums(df_c))
)
test4
coef(summary(test4))

#Neg binom with batch as a covariate - check first gene

test5 <- glm.nb(
  formula = df_c[1,] ~ batch_FCT + treat_FCT + offset(log(s_a*colSums(df_c)))
)
test5
coef(summary(test5))

#Neg binom without batch as a covariate - check first gene

test6 <- glm.nb(
  formula = df_c[1,] ~ treat_FCT + offset(log(s_a*colSums(df_c)))
)
test6
coef(summary(test6))

################################################################################
#Loop through each gene and perform GLMs

#First make template vectors of zeros 
#Loop to go through all the rows (genes):
#...for aic score (lower indicates better fit) and treatment/ batch p-values
#...extract the values and put them in the vectors


#Poisson with batch coefficient
aic_po1_v <- rep(0, dim(df_c)[1])
p_tr_po1_v <- rep(0, dim(df_c)[1]) 
p_ba_po1_v <- rep(0, dim(df_c)[1])

for (i in 1:dim(df_c)[1]) {
  z <- glm(
    formula = df_c[i,] ~ batch_FCT + treat_FCT,
    family = "poisson",
    offset = log(s_a*colSums(df_c))
  )
  
  aic_po1_v[i] <- z[["aic"]]
  p_tr_po1_v[i] <- coef(summary(z))[[3,4]]
  p_ba_po1_v[i] <- coef(summary(z))[[2,4]]
}
rm(z)


#Poisson without batch coefficient
aic_po2_v <- rep(0, dim(df_c)[1])
p_tr_po2_v <- rep(0, dim(df_c)[1]) 

for (i in 1:dim(df_c)[1]) {
  z <- glm(
    formula = df_c[i,] ~ treat_FCT,
    family = "poisson",
    offset = log(s_a*colSums(df_c))
  )
  
  aic_po2_v[i] <- z[["aic"]]
  p_tr_po2_v[i] <- coef(summary(z))[[2,4]]
}
rm(z)


#Negative binomial with batch coefficient
aic_nb1_v <- rep(0, dim(df_c)[1])
p_tr_nb1_v <- rep(0, dim(df_c)[1]) 
p_ba_nb1_v <- rep(0, dim(df_c)[1])

for (i in 1:dim(df_c)[1]) {
  z <- glm.nb(
    formula = df_c[i,] ~ batch_FCT + treat_FCT + offset(log(s_a*colSums(df_c)))
  )
  
  aic_nb1_v[i] <- z[["aic"]]
  p_tr_nb1_v[i] <- coef(summary(z))[[3,4]]
  p_ba_nb1_v[i] <- coef(summary(z))[[2,4]]
}
rm(z)
  

#Negative binomial without batch coefficient
aic_nb2_v <- rep(0, dim(df_c)[1])
p_tr_nb2_v <- rep(0, dim(df_c)[1]) 

for (i in 1:dim(df_c)[1]) {
  z <- glm.nb(
    formula = df_c[i,] ~ treat_FCT + offset(log(s_a*colSums(df_c)))
  )
  
  aic_nb2_v[i] <- z[["aic"]]
  p_tr_nb2_v[i] <- coef(summary(z))[[2,4]]
}
rm(z)


#plot of the both sets of aic densities (pois vs neg-binom)

#plots of p values for batch / treat coefficients (pois vs neg-binom)
