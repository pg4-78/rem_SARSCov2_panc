################################################################################
#Script-file: glm_separate3a.R
#Description: PART A of poisson vs negative binomial
################################################################################

#Re-run the setup script (TRUE)
#...alternatively, load a save(FALSE)
if (FALSE) {
  source("Scripts/setup2.R")
  #or setup1 for no average log CPM filter
} else {
  
  #Clear
  #...variables
  rm(list=ls())
  #...console
  cat("\014\n")
  #...graphs
  tryCatch(dev.off(), error = function(e) {NULL})
  dev.new()
  
  load(file = "./Data/setup2_co_2p01.RData")
  
  #Weak filter only
  #file = "./Data/setup2_wo.RData"
  
  #Weak filter and cutoff with opt_aveCPM_thresh
  #file = "./Data/setup2_wc_2pn3.RData"
  #file = "./Data/setup2_wc_2p01.RData"
  #file = "./Data/setup2_wc_2p03.RData"
  #file = "./Data/setup2_wc_2p05.RData"
  #average CPM cutoff at 2 to the power of (...)
  #-3 (below the lowest, no filter); 1, 3, 5
  
  #Cutoff with opt_aveCPM_thresh only; no weak filter
  #file = "./Data/setup2_co_2p01.RData"
  
}

library(data.table)
library(tidyverse)
library(edgeR)
library(MASS)


################################################################################

#Make factor variables
batch_FCT_ <- factor(c(0,0,0,1,1,1))

#...levels: baseline should go first: baseline is uninfected:untreated

infect_FCT_ <- factor(c(1,1,0,1,1,0))

treat_FCT_ <- factor(c(1,0,0,1,0,0))

tibble(sample = colnames(df_c), batch_FCT_, infect_FCT_, treat_FCT_)

#Check design matrices
design1 <- model.matrix(~batch_FCT_ + infect_FCT_ + treat_FCT_)
rownames(design1) <- colnames(df_c)
design1

design2 <- model.matrix(~infect_FCT_ + treat_FCT_)
rownames(design2) <- colnames(df_c)
design2

################################################################################
#Convenient notation for skipping error/warning regressions (eg zeros)

#Include warnings in  tryCatch?
#Skip errors only (FALSE)
#Alternatively, skip errors and warnings (TRUE)

if (TRUE) {
  catch_e <- function(x) {
    tryCatch(x,
      error = function(e) {NA},
      warning = function(w) {NA}
    )
  }
} else {
  catch_e <- function(x) {
    tryCatch(x,
      error = function(e) {NA}
    )
  }
}

################################################################################
#For ease of troubleshooting: Before attempting to loop through every gene, 
#...try one gene

#Choose a gene between 1 and (...18302) (inclusive)
dim(df_c)[1]

df_c_tb <- as_tibble(df_c) %>% 
  add_column("gene_name" = df_c_nm) %>% 
  add_column("row_num" = as.numeric(1:dim(df_c)[1]))

#Error genes: 22: "KLHL17"
gene_num <- 1

#check offset

#glm:poisson and glm.nb require integer outcome values

####################
#______Poisson with batch as a covariate - check gene_num 

test1 <- catch_e(
  glm(
    formula = df_c[gene_num,] ~ batch_FCT_ + infect_FCT_ + treat_FCT_ + offset(log(s_a*og_lib_sz_v)),
    family = "poisson"
  )
)
test1
catch_e(coef(summary(test1)))

####################
#______Poisson without batch as a covariate - check gene_num 

test2 <- catch_e(
  glm(
    formula = df_c[gene_num,] ~ infect_FCT_ + treat_FCT_ + offset(log(s_a*og_lib_sz_v)),
    family = "poisson"
  )
)
test2
catch_e(coef(summary(test2)))

####################
#______Neg binom with batch as a covariate - check gene_num 

test3 <- catch_e(
  glm.nb(formula = df_c[gene_num,] ~ batch_FCT_ + infect_FCT_ + treat_FCT_ + offset(log(s_a*og_lib_sz_v)))
)
test3
catch_e(coef(summary(test3)))

#Test extraction of key values
test3_aic <- catch_e(test3[["aic"]])
test3_B_batch <- catch_e(coef(summary(test3))[["batch_FCT_1",1]])
test3_B_infect <- catch_e(coef(summary(test3))[["infect_FCT_1",1]])
test3_B_treat <- catch_e(coef(summary(test3))[["treat_FCT_1",1]])
test3_p_batch <- catch_e(coef(summary(test3))[["batch_FCT_1",4]])
test3_p_infect <- catch_e(coef(summary(test3))[["infect_FCT_1",4]])
test3_p_treat <- catch_e(coef(summary(test3))[["treat_FCT_1",4]])

test3_aic
test3_B_batch
test3_B_infect
test3_B_treat
test3_p_batch
test3_p_infect
test3_p_treat

####################
#______Neg binom without batch as a covariate - check gene_num 

test4 <- catch_e(
  glm.nb(formula = df_c[gene_num,] ~ infect_FCT_ + treat_FCT_ + offset(log(s_a*og_lib_sz_v)))
)
test4
catch_e(coef(summary(test4)))

################################################################################
#Either loop through each gene and perform GLMs
#...or just load the saved vectors
#First make template vectors of NAs
#Loop to go through all the rows (genes):
#...for aic score (lower indicates better fit) and treatment/ batch p-values
#...extract the values and put them in the vectors

#Count number of zeros in each column
#Using as.data.table for more convenient notation
df_c_dt <- as.data.table(df_c)
dim(df_c_dt[treat0 == 0, ]) [1]
dim(df_c_dt[cont0 == 0, ]) [1]
dim(df_c_dt[mock0 == 0, ]) [1]
dim(df_c_dt[treat1 == 0, ]) [1]
dim(df_c_dt[cont1 == 0, ]) [1]
dim(df_c_dt[mock1 == 0, ]) [1]

#Possible option: Only use rows with non-zero values for all 6 counts
zero_filter <- FALSE
if (zero_filter == TRUE) {
  df_c <- df_c[treat0*cont0*mock0*treat1*cont1*mock1!=0,]
}

####################
#______Poisson with batch coefficient

po1_aic <- as.numeric(rep(NA, dim(df_c)[1]))
po1_B_b <- as.numeric(rep(NA, dim(df_c)[1]))
po1_B_i <- as.numeric(rep(NA, dim(df_c)[1]))
po1_B_t <- as.numeric(rep(NA, dim(df_c)[1]))
po1_p_b <- as.numeric(rep(NA, dim(df_c)[1]))
po1_p_i <- as.numeric(rep(NA, dim(df_c)[1]))
po1_p_t <- as.numeric(rep(NA, dim(df_c)[1]))

for (i in 1:dim(df_c)[1]) {
  z <- catch_e(
    glm(
      formula = df_c[i,] ~ batch_FCT_ + infect_FCT_ + treat_FCT_ + offset(log(s_a*og_lib_sz_v)),
      family = "poisson"
    )
  )
  po1_aic[i] <- catch_e(z[["aic"]])
  po1_B_b[i] <- catch_e(coef(summary(z))[["batch_FCT_1",1]])
  po1_B_i[i] <- catch_e(coef(summary(z))[["infect_FCT_1",1]])
  po1_B_t[i] <- catch_e(coef(summary(z))[["treat_FCT_1",1]])
  po1_p_b[i] <- catch_e(coef(summary(z))[["batch_FCT_1",4]])
  po1_p_i[i] <- catch_e(coef(summary(z))[["infect_FCT_1",4]])
  po1_p_t[i] <- catch_e(coef(summary(z))[["treat_FCT_1",4]])
}
rm(z)

####################
#______Poisson without batch coefficient
po2_aic <- as.numeric(rep(NA, dim(df_c)[1]))
po2_B_i <- as.numeric(rep(NA, dim(df_c)[1]))
po2_B_t <- as.numeric(rep(NA, dim(df_c)[1]))
po2_p_i <- as.numeric(rep(NA, dim(df_c)[1]))
po2_p_t <- as.numeric(rep(NA, dim(df_c)[1]))

for (i in 1:dim(df_c)[1]) {
  z <- catch_e(
    glm(
      formula = df_c[i,] ~ infect_FCT_ + treat_FCT_ + offset(log(s_a*og_lib_sz_v)),
      family = "poisson"
    )
  )
  po2_aic[i] <- catch_e(z[["aic"]])
  po2_B_i[i] <- catch_e(coef(summary(z))[["infect_FCT_1",1]])
  po2_B_t[i] <- catch_e(coef(summary(z))[["treat_FCT_1",1]])
  po2_p_i[i] <- catch_e(coef(summary(z))[["infect_FCT_1",4]])
  po2_p_t[i] <- catch_e(coef(summary(z))[["treat_FCT_1",4]])
}
rm(z)


####################
#______Negative binomial with batch coefficient
nb1_aic <- as.numeric(rep(NA, dim(df_c)[1]))
nb1_B_b <- as.numeric(rep(NA, dim(df_c)[1]))
nb1_B_i <- as.numeric(rep(NA, dim(df_c)[1]))
nb1_B_t <- as.numeric(rep(NA, dim(df_c)[1]))
nb1_p_b <- as.numeric(rep(NA, dim(df_c)[1]))
nb1_p_i <- as.numeric(rep(NA, dim(df_c)[1]))
nb1_p_t <- as.numeric(rep(NA, dim(df_c)[1]))

for (i in 1:dim(df_c)[1]) {
  z <- catch_e(
    glm.nb(formula = df_c[i,] ~ batch_FCT_ + infect_FCT_ + treat_FCT_ + offset(log(s_a*og_lib_sz_v)))
  )
  
  nb1_aic[i] <- catch_e(z[["aic"]])
  nb1_B_b[i] <- catch_e(coef(summary(z))[["batch_FCT_1",1]])
  nb1_B_i[i] <- catch_e(coef(summary(z))[["infect_FCT_1",1]])
  nb1_B_t[i] <- catch_e(coef(summary(z))[["treat_FCT_1",1]])
  nb1_p_b[i] <- catch_e(coef(summary(z))[["batch_FCT_1",4]])
  nb1_p_i[i] <- catch_e(coef(summary(z))[["infect_FCT_1",4]])
  nb1_p_t[i] <- catch_e(coef(summary(z))[["treat_FCT_1",4]])
  
}
rm(z)
  

####################
#______Negative binomial without batch coefficient
nb2_aic <- as.numeric(rep(NA, dim(df_c)[1]))
nb2_B_i <- as.numeric(rep(NA, dim(df_c)[1]))
nb2_B_t <- as.numeric(rep(NA, dim(df_c)[1]))
nb2_p_i <- as.numeric(rep(NA, dim(df_c)[1]))
nb2_p_t <- as.numeric(rep(NA, dim(df_c)[1]))

for (i in 1:dim(df_c)[1]) {
  z <- catch_e(
      glm.nb(formula = df_c[i,] ~ infect_FCT_ + treat_FCT_ + offset(log(s_a*og_lib_sz_v)))
  )
  
  nb2_aic[i] <- catch_e(z[["aic"]])
  nb2_B_i[i] <- catch_e(coef(summary(z))[["infect_FCT_1",1]])
  nb2_B_t[i] <- catch_e(coef(summary(z))[["treat_FCT_1",1]])
  nb2_p_i[i] <- catch_e(coef(summary(z))[["infect_FCT_1",4]])
  nb2_p_t[i] <- catch_e(coef(summary(z))[["treat_FCT_1",4]])
  
}
rm(z)

################################################################################

rm(
  test1, test2, test3, test4,
  test3_aic,
  test3_B_batch,
  test3_B_infect,
  test3_B_treat,
  test3_p_batch,
  test3_p_infect,
  test3_p_treat,
  i
)

#####
#The previous calculations can take several minutes
#The global environment variables are saved
#Refresh the save?

#Filter out errors regressions only; catch_e: FALSE option
if (0) {
  save.image(file = "./Data/glm_separate3a_eo.RData")
}

#Filter out both error and warning regressions; catch_e: TRUE option
if (0) {
  save.image(file = "./Data/glm_separate3a_ew.RData")
}

#####
#The graphing steps are in "./Scripts/glm_separate3b"
