library(data.table)
library(tidyverse)
library(edgeR)
library(MASS)

################################################################################
#Run the setup script
source("Scripts/setup1.R")

################################################################################

#Make factor variables
batch_FCT_ <- factor(c(0,0,0,1,1,1))

#...levels: baseline should go first
expos_FCT_ <- factor(
  c("treat","cont","mock","treat","cont","mock"),
  levels = c("cont","treat","mock")
)

tibble(sample = colnames(df_c), batch_FCT_, expos_FCT_)

#Check design matrices
design1 <- model.matrix(~batch_FCT_ + expos_FCT_)
rownames(design1) <- colnames(df_c)
design1

design2 <- model.matrix(~expos_FCT_)
rownames(design2) <- colnames(df_c)
design2

################################################################################
#For ease of troubleshooting: Before attempting to loop through every gene, 
#...try one gene

#Choose a gene between 1 and (...18302) (inclusive)
dim(df_c)[1]

df_c_tb <- as_tibble(df_c) %>% 
  add_column("gene_name" = df_c_nm) %>% 
  add_column("row_num" = as.numeric(1:dim(df_c)[1]))

#Error genes: 22: "KLHL17"
gene_num <- 22

#check offset

#glm:poisson and glm.nb require integer outcome values

####################
#______Poisson with batch as a covariate - check gene_num 
#test1b: also check if offset option can be written 
#...as part of the formula instead

test1 <- glm(
  formula = df_c[gene_num,] ~ batch_FCT_ + expos_FCT_,
  family = "poisson",
  offset = log(s_a*og_lib_sz_v)
)
test1
coef(summary(test1))

test1b <- glm(
  formula = df_c[gene_num,] ~ batch_FCT_ + expos_FCT_ + offset(log(s_a*og_lib_sz_v)),
  family = "poisson"
)
test1b
coef(summary(test1b))

####################
#______Poisson without batch as a covariate - check gene_num 
#test2b: also check if offset option can be written 
#...as part of the formula instead

test2 <- glm(
  formula = df_c[gene_num,] ~ expos_FCT_,
  family = "poisson",
  offset = log(s_a*og_lib_sz_v)
)
test2
coef(summary(test2))

test2b <- glm(
  formula = df_c[gene_num,] ~ expos_FCT_ + offset(log(s_a*og_lib_sz_v)),
  family = "poisson"
)
test2b
coef(summary(test2b))

####################
#______Neg binom with batch as a covariate - check gene_num 

test3 <- tryCatch(
  glm.nb(formula = df_c[gene_num,] ~ batch_FCT_ + expos_FCT_ + offset(log(s_a*og_lib_sz_v))), 
  error = function(e) NULL
)
test3
tryCatch(coef(summary(test3)), error = function(e) NA)
test3_aic <- tryCatch(test3[["aic"]], error = function(e) NA)
test3_treat <- tryCatch(coef(summary(test3))[["expos_FCT_treat",4]], error = function(e) NA)
test3_batch <- tryCatch(coef(summary(test3))[["batch_FCT_1",4]], error = function(e) NA)
test3_aic
test3_treat
test3_batch

####################
#______Neg binom without batch as a covariate - check gene_num 

test4 <- tryCatch(
  glm.nb(formula = df_c[gene_num,] ~ expos_FCT_ + offset(log(s_a*og_lib_sz_v))), 
  error = function(e) NULL
)
test4
tryCatch(coef(summary(test4)), error = function(e) NA)
test4_aic <- tryCatch(test4[["aic"]], error = function(e) NA)
test4_treat <- tryCatch(coef(summary(test4))[["expos_FCT_treat",4]], error = function(e) NA)
test4_aic
test4_treat

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
aic_po1_v <- as.numeric(rep(NA, dim(df_c)[1]))
p_tr_po1_v <- as.numeric(rep(NA, dim(df_c)[1]))
p_ba_po1_v <- as.numeric(rep(NA, dim(df_c)[1]))

for (i in 1:dim(df_c)[1]) {
  z <- glm(
    formula = df_c[i,] ~ batch_FCT_ + expos_FCT_,
    family = "poisson",
    offset = log(s_a*og_lib_sz_v)
  )
  
  aic_po1_v[i] <- z[["aic"]]
  p_tr_po1_v[i] <- coef(summary(z))[["expos_FCT_treat",4]]
  p_ba_po1_v[i] <- coef(summary(z))[["batch_FCT_1",4]]
}
rm(z)


####################
#______Poisson without batch coefficient
aic_po2_v <- as.numeric(rep(NA, dim(df_c)[1]))
p_tr_po2_v <- as.numeric(rep(NA, dim(df_c)[1]))

for (i in 1:dim(df_c)[1]) {
  z <- glm(
    formula = df_c[i,] ~ expos_FCT_,
    family = "poisson",
    offset = log(s_a*og_lib_sz_v)
  )
  
  aic_po2_v[i] <- z[["aic"]]
  p_tr_po2_v[i] <- coef(summary(z))[["expos_FCT_treat",4]]
}
rm(z)


####################
#______Negative binomial with batch coefficient
aic_nb1_v <- as.numeric(rep(NA, dim(df_c)[1]))
p_tr_nb1_v <- as.numeric(rep(NA, dim(df_c)[1]))
p_ba_nb1_v <- as.numeric(rep(NA, dim(df_c)[1]))

for (i in 1:dim(df_c)[1]) {
  z <- tryCatch(
    glm.nb(formula = df_c[i,] ~ batch_FCT_ + expos_FCT_ + offset(log(s_a*og_lib_sz_v))), 
    error = function(e) {NA}
)
  
  aic_nb1_v[i] <- tryCatch(z[["aic"]], error = function(e){NA})
  p_tr_nb1_v[i] <- tryCatch(coef(summary(z))[["expos_FCT_treat",4]], error = function(e){NA})
  p_ba_nb1_v[i] <- tryCatch(coef(summary(z))[["batch_FCT_1",4]], error = function(e){NA})
}
rm(z)
  

####################
#______Negative binomial without batch coefficient
aic_nb2_v <- as.numeric(rep(NA, dim(df_c)[1]))
p_tr_nb2_v <- as.numeric(rep(NA, dim(df_c)[1]))

for (i in 1:dim(df_c)[1]) {
  z <- tryCatch(
    glm.nb(formula = df_c[i,] ~ expos_FCT_ + offset(log(s_a*og_lib_sz_v))), 
    error = function(e){NA}
)
  
  aic_nb2_v[i] <- tryCatch(z[["aic"]], error = function(e){NA})
  p_tr_nb2_v[i] <- tryCatch(coef(summary(z))[["expos_FCT_treat",4]], error = function(e){NA})
}
rm(z)

################################################################################
#The previous calculations can take several minutes
#The global environment variables are saved
#Refresh the save?
if (0) {
  save.image(file = "./Data/glm_separate1a.RData")
}

#The graphing steps are in "./Scripts/glm_separate1b"
