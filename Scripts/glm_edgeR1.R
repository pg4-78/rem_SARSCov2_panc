#Run the setup script
source("Scripts/setup1.R")

#Setup includes filtering and normalisation

#Just keep raw dataset and DGEList object y
rm(list=setdiff(ls(), c("raw", "y")))

library(tidyverse)
library(edgeR)

################################################################################
#Make factor variables

batch_FCT_ <- factor(c(0,0,0,1,1,1))
infect_FCT_ <- factor(c(1,1,0,1,1,0))
treat_FCT_ <- factor(c(1,0,0,1,0,0))

tibble(sample = colnames(y[["counts"]]), batch_FCT_, infect_FCT_, treat_FCT_)

#Check design matrices
design1 <- model.matrix(~batch_FCT_ + infect_FCT_ + treat_FCT_)
rownames(design1) <- colnames(df_c)
design1

design2 <- model.matrix(~infect_FCT_ + treat_FCT_)
rownames(design2) <- colnames(df_c)
design2

################################################################################
#Estimate the dispersion, other options are left default for now
y_spare <- y
########################################
#___with batch as a covariate

####################
#______Standard glmFit

#The glmFit does not work until dispersions are estimated

#Error in glmFit.DGEList(y, design = design1) : 
#  No dispersion values found in DGEList object.

#Dispersion
y <- estimateDisp(y, design = design1)

#The 4 coefficients are intercept, (batch, infect, treat)_FCT_1
design1

#Standard negative binomial glm
#Dispersion argument: If NULL will be extracted from y, with order of precedence: 
#...genewise dispersion, trended dispersions, common dispersion.
fit_st <- glmFit(y, design = design1)

#contrast takes priority over coef
#coef: to test if a coefficient is equal to zero (default last)
test_lr_i <- glmLRT(fit_st, coef = 3)
top_lr_i <- topTags(test_lr_i)
top_lr_i

#Made sure the alternate way of specifying is equivalent: OK
if (0) {
  test_lr_i_alt <- glmLRT(fit_st, contrast = c(0,0,1,0))
  top_lr_i_alt <- topTags(test_lr_i_alt)
  top_lr_i_alt
}

####################
#______Quasi-Likelihood dispersion: 

fit_ql <- glmQLFit(y, design = design1)
plotQLDisp(fit_ql)


test_ql_i <- glmLRT(fit_ql, coef = 3)

####################
#______Gene Ontology
ont_kegga_lr_i <- kegga(test_lr_i)
ont_goana_lr_i <- goana(test_lr_i)

ont_kegga_lr_i
ont_goana_lr_i
########################################
#___without batch as a covariate y_xb (just in case)]
if (0) {
  y_xb <- estimateDisp(y, design = design2)
  plotBCV(y_xb)
  
  #...
}
