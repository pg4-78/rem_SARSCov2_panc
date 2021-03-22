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

####################
#______Standard glmFit

#(The glmFit does not work until dispersions are estimated)
#Dispersion
y <- estimateDisp(y, design = design1)

#The 4 coefficients are intercept, (batch, infect, treat)_FCT_1
design1

#Standard negative binomial glm
#Dispersion argument: If NULL will be extracted from y, with order of precedence: 
#...genewise dispersion, trended dispersions, common dispersion.
fit_st <- glmFit(y, design = design1)

#Test INFECT coefficient
#...contrast takes priority over coef 
#...contrast = c(0,0,1,0) gives equivalent results
#...coef: to test if a coefficient is equal to zero (default last)
gene_test_lr_i <- glmLRT(fit_st, coef = 3)
#Show top genes for INFECT
#...default is to display n=10
gene_top_lr_i <- topTags(gene_test_lr_i) 

#over-represented gene GROUPS for INFECT
kegga_test_lr_i <- kegga(gene_test_lr_i, species="Hs")
kegga_top_lr_i <- topKEGG(kegga_test_lr_i, sort = "up")

goana_test_lr_i <- goana(gene_test_lr_i, species="Hs")
goana_top_lr_i <- topGO(goana_test_lr_i, sort = "up") 
#For reference:
#...Biological Process (BP)
#...Cellular Component (CC)
#...Molecular Function (MF)

cat("INFECT: Top genes\n")
gene_top_lr_i
cat("INFECT: Top KEGG pathways\n")
kegga_top_lr_i
cat("INFECT: Top GO\n")
goana_top_lr_i

#####

#Test TREAT coefficient
#...contrast takes priority over coef 
#...contrast = c(0,0,1,0) gives equivalent results
#...coef: to test if a coefficient is equal to zero (default last)
gene_test_lr_t <- glmLRT(fit_st, coef = 4)
#Show top genes for TREAT
#...default is to display n=10
gene_top_lr_t <- topTags(gene_test_lr_t) 

#over-represented gene GROUPS for TREAT
kegga_test_lr_t <- kegga(gene_test_lr_t, species="Hs")
kegga_top_lr_t <- topKEGG(kegga_test_lr_t, sort = "up")

goana_test_lr_t <- goana(gene_test_lr_t, species="Hs")
goana_top_lr_t <- topGO(goana_test_lr_t, sort = "up") 
#For reference:
#...Biological Process (BP)
#...Cellular Component (CC)
#...Molecular Function (MF)

cat("TREAT: Top genes\n")
gene_top_lr_t
cat("TREAT: Top KEGG pathways\n")
kegga_top_lr_t
cat("TREAT: Top GO\n")
goana_top_lr_t

####################
#______Quasi-Likelihood dispersion: 
# ... ???

fit_ql <- glmQLFit(y, design = design1)
plotQLDisp(fit_ql)


test_ql_i <- glmLRT(fit_ql, coef = 3)

########################################
#___without batch as a covariate y_xb (just in case)]
# ... ???

################################################################################

#INFECT: Volcano plots
#Extract the p-values
#Transform the p-values
#Extract the log? fold changes
#Transform the log? fold changes
#Figure out how to colour
as_tibble()

#INFECT: p-value Histograms/ density plots
#extract p-values
