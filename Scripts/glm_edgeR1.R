#Run the setup script
source("Scripts/setup2.R")

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
#...With batch coefficient
design1 <- model.matrix(~batch_FCT_ + infect_FCT_ + treat_FCT_)
rownames(design1) <- colnames(as.matrix(y))
design1
#...Alternatively, without batch coefficient
design2 <- model.matrix(~infect_FCT_ + treat_FCT_)
rownames(design2) <- colnames(as.matrix(y))
design2

#Number of included variables
n_include <- dim(as.matrix(y))[1]

################################################################################
#______Options
#Use the batch coefficient? 
#...(TRUE: yes; FALSE: no)
opt_batch <- TRUE

#Use QL methods?: accounts for uncertainty in dispersion estimation
#...(TRUE: yes; FALSE: no)
opt_ql <- FALSE

################################################################################
#EdgeR models and gene/kegg/ontology tests

####################
#______Design matrix 
#...(number coefficients depending on batch coef absent/present)
if (opt_batch==TRUE) {
  design_m <- design1
  coef_b <- 2
  coef_i <- 3
  coef_t <- 4
} else {
  design_m <- design2
  coef_b <- NA
  coef_i <- 2
  coef_t <- 3
}

####################
#______Dispersion
#(edgeR GLM does not work until dispersions are estimated)
y <- estimateDisp(y, design = design_m)
plotBCV(y)

####################
#______Fit negative binomial glm
#Dispersion argument: If NULL will be extracted from y, with order of precedence: 
#...genewise dispersion, trended dispersions, common dispersion.
#QL or standard likelihood ratio option above

if (opt_ql == TRUE) {fit <- glmQLFit(y, design = design_m)} else
{fit <- glmFit(y, design = design_m)}

####################
#______Test INFECT coefficient
#...contrast takes priority over coef 
#...coef: to test if a coefficient is equal to zero (default last)
if (opt_ql == TRUE) {gene_test_i <- glmQLFTest(fit, coef = coef_i)} else
{gene_test_i <- glmLRT(fit, coef = coef_i)}

#Show top genes for INFECT
#...default is to display n=10
gene_top_i <- topTags(gene_test_i)

#over-represented gene GROUPS for INFECT
kegga_test_i <- kegga(gene_test_i, species="Hs")
kegga_top_i <- topKEGG(kegga_test_i, sort = "up")

goana_test_i <- goana(gene_test_i, species="Hs")
goana_top_i <- topGO(goana_test_i, sort = "up") 
#For reference:
#...Biological Process (BP)
#...Cellular Component (CC)
#...Molecular Function (MF)

cat("INFECT: Top genes\n")
gene_top_i
cat("INFECT: Top KEGG pathways\n")
kegga_top_i
cat("INFECT: Top GO\n")
goana_top_i

####################
#______Test TREAT coefficient
#...contrast takes priority over coef 
#...coef: to test if a coefficient is equal to zero (default last)
if (opt_ql == TRUE) {gene_test_t <- glmQLFTest(fit, coef = coef_t)} else
{gene_test_t <- glmLRT(fit, coef = coef_t)}

#Show top genes for TREAT
#...default is to display n=10
gene_top_t <- topTags(gene_test_t)

#over-represented gene GROUPS for TREAT
kegga_test_t <- kegga(gene_test_t, species="Hs")
kegga_top_t <- topKEGG(kegga_test_t, sort = "up")

goana_test_t <- goana(gene_test_t, species="Hs")
goana_top_t <- topGO(goana_test_t, sort = "up") 
#For reference:
#...Biological Process (BP)
#...Cellular Component (CC)
#...Molecular Function (MF)

cat("TREAT: Top genes\n")
gene_top_t
cat("TREAT: Top KEGG pathways\n")
kegga_top_t
cat("TREAT: Top GO\n")
goana_top_t

################################################################################
#Volcano plots

####################
#______INFECT:
#Extract the p-values
#They are stored in topTags: just increase the gene limit to all
temp_i <- topTags(gene_test_i, n = n_include)

#Extract the log fold changes
#Transform the log fold changes?
#Insert the p-values
#Transform the p-values
volc_i_tb <- as_tibble(fit[["coefficients"]]) %>% 
  dplyr::select("infect_FCT_1") %>% 
  add_column("infect_p" = temp_i[["table"]][["PValue"]]) %>% 
  mutate("infect_neg_log10_p" = -log10(infect_p))
colnames(volc_i_tb)[1] <- "infect_loge_FC"

#Volcano plot infect
ggplot(data = volc_i_tb) +
  geom_point(aes(x=infect_loge_FC, y=infect_neg_log10_p)) +
  geom_hline(yintercept = -log10(0.05), color = "red")
#Figure out how to colour and label

####################
#______TREAT:
#Extract the p-values
#They are stored in topTags: just increase the gene limit to all
temp_t <- topTags(gene_test_t, n = n_include)

#Extract the log fold changes
#Transform the log fold changes?
#Insert the p-values
#Transform the p-values
volc_t_tb <- as_tibble(fit[["coefficients"]]) %>% 
  dplyr::select("treat_FCT_1") %>% 
  add_column("treat_p" = temp_t[["table"]][["PValue"]]) %>% 
  mutate("treat_neg_log10_p" = -log10(treat_p))
colnames(volc_t_tb)[1] <- "treat_loge_FC"

#Volcano plot treat
ggplot(data = volc_t_tb) +
  geom_point(aes(x=treat_loge_FC, y=treat_neg_log10_p)) +
  geom_hline(yintercept = -log10(0.05), color = "red")
#Figure out how to colour and label

################################################################################
#P-value histograms

####################
#______INFECT:
ggplot(data = volc_i_tb, aes(x=infect_p)) +
  geom_histogram(binwidth=0.025, boundary=0) +
  ggtitle("p-value histogram infect")

####################
#______TREAT:
ggplot(data = volc_t_tb, aes(x=treat_p)) +
  geom_histogram(binwidth=0.025, boundary=0) +
  ggtitle("p-value histogram treat")