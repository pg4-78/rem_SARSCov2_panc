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
  
  load(file = "./Data/setup2_cut_2p01.RData")
  
  #No extra filter
  #file = "./Data/setup2_stock.RData"
  
  #Extra filter using opt_aveCPM_thresh
  #file = "./Data/setup2_cut_2pn3.RData"
  #file = "./Data/setup2_cut_2p01.RData"
  #file = "./Data/setup2_cut_2p03.RData"
  #file = "./Data/setup2_cut_2p05.RData"
  #average CPM cutoff at 2 to the power of (...)
  #-3 (below the lowest, no filter), 1, 3, 5
}

#Setup includes filtering and normalisation

#Just keep raw dataset, DGEList object y and aveLogCPM threshold
rm(list=setdiff(ls(), c("raw", "y", "opt_aveCPM_thresh")))

library(tidyverse)
library(ggrepel)
library(edgeR)
library(org.Hs.eg.db)

################################################################################
#Options

#Use the batch coefficient? 
#...(TRUE: yes; FALSE: no)
opt_batch <- TRUE

#Use QL methods?: accounts for uncertainty in dispersion estimation
#...(TRUE: yes; FALSE: no)
opt_ql <- FALSE

#Fold change threshold for being in top ranked genes
opt_fc_thresh <- 3

#Go to setup2 to change the average log CPM threshold
print(opt_aveCPM_thresh)

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
#EdgeR models and gene/kegg/ontology tests

####################
#______Design matrix 
#...(number coefficients depending on batch coef absent(FALSE)/ present(TRUE)
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
plotBCV(y, xlim = c(-2, 16))
axis(1, at = seq(-2, 16, by = 1), labels = NA)
axis(1, at = seq(0, 15, by = 5))

####################
#______Fit negative binomial glm
#Dispersion argument: If NULL will be extracted from y, with order of precedence: 
#...genewise dispersion, trended dispersions, common dispersion.
#QL or standard likelihood ratio option above

if (opt_ql == TRUE) {
  fit <- glmQLFit(y, design = design_m)
  plotQLDisp(fit) 
} else {
  fit <- glmFit(y, design = design_m)
}

####################
#______Test INFECT coefficient
#...contrast takes priority over coef 
#...coef: to test if a coefficient is equal to zero (default last)
if (opt_ql == TRUE) {gene_test_i <- glmQLFTest(fit, coef = coef_i)} else
{gene_test_i <- glmLRT(fit, coef = coef_i)}

#Show top genes for INFECT
#...default is to display n=10
gene_top_i <- topTags(gene_test_i, n = n_include)
dim(gene_top_i[["table"]])

temp_i <- gene_top_i
temp_i[["table"]] <- temp_i[["table"]][temp_i[["table"]][["logFC"]]>=log2(opt_fc_thresh),]
dim(temp_i[["table"]])

gene_top_i_tb <- as_tibble(gene_top_i[["table"]]) %>% 
  dplyr::filter(abs(logFC)>=log2(opt_fc_thresh)) %>% 
  dplyr::arrange(PValue)
head(gene_top_i_tb, n = 10)

#over-represented gene GROUPS for INFECT
kegga_test_i <- kegga(gene_test_i, species="Hs")
kegga_top_iu <- topKEGG(kegga_test_i, sort = "up")
kegga_top_id <- topKEGG(kegga_test_i, sort = "down")

goana_test_i <- goana(gene_test_i, species="Hs")
goana_top_iu <- topGO(goana_test_i, sort = "up") 
goana_top_id <- topGO(goana_test_i, sort = "down")
#For reference:
#...Biological Process (BP)
#...Cellular Component (CC)
#...Molecular Function (MF)

goana_top_iu <- goana_top_iu %>% mutate(P.Up = round(P.Up,3))##

cat("INFECT: Top genes\n")
gene_top_i
gene_top_i_tb
cat("INFECT: Top KEGG pathways\n")
kegga_top_iu
kegga_top_id
cat("INFECT: Top GO\n")
goana_top_iu
goana_top_id

####################
#______Test TREAT coefficient
#...contrast takes priority over coef 
#...coef: to test if a coefficient is equal to zero (default last)
if (opt_ql == TRUE) {gene_test_t <- glmQLFTest(fit, coef = coef_t)} else
{gene_test_t <- glmLRT(fit, coef = coef_t)}

#Show top genes for TREAT
#...default is to display n=10
gene_top_t <- topTags(gene_test_t, n = n_include)

temp_t <- gene_top_t
temp_t[["table"]] <- temp_t[["table"]][temp_t[["table"]][["logFC"]]>=log2(opt_fc_thresh),]
dim(temp_t[["table"]])

gene_top_t_tb <- as_tibble(gene_top_t[["table"]]) %>% 
  dplyr::filter(abs(logFC)>=log2(opt_fc_thresh)) %>% 
  dplyr::arrange(PValue) %>% 
head(gene_top_t_tb, n = 10)

#over-represented gene GROUPS for TREAT
kegga_test_t <- kegga(gene_test_t, species="Hs")
kegga_top_tu <- topKEGG(kegga_test_t, sort = "up")
kegga_top_td <- topKEGG(kegga_test_t, sort = "down")

goana_test_t <- goana(gene_test_t, species="Hs")
goana_top_tu <- topGO(goana_test_t, sort = "up")
goana_top_td <- topGO(goana_test_t, sort = "down")
#For reference:
#...Biological Process (BP)
#...Cellular Component (CC)
#...Molecular Function (MF)

cat("TREAT: Top genes\n")
gene_top_t
gene_top_t_tb
cat("TREAT: Top KEGG pathways\n")
kegga_top_tu
kegga_top_td
cat("TREAT: Top GO\n")
goana_top_tu
goana_top_td

################################################################################
#Volcano plots

#Extract the p-values
#Extract the log fold changes
#Transform the log fold changes?
#Insert the p-values
#Transform the p-values

####################
#______INFECT:

volc_i_tb <- as_tibble(gene_top_i[["table"]]) %>% 
  dplyr::mutate("neg_log10_p" = -log10(PValue))

ggplot(data = volc_i_tb, aes(x=logFC, y=neg_log10_p)) +
  geom_point(size = 1.5, alpha = 0.6, shape = 16) +
  geom_label_repel(data = volc_i_tb %>% filter(FDR<0.05 & (abs(logFC)>=log2(opt_fc_thresh))), #
    aes(label = gene_name), color = "red", alpha = 0.8, nudge_y = 1) +
  scale_x_continuous(breaks = seq(-10, 10, by = 5), minor_breaks = seq(-12, 12, by = 1)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), minor_breaks = NULL) +
  coord_cartesian(xlim = c(-12,12), ylim = c(0,10)) +
  theme_bw() +
  xlab("log2(fold change) infection") +
  ylab("-log10(p-value)") #+
  #ggtitle("volcano infect")

####################
#______TREAT:

volc_t_tb <- as_tibble(gene_top_t[["table"]]) %>% 
  dplyr::mutate("neg_log10_p" = -log10(PValue))

ggplot(data = volc_t_tb, aes(x=logFC, y=neg_log10_p)) +
  geom_point(size = 1.5, alpha = 0.6, shape = 16) +
  geom_label_repel(data = volc_t_tb %>% filter(FDR<0.06& (abs(logFC)>=log2(opt_fc_thresh))), # 
    aes(label = gene_name), color = "black", alpha = 0.8, nudge_y = 1) +
  scale_x_continuous(breaks = seq(-10, 10, by = 5), minor_breaks = seq(-12, 12, by = 1)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), minor_breaks = NULL) +
  coord_cartesian(xlim = c(-12,12), ylim = c(0,10)) +
  theme_bw() + 
  xlab("log2(fold change) treatment") +
  ylab("-log10(p-value)") #+
  #ggtitle("volcano treat") 

################################################################################
#P-value histograms

####################
#______INFECT:
ggplot(data = volc_i_tb, aes(x=PValue)) +
  geom_histogram(binwidth=0.025, boundary=0, fill="black", colour="white", size=0.2) +
  geom_density() +
  xlab("p-values of infect regression coefficients") +
  ylab("frequency") +
  coord_cartesian(x=c(0,1), y=c(0,500)) +
  scale_x_continuous(breaks = seq(0,1,0.2), minor_breaks = seq(0,1,0.025)) +
  theme_bw() #+ 
  #ggtitle("p-value histogram infect")

####################
#______TREAT: 
ggplot(data = volc_t_tb, aes(x=PValue)) +
  geom_histogram(binwidth=0.025, boundary=0, fill="black", colour="white", size=0.2) +
  xlab("p-values of treat regression coefficients") +
  ylab("frequency") +
  coord_cartesian(x=c(0,1), y=c(0,500)) +
  scale_x_continuous(breaks = seq(0,1,0.2), minor_breaks = seq(0,1,0.025)) +
  theme_bw() #+
  #ggtitle("p-value histogram treat")
