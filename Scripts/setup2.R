################################################################################
#Script-file: setup2.R
#Description: import data, filter, normalise, 
#... count matrices for PCA, RLE diagnostics
################################################################################

#Clear
#...variables
rm(list=ls())
#...console
cat("\014\n")
#...graphs
tryCatch(dev.off(), error = function(e) {NULL})
dev.new()

library(tidyverse)
library(edgeR)
library(org.Hs.eg.db)

################################################################################

#Data sub-folder inside the project
raw <- readr::read_tsv("./Data/rnaseq_deseq_5dpi_counts_raw.tsv")

head(raw)

#Treat (SARS-Cov-2 + Remdesivir), Control (SARS-Cov-2), Mock (Uninfected)
raw <- dplyr::rename(raw, "ensembl_gene_id" = gene_id,
  "treat0" = 7, "cont0" = 8, "mock0" = 9, 
  "treat1" = 10, "cont1" = 11, "mock1" = 12, 
)

#2 Replicates for each of the 3 conditions
group <- c(1:3,1:3)
#edgeR Bioconductor compatible object
y <- DGEList(counts=raw[,7:12], group=group, genes=raw[,c(1,3)])
print(dim(y))
#CPM raw: without any any filters
log2_cpm_raw_v <- aveLogCPM(y)

#Filter if expression is too low
if (FALSE) {
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  print(dim(y))
}

#Extra filter on top of default (?)
#Average Log(base2) CPM 
#Subset: require cpm >= opt_aveCPM_thresh
if (TRUE) {
  log2_cpm_pre_v <- aveLogCPM(y)
  opt_aveCPM_thresh <- 2^1

  keep_cpm <- log2_cpm_pre_v >= log2(opt_aveCPM_thresh)
  y <- y[keep_cpm, , keep.lib.sizes=FALSE]
  print(dim(y))
  log2_cpm_post_v <- aveLogCPM(y)
  
  #...check result
  cat(c("average log2 CPM before filter:", round(mean(log2_cpm_pre_v), 3), "\n"))
  cat(c("average log2 CPM after filter:", round(mean(log2_cpm_post_v), 3), "\n"))
  
  cat(c("min log2 CPM before filter:", round(min(log2_cpm_pre_v), 3), "\n"))
  cat(c("min log2 CPM after filter:", round(min(log2_cpm_post_v), 3), "\n"))
}

################################################################################
#Distribution of CPM before any genes were filtered
log2_cpm_raw_tb <- tibble("x" = log2_cpm_raw_v, "gene_biotype" = raw$gene_biotype)

#all genes
ggplot(data = log2_cpm_raw_tb, aes(x=x)) +
  geom_histogram(boundary = 0, binwidth = 0.5, fill="black", colour="white", size=0.2) +
  scale_x_continuous(breaks = seq(-6,16,2), minor_breaks = seq(-6,16,0.5)) +
  xlab("log2(Average CPM)") +
  ylab("frequency") +
  theme_bw() +
  geom_hline(yintercept=0, colour="black")

#only protein coding genes
ggplot(data = log2_cpm_raw_tb %>% filter(gene_biotype=="protein_coding"), aes(x=x)) +
  geom_histogram(boundary = 0, binwidth = 0.5, fill="black", colour="white", size=0.2) +
  scale_x_continuous(breaks = seq(-6,16,2), minor_breaks = seq(-6,16,0.5)) +
  xlab("log2(Average CPM)") +
  ylab("frequency") +
  theme_bw() +
  geom_hline(yintercept=0, colour="black")

################################################################################
#Filter out alternative splice forms (keep only predominant form)
#...No

#Normalisation factors
y <- calcNormFactors(y)
y$samples
y_alt <- calcNormFactors(y, method = "none")
y_alt$samples

egENSEMBL_tb <- toTable(org.Hs.egENSEMBL)
#Gene id numbers added to the edgeR DGEList object
m <- match(y[["genes"]][["ensembl_gene_id"]], egENSEMBL_tb$ensembl_id)
y$genes$EntrezGene <- egENSEMBL_tb$gene_id[m]

################################################################################

#Matrix of raw counts
df_c <- as.matrix(y$counts)
df_c_nm <- as_tibble(y[["genes"]][["gene_name"]])
colnames(df_c_nm) <- "gene_name"

#Vector of normalisation factors (accounting for proportion imbalance)
s_a <- y$samples$norm.factors

#Vector of scaling factors (accounting for library size difference)
og_lib_sz_v <- colSums(df_c)
s_b <- og_lib_sz_v / mean(og_lib_sz_v) #not: s_b / mean(df_c)

#Product of previous scaling factors for the final scaling factor
s <- s_a*s_b

#Matrix of counts normalised by column (_n)
#edgeR user guide section 2.8.3
#Norm factor below 1 
#...library scaled down
#...effectively up-scaling counts

#Get count matrix
df_n <- df_c

#Divide each column by the respective scaling factor
for (i in 1:dim(df_n)[2]) {
  df_n[,i] <- df_n[,i]/s[i]
}

################################################################################
#Save?
if (FALSE) {
  save.image(file = "./Data/setup2_co_2p01.RData")
  
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
