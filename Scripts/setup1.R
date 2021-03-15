library(tidyverse)
library(edgeR)

################################################################################

#Data sub-folder inside the project
raw <- readr::read_tsv("./Data/rnaseq_deseq_5dpi_counts_raw.tsv")

head(raw)

#Treat (SARS-Cov-2 + Remdesivir), Control (SARS-Cov-2), Mock (Uninfected)
raw <- dplyr::rename(raw, 
  "treat0" = 7, "cont0" = 8, "mock0" = 9, 
  "treat1" = 10, "cont1" = 11, "mock1" = 12, 
)

#2 Replicates for each of the 3 conditions
group <- c(1:3,1:3)
#edgeR Bioconductor compatible object
y <- DGEList(counts=raw[,7:12], group=group, genes=raw[,c(1,3)])

#Filter if expression is too low
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

#Filter out alternative splice forms (keep only predominant form)
#...Not yet

#Normalisation factors

y <- calcNormFactors(y)
y_alt <- calcNormFactors(y, method = "none")

################################################################################

#Matrix of raw counts (_c)
df_c <- y$counts
row.names(df_c) <- y[["genes"]][["gene_name"]]

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
for (i in 1:dim(df_n)[2]){
  df_n[,i] <- df_n[,i]/s[i]
}