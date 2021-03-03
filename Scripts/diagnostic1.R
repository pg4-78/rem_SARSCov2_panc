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
library(EDASeq)

################################################################################

#Data sub-folder inside the project
#CHANGE DEPENDING ON DATA LOCATION
raw <- readr::read_tsv("./Data/rnaseq_deseq_5dpi_counts_raw.tsv")

head(raw)

#Treat (SARS-Cov-2 + Remdesivir), Control (SARS-Cov-2), Mock (Uninfected)
raw <- dplyr::rename(raw, 
  "treat1" = 7, "cont1" = 8, "mock1" = 9, 
  "treat2" = 10, "cont2" = 11, "mock2" = 12, 
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

#MDS plot without effective library size normalisation
plotMDS(y, main="MDS without calcNormFactors")

#Normalisation factors
#Default
#y <- calcNormFactors(y)
y <- calcNormFactors(y, method = "none")

#MDS plot with effective library size normalisation
plotMDS(y, main="MDS with calcNormFactors")



################################################################################

#Matrix of raw counts (_c)
df_c <- y$counts

#Vector of normalisation factors (accounting for proportion imbalance)
s_a <- y$samples$norm.factors

#Vector of scaling factors (accounting for library size difference)
s_b <- colSums(df_c)
s_b <- s_b / mean(s_b) #not: s_b / mean(df_c)

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

#Perform the transformation log(count+1); raw & normalised versions
df_c_trans <- as_tibble(log(df_c+1))
df_n_trans <- as_tibble(log(df_n+1))

#Make vectors of the medians then append them to counts 
#...raw
med_c <- pmap_dbl(df_c_trans, 
  function(mock1, mock2, cont1, cont2, treat1, treat2) {
    median(c(mock1, mock2, cont1, cont2, treat1, treat2), na.rm=TRUE)
  }
)

#CHECK: The first 2 medians in the previous vector should be
median(as.numeric(df_c_trans[1,]))
median(as.numeric(df_c_trans[2,]))

#...normalised
med_n <- pmap_dbl(df_n_trans, 
  function(mock1, mock2, cont1, cont2, treat1, treat2) {
    median(c(mock1, mock2, cont1, cont2, treat1, treat2), na.rm=TRUE)
  }
)

df_c_trans <- df_c_trans %>% add_column("med" = med_c)
df_n_trans <- df_n_trans %>% add_column("med" = med_n)

#Sample log(count+1) value
#minus median of 6 ... log(count+1) values

#New dataframe for the differences from the transcript's count median

df_c_diff <- df_c_trans %>% 
  transmute(
    mock1 = mock1 - med,
    mock2 = mock2 - med,
    cont1 = cont1 - med,
    cont2 = cont2 - med,
    treat1 = treat1 - med,
    treat2 = treat2 - med,
  )

df_n_diff <- df_n_trans %>% 
  transmute(
    mock1 = mock1 - med,
    mock2 = mock2 - med,
    cont1 = cont1 - med,
    cont2 = cont2 - med,
    treat1 = treat1 - med,
    treat2 = treat2 - med,
  )

# function to plot rle (without outliers)
# y = rle matrix
rle.plot <- function(y, ...) {
  boxplot(y, outline=FALSE, ...)
  abline(h=0)
}

#RLE plot using the previously defined function
rle.plot(df_c_diff, main = "RLE un-normalised", 
  col = c("turquoise", "red", "turquoise", "red", "turquoise", "red"))
rle.plot(df_n_diff, main = "RLE normalised", 
  col = c("turquoise", "red", "turquoise", "red", "turquoise", "red"))

#PCA plots using log(count+1) matrix (without subtraction of median)
pca_mtrx_c <- as.matrix(df_c_trans %>% dplyr::select(!med))
pca_mtrx_n <- as.matrix(df_n_trans %>% dplyr::select(!med))

plotPCA(pca_mtrx_c, main = "PCA, un-normalised, log(count+1)")
plotPCA(pca_mtrx_n, main = "PCA, normalised, log(count+1)")