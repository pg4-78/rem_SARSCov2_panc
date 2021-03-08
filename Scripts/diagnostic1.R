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
#Run the setup script
source("Scripts/setup1.R")

################################################################################
#MDS plots

#MDS plot without effective library size normalisation
plotMDS(y_alt, main="MDS without calcNormFactors")

#MDS plot with effective library size normalisation
plotMDS(y, main="MDS with calcNormFactors")

################################################################################
#log(count+1) transformation for RLE box plots 

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

#Shuffle columns to show samples from the same batch next to each other
df_c_diff <- df_c_diff[, c("cont1", "treat1", "mock1", 
  "cont2", "treat2", "mock2")]
df_n_diff <- df_n_diff[, c("cont1", "treat1", "mock1", 
  "cont2", "treat2", "mock2")]

rle.plot(df_c_diff, main = "RLE un-normalised", 
  col = c(rep("turquoise", times = 3), rep("red", times = 3)))
rle.plot(df_n_diff, main = "RLE normalised", 
  col = c(rep("turquoise", times = 3), rep("red", times = 3)))

################################################################################

#PCA plots using log(count+1) matrix (without subtraction of median)
pca_mtrx_c <- as.matrix(df_c_trans %>% dplyr::select(!med))
pca_mtrx_n <- as.matrix(df_n_trans %>% dplyr::select(!med))

plotPCA(pca_mtrx_c, main = "PCA, un-normalised, log(count+1)")
plotPCA(pca_mtrx_n, main = "PCA, normalised, log(count+1)")