#Run the setup script
source("Scripts/setup2.R")

################################################################################

library(ggrepel)
library(tidyverse)
library(edgeR)
library(EDASeq)

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
  function(mock0, mock1, cont0, cont1, treat0, treat1) {
    median(c(mock0, mock1, cont0, cont1, treat0, treat1), na.rm=TRUE)
  }
)

#CHECK: The first 2 medians in the previous vector should be
median(as.numeric(df_c_trans[1,]))
median(as.numeric(df_c_trans[2,]))

#...normalised
med_n <- pmap_dbl(df_n_trans, 
  function(mock0, mock1, cont0, cont1, treat0, treat1) {
    median(c(mock0, mock1, cont0, cont1, treat0, treat1), na.rm=TRUE)
  }
)

df_c_trans <- df_c_trans %>% add_column("med" = med_c)
df_n_trans <- df_n_trans %>% add_column("med" = med_n)

#Sample log(count+1) value
#minus median of 6 ... log(count+1) values

#New dataframe for the differences from the transcript's count median

df_c_diff <- df_c_trans %>% 
  transmute(
    mock0 = mock0 - med,
    mock1 = mock1 - med,
    cont0 = cont0 - med,
    cont1 = cont1 - med,
    treat0 = treat0 - med,
    treat1 = treat1 - med,
  )

df_n_diff <- df_n_trans %>% 
  transmute(
    mock0 = mock0 - med,
    mock1 = mock1 - med,
    cont0 = cont0 - med,
    cont1 = cont1 - med,
    treat0 = treat0 - med,
    treat1 = treat1 - med,
  )

# function to plot rle (without outliers)
# y = rle matrix
rle.plot <- function(y, ...) {
  boxplot(y, outline=FALSE, ...)
  abline(h=0)
}

#RLE plot using the previously defined function

#Shuffle columns to show samples from the same batch next to each other
df_c_diff <- df_c_diff[, c("cont0", "treat0", "mock0", 
  "cont1", "treat1", "mock1")]
df_n_diff <- df_n_diff[, c("cont0", "treat0", "mock0", 
  "cont1", "treat1", "mock1")]

rle.plot(df_c_diff, main = "RLE un-normalised", 
  col = c(rep("turquoise", times = 3), rep("red", times = 3)))
rle.plot(df_n_diff, main = "RLE normalised", 
  col = c(rep("turquoise", times = 3), rep("red", times = 3)))

################################################################################

#PCA plots using log(count+1) matrix (without subtraction of median)
pca_mtrx_c <- as.matrix(df_c_trans %>% dplyr::select(!med))
pca_mtrx_n <- as.matrix(df_n_trans %>% dplyr::select(!med))

plotPCA(pca_mtrx_c, main = "PCA, un-normalised, log(count+1)", isLog= TRUE)
plotPCA(pca_mtrx_n, main = "PCA, normalised, log(count+1)", isLog= TRUE)

#####

#Getting the calculation step from findMethods(plotPCA) : matrix
#Pass the coordinates to ggplot instead default plot for easier formatting
coord_PCA_fn <- function (object, k = 2, labels = TRUE, isLog = TRUE, ...) {
  if (!isLog) {
      Y <- apply(log(object + 1), 1, function(y) scale(y, 
          center = TRUE, scale = FALSE))
  }
  else {
      Y <- apply(object, 1, function(y) scale(y, center = TRUE, scale = FALSE))
  }
  s <- svd(Y)
  percent <- s$d^2/sum(s$d^2) * 100
  labs <- sapply(
    seq_along(percent), 
    function(i) {
      paste("PC ", i, " (", round(percent[i], 2), "%)", sep = "")
    }
  )
  return(c("s" = s, "labs" = labs, "x" = s$u[, 1], "y" = s$u[, 2]))
}

#####
#PCA using ggplot raw count
pca_c <- coord_PCA_fn(pca_mtrx_c)

coord_PCA_c_tb <- tibble(
  samplea = colnames(pca_mtrx_c),
  sampleb = c("infect&treat", "treat only", "none", "infect&treat", "treat only", "none"),
  x = c(pca_c$x1, pca_c$x2, pca_c$x3, pca_c$x4, pca_c$x5, pca_c$x6),
  y = c(pca_c$y1, pca_c$y2, pca_c$y3, pca_c$y4, pca_c$y5, pca_c$y6)
)

ggplot(data = coord_PCA_c_tb, aes(x=x, y=y)) +
  geom_point(
    fill = c("red", "red", "red", "blue", "blue", "blue"),
    shape = c(22,21,24,22,21,24),
    alpha = 0.5,
    size = 3
  ) +
  xlab(pca_c[["labs1"]]) + ylab(pca_c[["labs2"]]) +
  theme_bw() +
  geom_label_repel(aes(label = sampleb), color = "black", nudge_x = 0.1, nudge_y = 0.3) +
  ggtitle("PCA, un-normalised, log(count+1)") +
  coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))

#####
#PCA using ggplot normalised count