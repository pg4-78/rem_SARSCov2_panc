################################################################################
#Script-file: diagnostic1.R
#Description: PCA, RLE diagnostic plots
################################################################################

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
  
  load(file = "./Data/setup2_co_2p01.RData")
  
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
  temp <- stack(as.data.frame(y))
  ggplot(temp) + 
  geom_boxplot(aes(x = ind, y = values), outlier.shape = NA, 
    fill = c(rep("#91bfdb", times = 3), rep("#fc8d59", times = 3))
  ) +
  theme_bw() +
  xlab("Sample") +
  ylab("Relative Log Expression") +
  geom_hline(yintercept=0, colour = "black", size = 0.5, linetype = "solid") +
  scale_y_continuous(breaks = seq(-2, 2, by = 0.5), minor_breaks = NULL) +
  scale_x_discrete(labels = c("infect&treat 1", "infect only 1", "none 1", "infect&treat 2", "infect only 2", "none 2")) +
  coord_cartesian(ylim = c(-2, 2))
}

#https://ggplot2.tidyverse.org/reference/geom_boxplot.html
#coef	Length of the whiskers as multiple of IQR. Defaults to 1.5.

#Shuffle columns to show samples from the same batch next to each other
df_c_diff <- df_c_diff[, c("treat0", "cont0", "mock0", 
  "treat1", "cont1", "mock1")]
df_n_diff <- df_n_diff[, c("treat0", "cont0", "mock0", 
  "treat1", "cont1", "mock1")]

#RLE plot using the previously defined function
rle.plot(df_c_diff) #un-normalised
rle.plot(df_n_diff) #normalised

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
  sampleb = c("infect&treat 1", "infect only 1", "none 1", "infect&treat 2", "infect only 2", "none 2"),
  x = c(pca_c$x1, pca_c$x2, pca_c$x3, pca_c$x4, pca_c$x5, pca_c$x6),
  y = c(pca_c$y1, pca_c$y2, pca_c$y3, pca_c$y4, pca_c$y5, pca_c$y6)
)

ggplot(data = coord_PCA_c_tb, aes(x=x, y=y)) +
  geom_point(
    fill = c("#91bfdb", "#91bfdb", "#91bfdb", "#fc8d59", "#fc8d59", "#fc8d59"),
    shape = c(22,21,24,22,21,24),
    alpha = 0.8,
    size = 3
  ) +
  xlab(pca_c[["labs1"]]) + ylab(pca_c[["labs2"]]) +
  theme_bw() +
  geom_label_repel(aes(label = sampleb), color = "black", nudge_x = 0.1, nudge_y = 0.3) +
  #ggtitle("PCA, un-normalised, log(count+1)") +
  coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))

#####
#PCA using ggplot normalised count

pca_n <- coord_PCA_fn(pca_mtrx_n)

coord_PCA_n_tb <- tibble(
  samplea = colnames(pca_mtrx_n),
  sampleb = c("infect&treat 1", "infect only 1", "none 1", "infect&treat 2", "infect only 2", "none 2"),
  x = c(pca_n$x1, pca_n$x2, pca_n$x3, pca_n$x4, pca_n$x5, pca_n$x6),
  y = c(pca_n$y1, pca_n$y2, pca_n$y3, pca_n$y4, pca_n$y5, pca_n$y6)
)

ggplot(data = coord_PCA_n_tb, aes(x=x, y=y)) +
  geom_point(
    fill = c("#91bfdb", "#91bfdb", "#91bfdb", "#fc8d59", "#fc8d59", "#fc8d59"),
    shape = c(22,21,24,22,21,24),
    alpha = 0.8,
    size = 3
  ) +
  xlab(pca_n[["labs1"]]) + ylab(pca_n[["labs2"]]) +
  theme_bw() +
  geom_label_repel(aes(label = sampleb), color = "black", nudge_x = 0.1, nudge_y = 0.3) +
  #ggtitle("PCA, normalised, log(count+1)") +
  coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))