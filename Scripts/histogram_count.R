################################################################################
#Script-file: histogram_count.R
#Description: histogram of counts (all / sample column specific)
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

library(tidyverse)

################################################################################
#Convenient summary for vectors
summary_v <- function(x) {
  mean <- mean(x, na.rm = TRUE)
  sd <- sd(x, na.rm = TRUE)
  var <- sd^2
  left <- quantile(x, probs = c(0, 0.001, 0.025, 0.25, 0.5, 0.75, 0.975, 0.999, 1), na.rm = TRUE)
  missing <- sum(is.na(x))
  present <- sum(!is.na(x))
  
  return(list("mean" = mean,
    "sd" = sd, "var"  = var, 
    "left" = left,
    "present" = present, "missing" = missing
    )
  )
}

################################################################################
#Summary of individual columns
summary_v(df_c)
for (i in 1:6) {
  z <- summary_v(df_c[,i])
  cat("column", i, "\n", sep = " ")
  print(z)
  rm(z)
}

#Summary of all columns combined
count_v <- as_tibble(c(df_c[,1], df_c[,2], df_c[,3], df_c[,4], df_c[,5], df_c[,6]))
colnames(count_v) <- "all_counts"
summary_v(as.matrix(count_v))

#Histogram of all counts
xran <- 2000
bwidth <- 25 
ggplot(data = count_v %>% filter(all_counts < xran)) +
  geom_histogram(aes(x = all_counts), binwidth = bwidth, boundary = 0)

#Histogram of column sepcific count
ggplot(data = as_tibble(df_c[,"mock0"]) %>% filter(value < xran)) +
  geom_histogram(aes(x = value), binwidth = bwidth, boundary = 0)
