################################################################################
#Script-file: power1.R
#Description: power calculation
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
library(RNASeqPower)

################################################################################
#Convenient summary for vectors
summary_v <- function(x) {
  mean <- mean(x, na.rm = TRUE)
  sd <- sd(x, na.rm = TRUE)
  var <- sd^2
  left <- quantile(
    x, probs = c(0, 0.001, 0.025, 0.25, 0.5, 0.75, 0.975, 0.999, 1), 
    na.rm = TRUE
  )
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

################################################################################
#Try a power calculation to check if the study is underpowered
#2 batches per exposure

depth_md <- median(y[["counts"]], na.rm = TRUE) #Check mean and median
depth_md 

depth_mn <- mean(y[["counts"]], na.rm = TRUE) #Check mean and median
depth_mn 

#Check: Sample size for effect 1 should be infinity
rnapower(depth=170, cv=0.4, effect=1.0, alpha=0.05, power=0.8)

#Check that the inverse works
rnapower(depth=170, cv=0.4, effect=2.0, alpha=0.05, power=0.8)
rnapower(depth=170, cv=0.4, n=5.419846, n2=5.419846, alpha=0.05, power=0.8)

#What sample size is required to detect effect (...) 
#... at alpha (...) and power (...)?

#...Median version
rnapower(
  depth = depth_md, 
  cv = 0.4, 
  effect = c(1.2, 2, 2.5, 3, 3.5, 4), 
  alpha = c(0.05, 0.01, 0.001), 
  power = c(0.8, 0.9)
)

#...Mean version
rnapower(
  depth = depth_mn, 
  cv = 0.4, 
  effect = c(1.2, 2, 2.5, 3, 3.5, 4), 
  alpha = c(0.05, 0.01, 0.001), 
  power = c(0.8, 0.9)
)

#Given sample size of 3 per group, what is the smallest effect size
#... that can be detected at alpha (...) and power (...)?
#...Median version
rnapower(
  depth = depth_md,
  cv = 0.4,
  n = 2,
  n2 = 2,
  alpha = c(0.05, 0.01, 0.001), 
  power = c(0.8, 0.9)
)

#...Mean version
rnapower(
  depth = depth_mn,
  cv = 0.4,
  n = 2,
  n2 = 2,
  alpha = c(0.05, 0.01, 0.001), 
  power = c(0.8, 0.9)
)

#depth is average count per gene (per sample?)
#effect is fold change 
#cv is usually estimated as 0.4 for humans
