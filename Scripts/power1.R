#Run the setup script
source("Scripts/setup1.R")

################################################################################

library(tidyverse)
library(RNASeqPower)

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

################################################################################
#Try a power calculation to check if the study is underpowered
#2 batches per exposure

depth <- median(y[["counts"]], na.rm = TRUE)
depth 

rnapower(depth=depth, cv=0.4, effect=1.2, alpha=0.05, power=0.8)
rnapower(depth=depth, cv=0.4, effect=2, alpha=0.05, power=0.8)
rnapower(depth=depth, cv=0.4, effect=2.5, alpha=0.05, power=0.8)
rnapower(depth=depth, cv=0.4, effect=3, alpha=0.05, power=0.8)


#depth is average count per gene (per sample?)
#effect is fold change 
#cv is usually estimated as 0.4 for humans
