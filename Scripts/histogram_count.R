#Run the setup script
source("Scripts/setup1.R")

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
