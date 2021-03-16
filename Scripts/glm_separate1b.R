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
library(MASS)

################################################################################
#Load the results of "./Scripts/glm_separate1a.R"
load(file = "./Data/glm_separate1a.RData")

################################################################################
#Combine regression results into a table
#Counts; Names; AIC/ p-values

regr_tb <- as_tibble(df_c)
regr_tb <- regr_tb %>% 
  add_column(df_c_nm, 
    aic_po1_v, p_tr_po1_v, p_ba_po1_v,
    aic_po2_v, p_tr_po2_v,
    aic_nb1_v, p_tr_nb1_v, p_ba_nb1_v,
    aic_nb2_v, p_tr_nb2_v
  )

#Just keep rows which don't have NA in any of AIC/ p-values
regr_tb_cut <- regr_tb %>% drop_na()

#Check how many rows are remaining after the drop
#Original
dim(regr_tb)
#Stays
dim(regr_tb_cut)
#Leaves
dim(regr_tb)[1] - dim(regr_tb_cut)[1]

summary_v <- function(x, v = FALSE, ...) {
  mean <- mean(x, na.rm = TRUE)
  min <- min(x, na.rm = TRUE)
  max <- max(x, na.rm = TRUE)
  left <- quantile(x, probs = c(0,  0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE)
  missing <- sum(is.na(x))
  present <- sum(!is.na(x))
  
  return(list("mean" = mean,
    "min" = min, "max"  = max, 
    "left" = left,
    "present" = present, "missing" = missing
    )
  )
}

summary_v(aic_po1_v)
summary_v(aic_nb1_v)

####################
#______plot of the both sets of aic densities (pois vs neg-binom)
ggplot(data = regr_tb_cut) +
  geom_density(aes(x = aic_po1_v), color = "red") + 
  geom_density(aes(x = aic_nb1_v), color = "blue") +
  xlab("aic") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("AIC")

#Try again but with a narrower x-range
ggplot(data = regr_tb_cut) +
  geom_density(aes(x = aic_po1_v), color = "red") + 
  geom_density(aes(x = aic_nb1_v), color = "blue") +
  xlab("aic") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("AIC: zoom in") +
  coord_cartesian(xlim = c(0, 1000))

#
ggplot(data = regr_tb_cut) +
  geom_density(aes(x = aic_po1_v), color = "red") + 
  xlab("aic") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("AIC: poisson only") +
    coord_cartesian(xlim = c(0, 1000))

#
ggplot(data = regr_tb_cut) +
  geom_density(aes(x = aic_nb1_v), color = "blue") + 
  xlab("aic") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("AIC: neg binom only") +
    coord_cartesian(xlim = c(0, 170))

####################
#______plots of p values for treat coefficients 
#...(pois vs neg-binom) 
#...(batch coef vs none)

summary_v(p_tr_po1_v)
summary_v(p_tr_po2_v)
summary_v(p_tr_nb1_v)
summary_v(p_tr_nb2_v)

ggplot(data = regr_tb_cut) +
  geom_density(aes(x = p_tr_po1_v), color = "red") + 
  geom_density(aes(x = p_tr_nb1_v), color = "blue") +
  geom_density(aes(x = p_tr_po2_v), color = "red", linetype="dotted") + 
  geom_density(aes(x = p_tr_nb2_v), color = "blue", linetype="dotted") +
  xlab("p-value") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("Treat coefficient")

#Try again but with a narrower x-range
ggplot(data = regr_tb_cut) +
  geom_density(aes(x = p_tr_po1_v), color = "red") + 
  geom_density(aes(x = p_tr_nb1_v), color = "blue") +
  geom_density(aes(x = p_tr_po2_v), color = "red", linetype="dotted") + 
  geom_density(aes(x = p_tr_nb2_v), color = "blue", linetype="dotted") +
  xlab("p-value") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("Treat coefficient: zoom in") +
  coord_cartesian(xlim = c(0, 0.15))

#Neg binom only
ggplot(data = regr_tb_cut) +
  geom_density(aes(x = p_tr_nb1_v), color = "blue") +
  geom_density(aes(x = p_tr_nb2_v), color = "blue", linetype="dotted") +
  xlab("p-value") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("Treat coefficient: neg binom only") +
  coord_cartesian(xlim = c(0, 1))

####################
#______plots of p values for batch coefficients 
#...(pois vs neg-binom)

summary_v(p_ba_po1_v)
summary_v(p_ba_nb1_v)

ggplot(data = regr_tb_cut) +
  geom_density(aes(x = p_ba_po1_v), color = "red") + 
  geom_density(aes(x = p_ba_nb1_v), color = "blue") +
  xlab("p-value") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("Batch coefficient")

#Try again but with a narrower x-range
ggplot(data = regr_tb_cut) +
  geom_density(aes(x = p_ba_po1_v), color = "red") + 
  geom_density(aes(x = p_ba_nb1_v), color = "blue") +
  xlab("p-value") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("Batch coefficient: zoom in") +
  coord_cartesian(xlim = c(0, 0.15))

#Neg binom only
ggplot(data = regr_tb_cut) +
  geom_density(aes(x = p_ba_nb1_v), color = "blue") +
  xlab("p-value") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("Batch coefficient: neg binom only") +
  coord_cartesian(xlim = c(0, 1))
