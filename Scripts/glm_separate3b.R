################################################################################
#Script-file: glm_separate3b.R
#Description: PART B of poisson vs negative binomial
################################################################################

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

#Two new dataframes/ tibbles extracting regression 
#...(poisson/neg binom) (batch coeff present/absent)
#...1st dataframe excludes genes with at least 1 error/warning regression
#...2nd dataframe excludes genes with at least 1 error regression (warnings allowed)

################################################################################
#Load the results of "./Scripts/glm_separate3a.R"
#Filter out error and warning regressions
load(file = "./Data/glm_separate3a_ew.RData")

#Combine regression results into a table
#Counts; Names; AIC; batch/infect/treat: p-values/coefficients

tb_ew <- as_tibble(df_c)
tb_ew <- tb_ew %>% 
  add_column(df_c_nm, log2_cpm_post_v,
    po1_aic, po1_B_b, po1_B_i, po1_B_t, po1_p_b, po1_p_i, po1_p_t,
    po2_aic, po2_B_i, po2_B_t, po2_p_i, po2_p_t,
    nb1_aic, nb1_B_b, nb1_B_i, nb1_B_t, nb1_p_b, nb1_p_i, nb1_p_t,
    nb2_aic, nb2_B_i, nb2_B_t, nb2_p_i, nb2_p_t
  )

#Just keep rows which don't have NA in any of AIC/ p-values/ coefficients
tb_ew_keep <- tb_ew %>% drop_na()
tb_ew_drop <- anti_join(tb_ew, tb_ew_keep)

################################################################################
rm(list=setdiff(ls(), c("tb_ew", "tb_ew_keep", "tb_ew_drop")))
#Load the results of "./Scripts/glm_separate3a.R"
#Filter out error regressions
load(file = "./Data/glm_separate3a_eo.RData")

#Combine regression results into a table
#Counts; Names; AIC; batch/infect/treat: p-values/coefficients

tb_eo <- as_tibble(df_c)
tb_eo <- tb_eo %>% 
  add_column(df_c_nm, log2_cpm_post_v,
    po1_aic, po1_B_b, po1_B_i, po1_B_t, po1_p_b, po1_p_i, po1_p_t,
    po2_aic, po2_B_i, po2_B_t, po2_p_i, po2_p_t,
    nb1_aic, nb1_B_b, nb1_B_i, nb1_B_t, nb1_p_b, nb1_p_i, nb1_p_t,
    nb2_aic, nb2_B_i, nb2_B_t, nb2_p_i, nb2_p_t
  )

#Just keep rows which don't have NA in any of AIC/ p-values/ coefficients
tb_eo_keep <- tb_eo %>% drop_na()
tb_eo_drop <- anti_join(tb_eo, tb_eo_keep)

################################################################################
#Check how many rows are remaining after the drop
#18302, 11364, 6938: When filtering errors and warnings
#Original
dim(tb_ew)[1]
#Stays
dim(tb_ew_keep)[1]
#Leaves
dim(tb_ew)[1] - dim(tb_ew_keep)[1]

#18302, 18275, 27: When filtering errors only
#Original
dim(tb_eo)[1]
#Stays
dim(tb_eo_keep)[1]
#Leaves
dim(tb_eo)[1] - dim(tb_eo_keep)[1]

#Clean up
rm(list=setdiff(ls(), c("tb_ew", "tb_ew_keep", "tb_ew_drop", 
  "tb_eo", "tb_eo_keep", "tb_eo_drop")))

#Save these regression result tables
#Replace current save?
if (0) {
  save.image(file = "./Data/glm_separate3b_tb.RData")
}

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
#Summarise average log2 CPM by keep/drop status
#for each gene check if any of its regressions has errors/warnings
print("excluded due to error:")
summary_v(tb_eo_drop$log2_cpm_post_v)
print("----------")

print("excluded due to error/warning:")
summary_v(tb_ew_drop$log2_cpm_post_v)
print("----------")

print("genes without error regressions:")
summary_v(tb_eo_keep$log2_cpm_post_v)
print("----------")

print("genes without errors/warning regressions:")
summary_v(tb_ew_keep$log2_cpm_post_v)
print("----------")

ggplot() +
  geom_histogram(data = tb_ew_keep, aes(x = log2_cpm_post_v),
    boundary = 0, binwidth = 0.5) +
  scale_x_continuous(minor_breaks = seq(-2, 15, 1), breaks = seq(0, 15, 5)) +
  coord_cartesian(xlim = c(-2,15)) + 
  ggtitle("average log2 CPM: without error/warning regressions") +
  theme_bw()

ggplot() +
  geom_histogram(data = tb_ew_drop, aes(x = log2_cpm_post_v),
    boundary = 0, binwidth = 0.5) +
  scale_x_continuous(minor_breaks = seq(-2, 15, 1), breaks = seq(0, 15, 5)) +
  coord_cartesian(xlim = c(-2,15)) + 
  ggtitle("average log2 CPM: with error/warning regressions") +
  theme_bw()

ggplot() +
  geom_histogram(data = tb_eo_keep, aes(x = log2_cpm_post_v),
    boundary = 0, binwidth = 0.5) +
  scale_x_continuous(minor_breaks = seq(-2, 15, 1), breaks = seq(0, 15, 5)) +
  coord_cartesian(xlim = c(-2,15)) + 
  ggtitle("average log2 CPM: without error regressions") +
  theme_bw()

ggplot() +
  geom_histogram(data = tb_eo_drop, aes(x = log2_cpm_post_v),
    boundary = 0, binwidth = 0.5) +
  scale_x_continuous(minor_breaks = seq(-2, 15, 1), breaks = seq(0, 15, 5)) +
  coord_cartesian(xlim = c(-2,15)) + 
  ggtitle("average log2 CPM: with error regressions") +
  theme_bw()

################################################################################
#______Kernel density plots of AIC

#Exclude errors only: poisson/ neg binom: with/out batch coefficient
#Dashed excludes batch coefficient
map(list(tb_eo_keep$po1_aic, tb_eo_keep$po2_aic, tb_eo_keep$nb1_aic, tb_eo_keep$nb2_aic), summary_v)
xran <- 1000
ggplot() +
  geom_density(data = tb_eo_keep %>% filter(po1_aic<xran), aes(x = po1_aic), color = "orange", size = 0.8) +
  geom_density(data = tb_eo_keep %>% filter(po2_aic<xran), aes(x = po2_aic), color = "blue", linetype = "dotted", size = 0.8) +
  geom_density(data = tb_eo_keep %>% filter(nb1_aic<xran), aes(x = nb1_aic), color = "black", size = 0.8) +
  geom_density(data = tb_eo_keep %>% filter(nb2_aic<xran), aes(x = nb2_aic), color = "red", linetype = "dotted", size = 0.8) +
  xlab("AIC") +
  ylab("kernel probability density est.") +
  theme_bw() +
  #ggtitle("AIC density: exclude error regr. only") +
  coord_cartesian(xlim = c(0, xran))

#####

#Exclude errors and warnings: poisson/ neg binom: with/out batch coefficient
map(list(tb_ew_keep$po1_aic, tb_ew_keep$po2_aic, tb_ew_keep$nb1_aic, tb_ew_keep$nb2_aic), summary_v)
xran <- 1000
ggplot() +
  geom_density(data = tb_ew_keep %>% filter(po1_aic<xran), aes(x = po1_aic), color = "orange", size = 0.8) +
  geom_density(data = tb_ew_keep %>% filter(po2_aic<xran), aes(x = po2_aic), color = "blue", linetype = "dashed", size = 0.8) +
  geom_density(data = tb_ew_keep %>% filter(nb1_aic<xran), aes(x = nb1_aic), color = "black", size = 0.8) +
  geom_density(data = tb_ew_keep %>% filter(nb2_aic<xran), aes(x = nb2_aic), color = "red", linetype = "dashed", size = 0.8) +
  xlab("AIC") +
  ylab("Kernel Probability Density Est.") +
  theme_bw() +
  #ggtitle("AIC density: exclude error & warning regr.") +
  coord_cartesian(xlim = c(0,xran))


################################################################################
#______Two way scatterplot of AIC

#negative binomial vs poisson
#Errors and warnings excluded
#Batch covariate used
ggplot(data = tb_ew_keep %>% 
    dplyr::mutate(ln_po1_aic = log10(po1_aic)) %>% 
    dplyr::mutate(ln_nb1_aic = log10(nb1_aic)) 
  ) +
  geom_point(aes(x=ln_po1_aic, y=ln_nb1_aic), size=0.5, alpha=0.4, pch = 16) +
  geom_abline(slope = 1, color = "red") +
  xlab("poisson log10(AIC)") +
  ylab("negative binomial log10(AIC)") +
  coord_cartesian(xlim = c(1,5), ylim = c(1,5)) +
  theme_bw()

################################################################################
#______plots of p values for infect coefficients 
#...(pois vs neg-binom) 
#...(batch coef vs none)
map(list(tb_eo_keep$po1_p_i, tb_eo_keep$po2_p_i, tb_eo_keep$nb1_p_i, tb_eo_keep$nb2_p_i), summary_v)

ggplot(data = tb_eo_keep) +
  geom_density(aes(x = po1_p_i), color = "red") +
  geom_density(aes(x = po2_p_i), color = "red", linetype = "dotted") +
  xlab("p-value") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("p-value infect: exclude error regr.") +
  coord_cartesian(ylim = NULL)

ggplot(data = tb_eo_keep) +
  geom_density(aes(x = nb1_p_i), color = "blue") +
  geom_density(aes(x = nb2_p_i), color = "blue", linetype = "dotted") +
  xlab("p-value") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("p-value infect: exclude error regr.") +
  coord_cartesian(ylim = NULL)

xran <- 0.25
ggplot(data = tb_eo_keep %>% filter(nb1_p_i < xran)) +
  geom_histogram(aes(x = nb1_p_i), binwidth = 0.01, boundary=0) +
  xlab("p-value") +
  ylab("count") +
  theme_bw() +
  ggtitle("p-value infect: exclude error regr.") +
  coord_cartesian(ylim = NULL)

# !!
ggplot(data = tb_ew_keep) +
  geom_histogram(aes(x = po1_p_i), binwidth = 0.025, boundary=0) +
  xlab("p-value infect (poisson)") +
  ylab("count") +
  theme_bw() +
  #ggtitle("p-value treat: exclude error regr.") +
  coord_cartesian(ylim = NULL)


####################
#______plots of p values for treat coefficients 
#...(pois vs neg-binom) 
#...(batch coef vs none)
map(list(tb_eo_keep$po1_p_t, tb_eo_keep$po2_p_t, tb_eo_keep$nb1_p_t, tb_eo_keep$nb2_p_t), summary_v)

ggplot(data = tb_eo_keep) +
  geom_density(aes(x = po1_p_t), color = "red") +
  geom_density(aes(x = po2_p_t), color = "red", linetype = "dotted") +
  xlab("p-value") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("p-value treat: exclude error regr.") +
  coord_cartesian(ylim = NULL)

ggplot(data = tb_eo_keep) +
  geom_density(aes(x = nb1_p_t), color = "blue") +
  geom_density(aes(x = nb2_p_t), color = "blue", linetype = "dotted") +
  xlab("p-value") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("p-value treat: exclude error regr.") +
  coord_cartesian(ylim = NULL)

xran <- 0.25
ggplot(data = tb_eo_keep %>% filter(nb1_p_t < xran)) +
  geom_histogram(aes(x = nb1_p_t), binwidth = 0.01, boundary=0) +
  xlab("p-value") +
  ylab("count") +
  theme_bw() +
  ggtitle("p-value treat: exclude error regr.") +
  coord_cartesian(ylim = NULL)

# !!
ggplot(data = tb_ew_keep) +
  geom_histogram(aes(x = po1_p_t), binwidth = 0.025, boundary=0) +
  xlab("p-value treat (poisson)") +
  ylab("count") +
  theme_bw() +
  #ggtitle("p-value treat: exclude error regr.") +
  coord_cartesian(ylim = NULL)

####################
#______plots of p values for batch coefficients 
#...(pois vs neg-binom) 
map(list(tb_eo_keep$po1_p_b, tb_eo_keep$nb1_p_b), summary_v)

ggplot(data = tb_eo_keep) +
  geom_density(aes(x = po1_p_b), color = "red") +
  xlab("p-value") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("p-value batch: exclude error regr.") +
  coord_cartesian(ylim = NULL)

ggplot(data = tb_eo_keep) +
  geom_density(aes(x = nb1_p_b), color = "blue") +
  xlab("p-value") +
  ylab("kernel probability density est.") +
  theme_bw() +
  ggtitle("p-value batch: exclude error regr.") +
  coord_cartesian(ylim = NULL)

xran <- 0.25
ggplot(data = tb_eo_keep %>% filter(nb1_p_b < xran)) +
  geom_histogram(aes(x = nb1_p_b), binwidth = 0.01, boundary=0) +
  xlab("p-value") +
  ylab("count") +
  theme_bw() +
  ggtitle("p-value batch: exclude error regr.") +
  coord_cartesian(ylim = NULL)

################################################################################
#______plots of treat coefficients 

#...(pois vs neg-binom) 
map(list(tb_eo_keep$po1_B_t, tb_eo_keep$po2_B_t, tb_eo_keep$nb1_B_t, tb_eo_keep$nb2_B_t), summary_v)

ggplot(data = tb_eo_keep) +
  geom_point(aes(x = po1_B_t, y = nb1_B_t), alpha = 0.1) +
  xlab("treat coeff in poisson") +
  ylab("treat coeff in neg binom") +
  theme_bw() +
  ggtitle("coeff treat: exclude error regr.") +
  coord_cartesian(ylim = NULL)

ggplot(data = tb_eo_keep) +
  geom_histogram(aes(x = nb1_B_t), binwidth = 2, boundary=0) +
  xlab("treat coeff in neg binom") +
  ylab("count") +
  theme_bw() +
  ggtitle("coeff treat: exclude error regr.")

####################
#______plots of infect coefficients 
#...(pois vs neg-binom) 

ggplot(data = tb_eo_keep) +
  geom_point(aes(x = po1_B_i, y = nb1_B_i), alpha = 0.1) +
  xlab("infect coeff in poisson") +
  ylab("infect coeff in neg binom") +
  theme_bw() +
  ggtitle("coeff infect: exclude error regr.") +
  coord_cartesian(ylim = NULL)

ggplot(data = tb_eo_keep) +
  geom_histogram(aes(x = nb1_B_i), binwidth = 2, boundary=0) +
  xlab("infect coeff in neg binom") +
  ylab("count") +
  theme_bw() +
  ggtitle("coeff infect: exclude error regr.")

################################################################################
#Cleanup
xran <- 0
yran <- 0
xyran <- 0
rm(xran, yran, xyran)