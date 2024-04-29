# - Application based on Lalive (2008) data - # 
 
# Libraries 
library(haven)
library(ggplot2)
library(rdrobust)
library(rddtools)
library(rddensity)
library(rdd)
library(stargazer)
library(dplyr)
library(tidyr)

# Load Stata dataset from Lalive (2008)
data_ub = read_dta("/Users/larajoycekaser/Downloads/releaseData/releaseData.dta")

# Check pre-conditions (because we have only received the dataset without further explanations)
# a) all we want units to consider should be between 46 and 54
summary(data_ub$agei)
summary(data_ub$age)
# +

# b) all should have previous experience of at least 0.7
#    as defined by Lalive
summary(data_ub$previous_experience)
# +

# c) exclude steel sector from analysis
# +

# our main focus: people in treated region (where REBP is active)
# this means all with tr == 1
# as well as people that got unemployed during REBP being active
# this means all with period == 1 

# create subset for men, treated region, time period of interest and older than 50
subset_men_period_tr_age50 <- filter(data_ub, female == 0 & tr == 1 & period == 1 & age50 == 1)

# create subset for men, treated region, time period of interest and younger than 50
subset_men_period_tr_age49 <- filter(data_ub, female == 0 & tr == 1 & period == 1 & age50 == 0)

# for a robustness check compare with control regions (where REBP is not active)
# create subset for men, control region, time period of interest and older than 50
subset_men_period_nontr_age50 <- filter(data_ub, female == 0 & tr == 0 & period == 1 & age50 == 1)

# create subset for men, control region, time period of interest and younger than 50
subset_men_period_nontr_age49 <- filter(data_ub, female == 0 & tr == 0 & period == 1 & age50 == 0)


# Create summaries for age 50-53
summary_age50 <- subset_men_period_tr_age50 %>%
  summarize(
    age_mean = mean(age, na.rm = TRUE),
    lwage_ljob_mean = mean(lwage_ljob, na.rm = TRUE),
    single_mean = mean(single, na.rm = TRUE),
    marr_mean = mean(marr, na.rm = TRUE),
    white_collar_mean = mean(white_collar, na.rm = TRUE),
    foreign_mean = mean(foreign, na.rm = TRUE),
    educ_med_mean = mean(educ_med, na.rm = TRUE),
    educ_hi_mean = mean(educ_hi, na.rm = TRUE),
    unemployment_duration_mean = mean(unemployment_duration, na.rm = TRUE),
    previous_experience_mean = mean(previous_experience, na.rm = TRUE),
    rr_mean = mean(rr, na.rm = TRUE),
    n = sum(!is.na(age)) 
  )

# Create summaries for age 50-53
summary_age50_nontr <- subset_men_period_nontr_age50 %>%
  summarize(
    age_mean = mean(age, na.rm = TRUE),
    lwage_ljob_mean = mean(lwage_ljob, na.rm = TRUE),
    single_mean = mean(single, na.rm = TRUE),
    marr_mean = mean(marr, na.rm = TRUE),
    white_collar_mean = mean(white_collar, na.rm = TRUE),
    foreign_mean = mean(foreign, na.rm = TRUE),
    educ_med_mean = mean(educ_med, na.rm = TRUE),
    educ_hi_mean = mean(educ_hi, na.rm = TRUE),
    unemployment_duration_mean = mean(unemployment_duration, na.rm = TRUE),
    previous_experience_mean = mean(previous_experience, na.rm = TRUE),
    rr_mean = mean(rr, na.rm = TRUE),
    n = sum(!is.na(age))  
  )

# Create summaries for age 46-49
summary_age49 <- subset_men_period_tr_age49 %>%
  summarize(
    age_mean = mean(age, na.rm = TRUE),
    lwage_ljob_mean = mean(lwage_ljob, na.rm = TRUE),
    single_mean = mean(single, na.rm = TRUE),
    marr_mean = mean(marr, na.rm = TRUE),
    white_collar_mean = mean(white_collar, na.rm = TRUE),
    foreign_mean = mean(foreign, na.rm = TRUE),
    educ_med_mean = mean(educ_med, na.rm = TRUE),
    educ_hi_mean = mean(educ_hi, na.rm = TRUE),
    unemployment_duration_mean = mean(unemployment_duration, na.rm = TRUE),
    previous_experience_mean = mean(previous_experience, na.rm = TRUE),
    rr_mean = mean(rr, na.rm = TRUE),
    n = sum(!is.na(age))  
  )

# Create summaries for age 46-49
summary_age49_nontr <- subset_men_period_nontr_age49 %>%
  summarize(
    age_mean = mean(age, na.rm = TRUE),
    lwage_ljob_mean = mean(lwage_ljob, na.rm = TRUE),
    single_mean = mean(single, na.rm = TRUE),
    marr_mean = mean(marr, na.rm = TRUE),
    white_collar_mean = mean(white_collar, na.rm = TRUE),
    foreign_mean = mean(foreign, na.rm = TRUE),
    educ_med_mean = mean(educ_med, na.rm = TRUE),
    educ_hi_mean = mean(educ_hi, na.rm = TRUE),
    unemployment_duration_mean = mean(unemployment_duration, na.rm = TRUE),
    previous_experience_mean = mean(previous_experience, na.rm = TRUE),
    rr_mean = mean(rr, na.rm = TRUE),
    n = sum(!is.na(age)) 
  )

# Combine summaries into a single data frame
summary_combined <- bind_rows(
  mutate(summary_age50, age_group = "50-53", group = "treated"),
  mutate(summary_age49, age_group = "46-49", group = "treated"),
  mutate(summary_age50_nontr, age_group = "50-53", group = "control"),
  mutate(summary_age49_nontr, age_group = "46-49", group = "control")
)

# Display the combined summaries (Table 6)
summary_combined


# create datasets for plotting and estimates
# REBP region and period
subset_men_period_tr <- filter(data_ub, female == 0 & tr == 1 & period == 1)
# 9734 obs

# REBP region but period before REBP
subset_men_period_nontr <- filter(data_ub, female == 0 & tr == 0 & period == 1)
# 17572 obs

# non-REBP region but period of REBP
subset_men_nonperiod_tr <- filter(data_ub, female == 0 & tr == 1 & period == 0)
# 9726 obs


# plot unemployment_duration on age (linear regressions on both sides )

# age rounded to quarters (Figure 7)
rdplot(y = subset_men_period_tr$unemployment_duration, x = subset_men_period_tr$agei, c = 50, p = 1,  
       scale = NULL, kernel = "triangular", weights = NULL, 
       covs = NULL, covs_eval = "mean", covs_drop = TRUE, ginv.tol = 1e-20,
       support = NULL, subset = NULL, masspoints = "adjust",
       hide = FALSE, ci = NULL, shade = FALSE, , title = "",
       x.label = "age (years)", y.label = "unemployment duration (weeks)", x.lim = NULL, y.lim = c(10,35), 
       col.dots = NULL, col.lines = NULL) 
ggsave("rd.pdf")



# plot the same for unemployment duration and age 
# for observations in the control group of time before REBP but in 
# the later treated regions

# age rounded to quarters (Figure 8)
rdplot(y = subset_men_nonperiod_tr$unemployment_duration, x = subset_men_nonperiod_tr$agei, c = 50, p = 1,  
       scale = NULL, kernel = "triangular", weights = NULL, 
       covs = NULL, covs_eval = "mean", covs_drop = TRUE, ginv.tol = 1e-20,
       support = NULL, subset = NULL, masspoints = "adjust",
       hide = FALSE, ci = NULL, shade = FALSE, title = "", 
       x.label = "age (years)", y.label = "unemployment duration (weeks)", x.lim = NULL, y.lim = c(10,35), 
       col.dots = NULL, col.lines = NULL) 
ggsave("rd_nonperiod.pdf")

mean(subset_men_nonperiod_tr$unemployment_duration[subset_men_nonperiod_tr$age < 50], na.rm = TRUE)
mean(subset_men_nonperiod_tr$unemployment_duration[subset_men_nonperiod_tr$age >= 50], na.rm = TRUE)
# could endanger the identification strategy !!

# in control region (no REBP)
# age rounded to quarters (Figure 9)
rdplot(y = subset_men_period_nontr$unemployment_duration, x = subset_men_period_nontr$agei, c = 50, p = 1,  
       scale = NULL, kernel = "triangular", weights = NULL, 
       covs = NULL, covs_eval = "mean", covs_drop = TRUE, ginv.tol = 1e-20,
       support = NULL, subset = NULL, masspoints = "adjust",
       hide = FALSE, ci = NULL, shade = FALSE, title = "", 
       x.label = "age (years)", y.label = "unemployment duration (weeks)", x.lim = NULL, y.lim = c(10,35), 
       col.dots = NULL, col.lines = NULL)
ggsave("rd_nontr.pdf")



# further identification checks 
# perform rdd analysis with covariates as the outcome
#  single
h_IK_rdrobust_2014_single <- rdbwselect_2014(subset_men_period_tr$single,subset_men_period_tr$age,c=50,p=1, kernel = "tri", bwselect = "IK")
h_IK_single <- h_IK_rdrobust_2014_single$bws[1]
rdd_result_IK_single <- rdrobust(subset_men_period_tr$single,subset_men_period_tr$age, c = 50, h = h_IK_single, p =1)
tau_hat_IK_single <- rdd_result_IK_single$coef[1]
print(h_IK_single)
print(tau_hat_IK_single)

rdplot(y = subset_men_period_tr$single, x = subset_men_period_tr$age, c = 50, p = 1,  
       scale = NULL, kernel = "triangular", weights = NULL, 
       covs = NULL, covs_eval = "mean", covs_drop = TRUE, ginv.tol = 1e-20,
       support = NULL, subset = NULL, masspoints = "adjust",
       hide = FALSE, ci = NULL, shade = FALSE, title = "", 
       x.label = "age (years)", y.label = "single (share)", x.lim = NULL, y.lim =NULL, 
       col.dots = NULL, col.lines = NULL)
ggsave("rd_single.pdf") # Figure 10(a)

# foreign
h_IK_rdrobust_2014_foreign <- rdbwselect_2014(subset_men_period_tr$foreign,subset_men_period_tr$age,c=50,p=1, kernel = "tri", bwselect = "IK")
h_IK_foreign <- h_IK_rdrobust_2014_foreign$bws[1]
rdd_result_IK_foreign <- rdrobust(subset_men_period_tr$foreign,subset_men_period_tr$age, c = 50, h = h_IK_foreign, p =1)
tau_hat_IK_foreign <- rdd_result_IK_foreign$coef[1]
print(h_IK_foreign)
print(tau_hat_IK_foreign)

rdplot(y = subset_men_period_tr$foreign, x = subset_men_period_tr$age, c = 50, p = 1,  
       scale = NULL, kernel = "triangular", weights = NULL, 
       covs = NULL, covs_eval = "mean", covs_drop = TRUE, ginv.tol = 1e-20,
       support = NULL, subset = NULL, masspoints = "adjust",
       hide = FALSE, ci = NULL, shade = FALSE, title = "", 
       x.label = "age (years)", y.label = "foreign (share)", x.lim = NULL, y.lim = NULL, 
       col.dots = NULL, col.lines = NULL)
ggsave("rd_foreign.pdf") # Figure 10(b)

# white collar 
h_IK_rdrobust_2014_white_collar <- rdbwselect_2014(subset_men_period_tr$white_collar,subset_men_period_tr$age,c=50,p=1, kernel = "tri", bwselect = "IK")
h_IK_white_collar <- h_IK_rdrobust_2014_white_collar$bws[1]
rdd_result_IK_white_collar <- rdrobust(subset_men_period_tr$white_collar,subset_men_period_tr$age, c = 50, h = h_IK_white_collar, p =1)
tau_hat_IK_white_collar <- rdd_result_IK_white_collar$coef[1]
print(h_IK_white_collar)
print(tau_hat_IK_white_collar)

rdplot(y = subset_men_period_tr$white_collar, x = subset_men_period_tr$age, c = 50, p = 1,  
       scale = NULL, kernel = "triangular", weights = NULL, 
       covs = NULL, covs_eval = "mean", covs_drop = TRUE, ginv.tol = 1e-20,
       support = NULL, subset = NULL, masspoints = "adjust",
       hide = FALSE, ci = NULL, shade = FALSE, title = "", 
       x.label = "age (years)", y.label = "white_collar (share)", x.lim = NULL, y.lim = NULL, 
       col.dots = NULL, col.lines = NULL)
ggsave("rd_whitecollar.pdf") # Figure 10(c)

# intermediate education
h_IK_rdrobust_2014_educ_med <- rdbwselect_2014(subset_men_period_tr$educ_med,subset_men_period_tr$age,c=50,p=1, kernel = "tri", bwselect = "IK")
h_IK_educ_med <- h_IK_rdrobust_2014_educ_med$bws[1]
rdd_result_IK_educ_med <- rdrobust(subset_men_period_tr$educ_med,subset_men_period_tr$age, c = 50, h = h_IK_educ_med, p =1)
tau_hat_IK_educ_med <- rdd_result_IK_educ_med$coef[1]
print(h_IK_educ_med)
print(tau_hat_IK_educ_med)

rdplot(y = subset_men_period_tr$educ_med, x = subset_men_period_tr$age, c = 50, p = 1,  
       scale = NULL, kernel = "triangular", weights = NULL, 
       covs = NULL, covs_eval = "mean", covs_drop = TRUE, ginv.tol = 1e-20,
       support = NULL, subset = NULL, masspoints = "adjust",
       hide = FALSE, ci = NULL, shade = FALSE, title = "", 
       x.label = "age (years)", y.label = "educ_med (share)", x.lim = NULL, y.lim = NULL, 
       col.dots = NULL, col.lines = NULL)
ggsave("rd_educmed.pdf") # Figure 10(d)

# high education
h_IK_rdrobust_2014_educ_hi <- rdbwselect_2014(subset_men_period_tr$educ_hi,subset_men_period_tr$age,c=50,p=1, kernel = "tri", bwselect = "IK")
h_IK_educ_hi <- h_IK_rdrobust_2014_educ_hi$bws[1]
rdd_result_IK_educ_hi <- rdrobust(subset_men_period_tr$educ_hi,subset_men_period_tr$age, c = 50, h = h_IK_educ_hi, p =1)
tau_hat_IK_educ_hi <- rdd_result_IK_educ_hi$coef[1]
print(h_IK_educ_hi)
print(tau_hat_IK_educ_hi)

rdplot(y = subset_men_period_tr$educ_hi, x = subset_men_period_tr$age, c = 50, p = 1,  
       scale = NULL, kernel = "triangular", weights = NULL, 
       covs = NULL, covs_eval = "mean", covs_drop = TRUE, ginv.tol = 1e-20,
       support = NULL, subset = NULL, masspoints = "adjust",
       hide = FALSE, ci = NULL, shade = FALSE, title = "", 
       x.label = "age (years)", y.label = "educ_hi (share)", x.lim = NULL, y.lim = NULL, 
       col.dots = NULL, col.lines = NULL)
ggsave("rd_educhi.pdf") # Figure 10(e)

# log(wage last job)
h_IK_rdrobust_2014_lwage_ljob <- rdbwselect_2014(subset_men_period_tr$lwage_ljob,subset_men_period_tr$age,c=50,p=1, kernel = "tri", bwselect = "IK")
h_IK_lwage_ljob <- h_IK_rdrobust_2014_lwage_ljob$bws[1]
rdd_result_IK_lwage_ljob <- rdrobust(subset_men_period_tr$lwage_ljob,subset_men_period_tr$age, c = 50, h = h_IK_lwage_ljob, p =1)
tau_hat_IK_lwage_ljob <- rdd_result_IK_lwage_ljob$coef[1]
print(h_IK_lwage_ljob)
print(tau_hat_IK_lwage_ljob)

rdplot(y = subset_men_period_tr$lwage_ljob, x = subset_men_period_tr$age, c = 50, p = 1,  
       scale = NULL, kernel = "triangular", weights = NULL, 
       covs = NULL, covs_eval = "mean", covs_drop = TRUE, ginv.tol = 1e-20,
       support = NULL, subset = NULL, masspoints = "adjust",
       hide = FALSE, ci = NULL, shade = FALSE, title = "", 
       x.label = "age (years)", y.label = "lwage_ljob (share)", x.lim = NULL, y.lim = NULL, 
       col.dots = NULL, col.lines = NULL)
ggsave("rd_lwageljob.pdf") # Figure 10(f)




# is the treatment-determining variable (age) continuous at c=50?
# perform McCrary Test (Figure 11)

# REBP treated (a)
DCdensity(subset_men_period_tr$age, cutpoint = 50, bin = NULL, bw = NULL, verbose = FALSE,
          plot = TRUE, ext.out = FALSE, htest = TRUE)

# control (b)
DCdensity(subset_men_period_nontr$age, cutpoint = 50, bin = NULL, bw = NULL, verbose = FALSE,
          plot = TRUE, ext.out = FALSE, htest = TRUE)

# control (c)
DCdensity(subset_men_nonperiod_tr$age, cutpoint = 50, bin = NULL, bw = NULL, verbose = FALSE,
          plot = TRUE, ext.out = FALSE, htest = TRUE)

# all results are significant, consider this in analysis and interpretation


# relative inflow = ratio of the density of age in treated region to density of age in the control region 
# -> should be around 1 so no concentration around the threshold -> no reaction?
# density = ratio of number of unemployed observed in each age bracket relative to total number of spells in treated regions 
# then ratio between these two ratios 

# subset of only men 
subset_men <- filter(data_ub, female == 0)

# calculate density in treated and control region 
density_tr <- table(subset_men[subset_men$tr == 1, "age"]) / sum(subset_men$tr == 1)
density_control <- table(subset_men[subset_men$tr == 0, "age"]) / sum(subset_men$tr == 0)

# calculate relative inflow
relative_inflow <- density_tr / density_control
df_rel_inflow <- data.frame(age = as.numeric(names(relative_inflow)),
                            relative_inflow = unname(relative_inflow))
df_rel_inflow

# check for missing values in df_rel_inflow
missing_values <- sum(is.na(df_rel_inflow))

# print the number of missing values
print(missing_values)
# no missings 

# Plot with rdplot(Figure 12)
rdplot(y = df_rel_inflow$relative_inflow.Freq, x = df_rel_inflow$age, c = 50, p = 1,  
       scale = NULL, kernel = "triangular", weights = NULL, 
       covs = NULL, covs_eval = "mean", covs_drop = TRUE, ginv.tol = 1e-20,
       support = NULL, subset = NULL, masspoints = "adjust",
       hide = FALSE, ci = NULL, shade = FALSE, title = "", 
       x.label = "age (years)",
       y.label = "relative inflow", x.lim = NULL, y.lim = c(0,2), 
       col.dots = NULL, col.lines = NULL)

# estimate treatment effect at threshold for y = relative inflow 
h_IK_select_check <- rdbwselect_2014(df_rel_inflow$relative_inflow.Freq, df_rel_inflow$age,c=50,p=1, kernel = "tri", bwselect = "IK")
h_IK_check <- h_IK_select_check$bws[1]
check <- rdrobust(df_rel_inflow$relative_inflow.Freq, df_rel_inflow$age, c = 50, h = h_IK_check, p =1)



# rdd analysis 

# IK bandwidth
h_IK_rdrobust_2014_app <- rdbwselect_2014(subset_men_period_tr$unemployment_duration,subset_men_period_tr$age,c=50,p=1, kernel = "tri", bwselect = "IK")
h_IK_app <- h_IK_rdrobust_2014_app$bws[1]

# robustness checks bandwidths 
h_alt1 <- 1
h_alt2 <- 2
h_alt3 <- 3

# CV bandwidth - not feasible 
h_CV_rdrobust_2014_app <- rdbwselect_2014(subset_men_period_tr$unemployment_duration,subset_men_period_tr$age,c=50,p=1, kernel = "tri", bwselect = "CV")
h_CV_app <- h_CV_rdrobust_2014_app$bws[1]
rdd_result_CV_app <- rdrobust(subset_men_period_tr$unemployment_duration,subset_men_period_tr$age, c = 50, h = h_CV_app, p =1)
# does not work ? Error in qr.default(XX_CV_l * sqrt(w_CV_l), tol = 1e-10) : 
# NA/NaN/Inf in foreign function call (arg 1)
anyNA(subset_men_period_tr$unemployment_duration)
anyNA(subset_men_period_tr$age)
# no missings 


# Estimation of ATE using RDD

# ...with IK bandwidth (Table 7, Column (1))
rdd_result_IK_app <- rdrobust(subset_men_period_tr$unemployment_duration,subset_men_period_tr$age, c = 50, h = h_IK_app, p =1)
tau_hat_IK_app <- rdd_result_IK_app$coef[1]
print(rdd_result_IK_app)
print(tau_hat_IK_app)

#plot only within interval of bandwidth (Figure 13)
rdplot(y = subset_men_nonperiod_tr$unemployment_duration, x = subset_men_nonperiod_tr$age, c = 50, p = 1,  
       scale = NULL, kernel = "triangular", h=h_IK_app, weights = NULL, 
       covs = NULL, covs_eval = "mean", covs_drop = TRUE, ginv.tol = 1e-20,
       support = NULL, subset = NULL, masspoints = "adjust",
       hide = FALSE, ci = NULL, shade = FALSE, title = "", 
       x.label = "age (years)", y.label = "unemployment duration (weeks)", x.lim = c(50-rdd_result_IK_app$bws[1],50+rdd_result_IK_app$bws[1]),, y.lim = NULL, 
       col.dots = NULL, col.lines = NULL)


# ...with bandwidth 1 (Table 7, Column (2))
rdd_result_alt1 <- rdrobust(subset_men_period_tr$unemployment_duration,subset_men_period_tr$age, c = 50, h = h_alt1, p =1 )
tau_hat_alt1 <- rdd_result_alt1$coef[1]


# ...with bandwidth 2 (Table 7, Column (3))
rdd_result_alt2 <- rdrobust(subset_men_period_tr$unemployment_duration,subset_men_period_tr$age, c = 50, h = h_alt2, p =1 )
tau_hat_alt2 <- rdd_result_alt2$coef[1]


# ...with bandwidth 3 (Table 7, Column (4))
rdd_result_alt3 <- rdrobust(subset_men_period_tr$unemployment_duration,subset_men_period_tr$age, c = 50, h = h_alt3, p =1 )
tau_hat_alt3 <- rdd_result_alt3$coef[1]


# estimates with three significant control variables (use IK bandwidth) (Table 7, Column (5))
rdd_result_ik_control_3 <- rdrobust(subset_men_period_tr$unemployment_duration,subset_men_period_tr$age, c = 50, h = h_IK_app, p =1, 
                                    covs = cbind(subset_men_period_tr$single, subset_men_period_tr$white_collar, subset_men_period_tr$educ_med))


# estimates with all presented control variables (use IK bandwidth) (Table 7, Column (6))
rdd_result_ik_control <- rdrobust(subset_men_period_tr$unemployment_duration,subset_men_period_tr$age, c = 50, h = h_IK_app, p =1, 
                                  covs = cbind(subset_men_period_tr$single, subset_men_period_tr$foreign, subset_men_period_tr$white_collar, subset_men_period_tr$lwage_ljob, subset_men_period_tr$educ_hi, subset_men_period_tr$educ_med))



# test at alternative cutoffs 

# use subsample subset_men_period_tr_age49 to test at median
median1 <- median(subset_men_period_tr_age49$age)
h_ik_med1 <- rdbwselect_2014(subset_men_period_tr_age49$unemployment_duration,subset_men_period_tr_age49$age,c=median1,p=1, kernel = "tri", bwselect = "IK")
rd_ik_med1 <- rdrobust(subset_men_period_tr_age49$unemployment_duration,subset_men_period_tr_age49$age, c = median1, h = h_ik_med1$bws[1], p =1)
tau_ik_med1 <- rd_ik_med1$coef[1]

# use subsample subset_men_period_tr_age50 to test at median
median2 <- median(subset_men_period_tr_age50$age)
h_ik_med2 <- rdbwselect_2014(subset_men_period_tr_age50$unemployment_duration,subset_men_period_tr_age50$age,c=median2,p=1, kernel = "tri", bwselect = "IK")
rd_ik_med2 <- rdrobust(subset_men_period_tr_age50$unemployment_duration,subset_men_period_tr_age50$age, c = median2, h = h_ik_med2$bws[1], p =1)
tau_ik_med2 <- rd_ik_med2$coef[1]









## Further analysis not in the thesis ##

# test clustering for local linear regression 
rdd_result_ik_cluster <- rdrobust(subset_men_period_tr$unemployment_duration,subset_men_period_tr$age, c = 50, h = h_IK_app, p =1, cluster = subset_men_period_tr$age)

# Replication of Lalive (2008) results (parametric regression)
# Define the control variables
controls <- c("marr", "single", "educ_med", "educ_hi", "foreign", "rr", "lwage_ljob",
              "previous_experience", "white_collar", "landw", "versorg", "nahrung",
              "textil", "holzind", "elmasch", "andfabr", "bau", "gasthand", 
              "verkehr", "dienstl")

# Replicate result in Table 2 Column (2)
model1 <- lm(unemployment_duration ~ age50 + dage_1 + age50_dage_1, 
             data = data_ub,
             subset = female == 0 & period == 1 & tr == 1)
summary(model1)

# Replicate result in Table 2 Column (6)
model2 <- lm(unemployment_duration ~ age50 + dage_1 + age50_dage_1 + marr + single + educ_med + educ_hi + foreign + rr + lwage_ljob +
               previous_experience + white_collar + landw + versorg + nahrung +
               textil + holzind + elmasch + andfabr + bau + gasthand + 
               verkehr + dienstl, 
             data = data_ub,
             subset = female == 0 & period == 1 & tr == 1)
summary(model2)

