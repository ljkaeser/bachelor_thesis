# - Monte Carlo Simulations - # 
## - low curvature setting - ##

# libraries 
library(ggplot2)
library(rdrobust)
library(rddtools)
library(rddensity)
library(rdd)
library(dplyr)
library(stargazer)

# (1) Define design setup #

# define cutoff point
c <- 0 

# define data generating functions 
m_pos <- function(x) 1/8*(x+2)^3 + 3*x 
m_neg <- function(x) 1/27*(x-3)^3 + 1.5

# generate x values to plot dgf
x_pos <- seq(0, 1, length.out = 100)  # Generate x values for m_pos (x >= 0)
x_neg <- seq(-1, 0, length.out = 100) # Generate x values for m_neg (x < 0)

# calculate corresponding y values
y_pos <- m_pos(x_pos)
y_neg <- m_neg(x_neg)

# create a data frame for ggplot
data_pos <- data.frame(x = x_pos, y = y_pos)
data_neg <- data.frame(x = x_neg, y = y_neg)

# plot using ggplot (Figure 2 (a))
ggplot() +
  geom_line(data = data_pos, aes(x = x, y = y), color = "blue", size = 1) +
  geom_line(data = data_neg, aes(x = x, y = y), color = "red", size = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey", size = 1) +
  labs(title = "",
       x = "x",
       y = "y", size = 8) +
  scale_y_continuous(limits = c(-5, 5)) +  # Set the y-axis limits
  theme_minimal() + theme(#panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"))  # Make axis lines black

# save the plot
ggsave("dgf.pdf")

# define the real treatment effect we try to estimate in the following 
true_ate = m_pos(0)-m_neg(0)
print(true_ate)

###############################################################################

# (2) Run 1000 simulations for sample size 100 #

# define sample size and number of repeats 
n <- 100
n_reps <- 1000

# generate vectors to store estimated treatment effects
ate_100_ik <- numeric(n_reps)
ate_100_cv <- numeric(n_reps)
ate_100_2ik <- numeric(n_reps)
ate_100_hik <- numeric(n_reps)

# generate vectors to store bandwidths
h_100_ik <- numeric(n_reps)
h_100_cv <- numeric(n_reps)
h_100_2ik <- numeric(n_reps)
h_100_hik <- numeric(n_reps)

# generate vectors to store sample size within bandwidth interval
n_100_h_ik <- numeric(n_reps)
n_100_h_2ik <- numeric(n_reps)
n_100_h_hik <- numeric(n_reps)
n_100_h_cv <- numeric(n_reps)

# loop for replicates
for (i in 1:n_reps) {
  # set seed for reproducibility
  set.seed(i)  
  
  # generate new data for each sample
  x <- rnorm(n) # generate observations from standard normal distribution (mean 0, sd 1)
  y <- ifelse(x < c, m_neg(x), m_pos(x) ) + rnorm(n, sd = 0.5) # add normal error term with variance 1/4
  w <- ifelse(x<c,0,1) # variables indicates treatment
  
  # IK bandwidth
  h_IK_rdrobust_2014 <- rdbwselect_2014(y,x,c=0,p=1, kernel = "tri", bwselect = "IK")
  h_IK <- h_IK_rdrobust_2014$bws[1]
  h_100_ik[i] <- h_IK 
  
  # robustness check bandwidth
  h_2IK <- 2*h_IK
  h_100_2ik[i] <- h_2IK 
  
  # robustness check bandwidth
  h_hIK <- 0.5*h_IK
  h_100_hik[i] <- h_hIK
  
  # CV bandwidth
  h_CV_rdrobust_2014 <- rdbwselect_2014(y=y, x = x, bwselect = "CV")
  h_CV <- h_CV_rdrobust_2014$bws[1]
  h_100_cv[i] <- h_CV
  
  # calculate sample size within interval 
  n_100_h_ik[i] <- sum(x >= (c - h_IK) & x <= (c + h_IK))
  n_100_h_2ik[i] <- sum(x >= (c - h_2IK) & x <= (c + h_2IK))
  n_100_h_hik[i] <- sum(x >= (c - h_hIK) & x <= (c + h_hIK))
  n_100_h_cv[i] <- sum(x >= (c - h_CV) & x <= (c + h_CV))
  
  # estimation of ATE using IK bandwidth
  rdd_result_IK <- rdrobust(y, x, c = 0, h = h_IK)
  ate_100_ik[i] <- rdd_result_IK$coef[1]
  
  # estimation of ATE using 2IK bandwidth
  rdd_result_2IK <- rdrobust(y, x, c = 0, h = h_2IK)
  ate_100_2ik[i] <- rdd_result_2IK$coef[1]
  
  # estimation of ATE using hIK bandwidth
  #rdd_result_hIK <- rdrobust(y, x, c = 0, h = h_hIK)
  #ate_100_hik[i] <- rdd_result_hIK$coef[1]
  # !! sample size is too low, not possible for n=100
  
  # estimation of ATE using CV bandwidth
  rdd_result_CV <- rdrobust(y, x, c = 0, h = h_CV)
  ate_100_cv[i] <- rdd_result_CV$coef[1]
}

# average of 1000 derived estimates 
mean_ate_100_ik <- mean(ate_100_ik)
mean_ate_100_2ik <- mean(ate_100_2ik)
mean_ate_100_hik <- mean(ate_100_hik)
mean_ate_100_cv <- mean(ate_100_cv)

# bias of 1000 derived estimates 
bias_ate_100_ik <- mean_ate_100_ik - true_ate
bias_ate_100_2ik <- mean_ate_100_2ik - true_ate
bias_ate_100_hik <- mean_ate_100_hik - true_ate
bias_ate_100_cv <- mean_ate_100_cv - true_ate

# variance of 1000 derived estimates 
var_ate_100_ik <- var(ate_100_ik)
var_ate_100_2ik <- var(ate_100_2ik)
var_ate_100_hik <- var(ate_100_hik)
var_ate_100_cv <- var(ate_100_cv)

# mse of 1000 derived estimates 
mse_ate_100_ik <- bias_ate_100_ik^2 + var_ate_100_ik
mse_ate_100_2ik <- bias_ate_100_2ik^2 + var_ate_100_2ik
mse_ate_100_hik <- bias_ate_100_hik^2 + var_ate_100_hik
mse_ate_100_cv <- bias_ate_100_cv^2 + var_ate_100_cv

# mean of 1000 calculated bandwidths
mean_h_100_ik <- mean(h_100_ik)
mean_h_100_2ik <- mean(h_100_2ik)
mean_h_100_hik <- mean(h_100_hik)
mean_h_100_cv <- mean(h_100_cv)

# mean of effectively used sample sizes
mean_n_h_100_ik <- mean(n_100_h_ik)
mean_n_h_100_2ik <- mean(n_100_h_2ik)
mean_n_h_100_hik <- mean(n_100_h_hik)
mean_n_h_100_cv <- mean(n_100_h_cv)

# visualization of MSE in dependence of bandwidth:
# generate data set
data100 = data.frame(sample_size = c(100, 100, 100),
                     sample_within_h = c(mean_n_h_100_ik,mean_n_h_100_2ik, mean_n_h_100_cv),
                     procedure = c("IK", "2IK", "CV" ),
                     bandwidth = c(mean_h_100_ik, mean_h_100_2ik, mean_h_100_cv),
                     bias = c(bias_ate_100_ik, bias_ate_100_2ik, bias_ate_100_cv),
                     var = c(var_ate_100_ik, var_ate_100_2ik, var_ate_100_cv),
                     mse = c(mse_ate_100_ik, mse_ate_100_2ik, mse_ate_100_cv))

# plot (Figure 3 (a))
ggplot(data100, aes(x = bandwidth)) +
  geom_line(aes(y = mse, color = "MSE"), size = 1) +
  geom_line(aes(y = bias, color = "Bias"), size = 1) +
  geom_line(aes(y = var, color = "Variance"), size = 1) +
  geom_text(aes(label = procedure, y=0), vjust = -0.5, hjust = 0.5, size = 5) +  # Add labels for 
  geom_vline(xintercept = data100$bandwidth, linetype = "dotted", color= "grey") +  # Add vertical lines at bandwidth points
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       x = "bandwidth", y="values"
  ) +
  theme_minimal() + theme(legend.position = "none",  # Turn off legend
                          panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          axis.line = element_line(color = "black"))  # Make axis lines black

# save the plot
ggsave("bvm_h_100.pdf")

# save data100 in latex table (Table 2 (a))
latex_100 <- stargazer(data100, summary = FALSE)
writeLines(latex_100, file.path("latex_100.txt"))

################################################################################
# (3) Run 1000 simulations for sample size 200 #

# define sample size and number of repeats 
n <- 200
n_reps <- 1000

# vector to store estimated treatment effects
ate_200_ik <- numeric(n_reps)
ate_200_cv <- numeric(n_reps)
ate_200_2ik <- numeric(n_reps)
ate_200_hik <- numeric(n_reps)

# vector to store bandwidths
h_200_ik <- numeric(n_reps)
h_200_cv <- numeric(n_reps)
h_200_2ik <- numeric(n_reps)
h_200_hik <- numeric(n_reps)

# vector to store effectively used sample size
n_200_h_ik <- numeric(n_reps)
n_200_h_2ik <- numeric(n_reps)
n_200_h_hik <- numeric(n_reps)
n_200_h_cv <- numeric(n_reps)

# loop 
for (i in 1:n_reps) {
  # set seed for reproducibility
  set.seed(i)  
  
  # generate new data for each sample
  x <- rnorm(n) # generate observations from standard normal distribution (mean 0, sd 1)
  y <- ifelse(x < c, m_neg(x), m_pos(x) ) + rnorm(n, sd = 0.5) # add normal error term with variance 1/4
  w <- ifelse(x<c,0,1) # variable indicates treatment
  
  # calculate and store IK bandwidth
  h_IK_rdrobust_2014 <- rdbwselect_2014(y,x,c=0,p=1, kernel = "tri", bwselect = "IK")
  h_IK <- h_IK_rdrobust_2014$bws[1]
  h_200_ik[i] <- h_IK 
  
  # calculate and store 2IK bandwidth
  h_2IK <- 2*h_IK
  h_200_2ik[i] <- h_2IK 
  
  # calculate and store hIK bandwidth
  h_hIK <- 0.5*h_IK
  h_200_hik[i] <- h_hIK
  
  # calculate and store CV bandwidth
  h_CV_rdrobust_2014 <- rdbwselect_2014(y=y, x = x, bwselect = "CV")
  h_CV <- h_CV_rdrobust_2014$bws[1]
  h_200_cv[i] <- h_CV
  
  # calculate and store effectively used sample sizes
  n_200_h_ik[i] <- sum(x >= (c - h_IK) & x <= (c + h_IK))
  n_200_h_2ik[i] <- sum(x >= (c - h_2IK) & x <= (c + h_2IK))
  n_200_h_hik[i] <- sum(x >= (c - h_hIK) & x <= (c + h_hIK))
  n_200_h_cv[i] <- sum(x >= (c - h_CV) & x <= (c + h_CV))
  
  # estimation of ATE using IK bandwidth
  rdd_result_IK <- rdrobust(y, x, c = 0, h = h_IK)
  ate_200_ik[i] <- rdd_result_IK$coef[1]
  
  # estimation of ATE using 2IK bandwidth
  rdd_result_2IK <- rdrobust(y, x, c = 0, h = h_2IK)
  ate_200_2ik[i] <- rdd_result_2IK$coef[1]
  
  # estimation of ATE using hIK bandwidth
  rdd_result_hIK <- rdrobust(y, x, c = 0, h = h_hIK)
  ate_200_hik[i] <- rdd_result_hIK$coef[1]
  
  # estimation of ATE using CV bandwidth
  rdd_result_CV <- rdrobust(y, x, c = 0, h = h_CV)
  ate_200_cv[i] <- rdd_result_CV$coef[1]
}

# average of 1000 derived estimates 
mean_ate_200_ik <- mean(ate_200_ik)
mean_ate_200_2ik <- mean(ate_200_2ik)
mean_ate_200_hik <- mean(ate_200_hik)
mean_ate_200_cv <- mean(ate_200_cv)

# bias of 1000 derived estimates 
bias_ate_200_ik <- mean_ate_200_ik - true_ate
bias_ate_200_2ik <- mean_ate_200_2ik - true_ate
bias_ate_200_hik <- mean_ate_200_hik - true_ate
bias_ate_200_cv <- mean_ate_200_cv - true_ate

# variance of 1000 derived estimates 
var_ate_200_ik <- var(ate_200_ik)
var_ate_200_2ik <- var(ate_200_2ik)
var_ate_200_hik <- var(ate_200_hik)
var_ate_200_cv <- var(ate_200_cv)

# mse of 1000 derived estimates 
mse_ate_200_ik <- bias_ate_200_ik^2 + var_ate_200_ik
mse_ate_200_2ik <- bias_ate_200_2ik^2 + var_ate_200_2ik
mse_ate_200_hik <- bias_ate_200_hik^2 + var_ate_200_hik
mse_ate_200_cv <- bias_ate_200_cv^2 + var_ate_200_cv

# average of bandwidths
mean_h_200_ik <- mean(h_200_ik)
mean_h_200_2ik <- mean(h_200_2ik)
mean_h_200_hik <- mean(h_200_hik)
mean_h_200_cv <- mean(h_200_cv)

# average of effectively used sample sizes
mean_n_h_200_ik <- mean(n_200_h_ik)
mean_n_h_200_2ik <- mean(n_200_h_2ik)
mean_n_h_200_hik <- mean(n_200_h_hik)
mean_n_h_200_cv <- mean(n_200_h_cv)

# visualization of MSE in dependence of bandwidth 
data200 = data.frame(sample_size = c(200, 200, 200, 200),
                     sample_within_h = c(mean_n_h_200_ik,mean_n_h_200_2ik, mean_n_h_200_hik, mean_n_h_200_cv),
                     procedure = c("IK", "2IK", "hIK", "CV" ),
                     bandwidth = c(mean_h_200_ik, mean_h_200_2ik, mean_h_200_hik, mean_h_200_cv),
                     bias = c(bias_ate_200_ik, bias_ate_200_2ik, bias_ate_200_hik, bias_ate_200_cv),
                     var = c(var_ate_200_ik, var_ate_200_2ik,var_ate_200_hik, var_ate_200_cv),
                     mse = c(mse_ate_200_ik, mse_ate_200_2ik, mse_ate_200_hik, mse_ate_200_cv))

# plot (Figure 3(b))
ggplot(data200, aes(x = bandwidth)) +
  geom_line(aes(y = mse, color = "MSE"), size =1 ) +
  geom_line(aes(y = bias, color = "Bias"), size = 1) +
  geom_line(aes(y = var, color = "Variance"), size = 1) +
  geom_text(aes(label = procedure, y=0), vjust = -0.5, hjust = 0.5, size = 5) +  # Add labels for 
  geom_vline(xintercept = data200$bandwidth, linetype = "dotted", color= "grey") +  # Add vertical lines at bandwidth points
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       x = "bandwidth", y="values"
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          axis.line = element_line(color = "black"),
                          legend.position = "none" )# Turn off legend)  # Make axis lines black

# save the plot
ggsave("bvm_h_200.pdf")


# save data200 in Latex Table (Table 2(b))
latex_200 <- stargazer(data200, summary = FALSE)
writeLines(latex_200, file.path("latex_200.txt"))

################################################################################
# (4) Run 1000 simulations for sample size 500 #

# define sample size and number of repeats 
n <- 500
n_reps <- 1000

# vector to store estimated treatment effects
ate_500_ik <- numeric(n_reps)
ate_500_cv <- numeric(n_reps)
ate_500_2ik <- numeric(n_reps)
ate_500_hik <- numeric(n_reps)

# vector to store bandwidths
h_500_ik <- numeric(n_reps)
h_500_cv <- numeric(n_reps)
h_500_2ik <- numeric(n_reps)
h_500_hik <- numeric(n_reps)

# vector to store effective sample sizes
n_500_h_ik <- numeric(n_reps)
n_500_h_2ik <- numeric(n_reps)
n_500_h_hik <- numeric(n_reps)
n_500_h_cv <- numeric(n_reps)

# loop 
for (i in 1:n_reps) {
  # set seed for reproducibility
  set.seed(i) 
  
  # generate new data for each sample
  x <- rnorm(n) # generate observations from standard normal distribution (mean 0, sd 1)
  y <- ifelse(x < c, m_neg(x), m_pos(x) ) + rnorm(n, sd = 0.5) # add normal error term with variance 1/4
  w <- ifelse(x<c,0,1) # indicates treatment
  
  # calculate and store IK bandwidth
  h_IK_rdrobust_2014 <- rdbwselect_2014(y,x,c=0,p=1, kernel = "tri", bwselect = "IK")
  h_IK <- h_IK_rdrobust_2014$bws[1]
  h_500_ik[i] <- h_IK 
  
  # calculate and store 2IK bandwidth
  h_2IK <- 2*h_IK
  h_500_2ik[i] <- h_2IK 
  
  # calculate and store hIK bandwidth
  h_hIK <- 0.5*h_IK
  h_500_hik[i] <- h_hIK
  
  # calculate and store CV bandwidth
  h_CV_rdrobust_2014 <- rdbwselect_2014(y=y, x = x, bwselect = "CV")
  h_CV <- h_CV_rdrobust_2014$bws[1]
  h_500_cv[i] <- h_CV
  
  # calculate and store effective sample sizes
  n_500_h_ik[i] <- sum(x >= (c - h_IK) & x <= (c + h_IK))
  n_500_h_2ik[i] <- sum(x >= (c - h_2IK) & x <= (c + h_2IK))
  n_500_h_hik[i] <- sum(x >= (c - h_hIK) & x <= (c + h_hIK))
  n_500_h_cv[i] <- sum(x >= (c - h_CV) & x <= (c + h_CV))
  
  # estimation of ATE using IK bandwidth
  rdd_result_IK <- rdrobust(y, x, c = 0, h = h_IK)
  ate_500_ik[i] <- rdd_result_IK$coef[1]
  
  # estimation of ATE using 2IK bandwidth
  rdd_result_2IK <- rdrobust(y, x, c = 0, h = h_2IK)
  ate_500_2ik[i] <- rdd_result_2IK$coef[1]
  
  # estimation of ATE using hIK bandwidth
  rdd_result_hIK <- rdrobust(y, x, c = 0, h = h_hIK)
  ate_500_hik[i] <- rdd_result_hIK$coef[1]
  
  # estimation of ATE using CV bandwidth
  rdd_result_CV <- rdrobust(y, x, c = 0, h = h_CV)
  ate_500_cv[i] <- rdd_result_CV$coef[1]
}

# average of 1000 derived estimates 
mean_ate_500_ik <- mean(ate_500_ik)
mean_ate_500_2ik <- mean(ate_500_2ik)
mean_ate_500_hik <- mean(ate_500_hik)
mean_ate_500_cv <- mean(ate_500_cv)

# bias 
bias_ate_500_ik <- mean_ate_500_ik - true_ate
bias_ate_500_2ik <- mean_ate_500_2ik - true_ate
bias_ate_500_hik <- mean_ate_500_hik - true_ate
bias_ate_500_cv <- mean_ate_500_cv - true_ate

# variance
var_ate_500_ik <- var(ate_500_ik)
var_ate_500_2ik <- var(ate_500_2ik)
var_ate_500_hik <- var(ate_500_hik)
var_ate_500_cv <- var(ate_500_cv)

# mse
mse_ate_500_ik <- bias_ate_500_ik^2 + var_ate_500_ik
mse_ate_500_2ik <- bias_ate_500_2ik^2 + var_ate_500_2ik
mse_ate_500_hik <- bias_ate_500_hik^2 + var_ate_500_hik
mse_ate_500_cv <- bias_ate_500_cv^2 + var_ate_500_cv

# mean of bandwidth
mean_h_500_ik <- mean(h_500_ik)
mean_h_500_2ik <- mean(h_500_2ik)
mean_h_500_hik <- mean(h_500_hik)
mean_h_500_cv <- mean(h_500_cv)

# mean of effective sample size 
mean_n_h_500_ik <- mean(n_500_h_ik)
mean_n_h_500_2ik <- mean(n_500_h_2ik)
mean_n_h_500_hik <- mean(n_500_h_hik)
mean_n_h_500_cv <- mean(n_500_h_cv)

# visualization of MSE in dependence of bandwidth 
data500 = data.frame(sample_size = c(500, 500, 500, 500),
                     sample_within_h = c(mean_n_h_500_ik,mean_n_h_500_2ik, mean_n_h_500_hik, mean_n_h_500_cv),
                     procedure = c("IK", "2IK", "hIK", "CV" ),
                     bandwidth = c(mean_h_500_ik, mean_h_500_2ik, mean_h_500_hik, mean_h_500_cv),
                     bias = c(bias_ate_500_ik, bias_ate_500_2ik, bias_ate_500_hik, bias_ate_500_cv),
                     var = c(var_ate_500_ik, var_ate_500_2ik,var_ate_500_hik, var_ate_500_cv),
                     mse = c(mse_ate_500_ik, mse_ate_500_2ik, mse_ate_500_hik, mse_ate_500_cv))

# plot (Figure 3(c))
ggplot(data500, aes(x = bandwidth)) +
  geom_line(aes(y = mse, color = "MSE"), size = 1) +
  geom_line(aes(y = bias, color = "Bias"), size = 1) +
  geom_line(aes(y = var, color = "Variance"), size = 1) +
  geom_text(aes(label = procedure, y=0), vjust = -0.5, hjust = 0.5, size = 5) +  # Add labels for 
  geom_vline(xintercept = data500$bandwidth, linetype = "dotted", color= "grey") +  # Add vertical lines at bandwidth points
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       x = "bandwidth", y="values"
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          legend.position = "none",  # Turn off legend
                          axis.line = element_line(color = "black"))  # Make axis lines black

# save the plot
ggsave("bvm_h_500.pdf")


# save data500 in latex table (Table 2(c))
latex_500 <- stargazer(data500, summary = FALSE)
writeLines(latex_500, file.path("latex_500.txt"))

#########################################################################
# (5) Run 1000 simulations for sample size 1000 #

# define sample size and number of repeats 
n <- 1000
n_reps <- 1000

# Vector to store estimated treatment effects
ate_1000_ik <- numeric(n_reps)
ate_1000_cv <- numeric(n_reps)
h_1000_ik <- numeric(n_reps)
h_1000_cv <- numeric(n_reps)
ate_1000_2ik <- numeric(n_reps)
ate_1000_hik <- numeric(n_reps)
h_1000_2ik <- numeric(n_reps)
h_1000_hik <- numeric(n_reps)

n_1000_h_ik <- numeric(n_reps)
n_1000_h_2ik <- numeric(n_reps)
n_1000_h_hik <- numeric(n_reps)
n_1000_h_cv <- numeric(n_reps)

# Loop over each sample
for (i in 1:n_reps) {
  set.seed(i)  # Set seed for reproducibility
  
  # Generate new data for each sample
  x <- rnorm(n) # generate observations from standard normal distribution (mean 0, sd 1)
  y <- ifelse(x < c, m_neg(x), m_pos(x) ) + rnorm(n, sd = 0.5) # add normal error term
  w <- ifelse(x<c,0,1) 
  
  
  h_IK_rdrobust_2014 <- rdbwselect_2014(y,x,c=0,p=1, kernel = "tri", bwselect = "IK")
  h_IK <- h_IK_rdrobust_2014$bws[1]
  h_1000_ik[i] <- h_IK 
  
  h_2IK <- 2*h_IK
  h_1000_2ik[i] <- h_2IK 
  
  h_hIK <- 0.5*h_IK
  h_1000_hik[i] <- h_hIK
  
  h_CV_rdrobust_2014 <- rdbwselect_2014(y=y, x = x, bwselect = "CV")
  h_CV <- h_CV_rdrobust_2014$bws[1]
  h_1000_cv[i] <- h_CV
  
  n_1000_h_ik[i] <- sum(x >= (c - h_IK) & x <= (c + h_IK))
  n_1000_h_2ik[i] <- sum(x >= (c - h_2IK) & x <= (c + h_2IK))
  n_1000_h_hik[i] <- sum(x >= (c - h_hIK) & x <= (c + h_hIK))
  n_1000_h_cv[i] <- sum(x >= (c - h_CV) & x <= (c + h_CV))
  
  # Estimation of ATE using RDD
  rdd_result_IK <- rdrobust(y, x, c = 0, h = h_IK)
  ate_1000_ik[i] <- rdd_result_IK$coef[1]
  
  rdd_result_2IK <- rdrobust(y, x, c = 0, h = h_2IK)
  ate_1000_2ik[i] <- rdd_result_2IK$coef[1]
  
  rdd_result_hIK <- rdrobust(y, x, c = 0, h = h_hIK)
  ate_1000_hik[i] <- rdd_result_hIK$coef[1]
  
  rdd_result_CV <- rdrobust(y, x, c = 0, h = h_CV)
  ate_1000_cv[i] <- rdd_result_CV$coef[1]
}

# Compare the estimated treatment effects to the true treatment effect
mean_ate_1000_ik <- mean(ate_1000_ik)
mean_ate_1000_2ik <- mean(ate_1000_2ik)
mean_ate_1000_hik <- mean(ate_1000_hik)
mean_ate_1000_cv <- mean(ate_1000_cv)

bias_ate_1000_ik <- mean_ate_1000_ik - true_ate
bias_ate_1000_2ik <- mean_ate_1000_2ik - true_ate
bias_ate_1000_hik <- mean_ate_1000_hik - true_ate
bias_ate_1000_cv <- mean_ate_1000_cv - true_ate

var_ate_1000_ik <- var(ate_1000_ik)
var_ate_1000_2ik <- var(ate_1000_2ik)
var_ate_1000_hik <- var(ate_1000_hik)
var_ate_1000_cv <- var(ate_1000_cv)

mse_ate_1000_ik <- bias_ate_1000_ik^2 + var_ate_1000_ik
mse_ate_1000_2ik <- bias_ate_1000_2ik^2 + var_ate_1000_2ik
mse_ate_1000_hik <- bias_ate_1000_hik^2 + var_ate_1000_hik
mse_ate_1000_cv <- bias_ate_1000_cv^2 + var_ate_1000_cv

mean_h_1000_ik <- mean(h_1000_ik)
mean_h_1000_2ik <- mean(h_1000_2ik)
mean_h_1000_hik <- mean(h_1000_hik)
mean_h_1000_cv <- mean(h_1000_cv)

mean_n_h_1000_ik <- mean(n_1000_h_ik)
mean_n_h_1000_2ik <- mean(n_1000_h_2ik)
mean_n_h_1000_hik <- mean(n_1000_h_hik)
mean_n_h_1000_cv <- mean(n_1000_h_cv)

# Visualization of MSE in dependence of bandwidth 
data1000 = data.frame(sample_size = c(1000, 1000, 1000, 1000),
                      sample_within_h = c(mean_n_h_1000_ik,mean_n_h_1000_2ik, mean_n_h_1000_hik, mean_n_h_1000_cv),
                      procedure = c("IK", "2IK", "hIK", "CV" ),
                      bandwidth = c(mean_h_1000_ik, mean_h_1000_2ik, mean_h_1000_hik, mean_h_1000_cv),
                      bias = c(bias_ate_1000_ik, bias_ate_1000_2ik, bias_ate_1000_hik, bias_ate_1000_cv),
                      var = c(var_ate_1000_ik, var_ate_1000_2ik,var_ate_1000_hik, var_ate_1000_cv),
                      mse = c(mse_ate_1000_ik, mse_ate_1000_2ik, mse_ate_1000_hik, mse_ate_1000_cv))

# plot (Figure 3 (d))
ggplot(data1000, aes(x = bandwidth)) +
  geom_line(aes(y = mse, color = "MSE"), size = 1) +
  geom_line(aes(y = bias, color = "Bias"), size = 1) +
  geom_line(aes(y = var, color = "Variance"), size = 1) +
  geom_text(aes(label = procedure, y=0), vjust = -0.5, hjust = 0.5, size = 5) +  # Add labels for 
  geom_vline(xintercept = data1000$bandwidth, linetype = "dotted", color= "grey") +  # Add vertical lines at bandwidth points
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey") +
  labs(title = "", colour = "",
       x = "bandwidth", y="values"
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          legend.position = "none",  # Turn off legend
                          axis.line = element_line(color = "black"))  # Make axis lines black

# Save the plot
ggsave("bvm_h_1000.pdf")


# Save data1000 in latex table (Table 2(d))
latex_1000 <- stargazer(data1000, summary = FALSE)
writeLines(latex_1000, file.path("latex_1000.txt"))


################################################################################
# (6) Plot MSE, bias and variance for increasing sample size #


data_ik = data.frame(sample_size = c(100, 200, 500, 1000),
                     bias_ik = c(bias_ate_100_ik, bias_ate_200_ik, bias_ate_500_ik, bias_ate_1000_ik),
                     var_ik = c(var_ate_100_ik, var_ate_200_ik, var_ate_500_ik, var_ate_1000_ik),
                     mse_ik = c(mse_ate_100_ik, mse_ate_200_ik, mse_ate_500_ik, mse_ate_1000_ik)
)


ggplot(data_ik, aes(x = sample_size)) +
  geom_line(aes(y = mse_ik, color = "MSE"), size = 0.5) +
  geom_line(aes(y = bias_ik, color = "Bias"), size =0.5) +
  geom_line(aes(y = var_ik, color = "Var"), size =0.5 ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "",
       x = "sample size", y= "values", color = "", size = 10
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          axis.line = element_line(color = "black"))  # Make axis lines black

ggsave("bvm_ik.pdf")

# Save in Latex Table
latex_ik <- stargazer(data_ik, summary = FALSE)
writeLines(latex_ik, file.path("latex_ik.txt"))

####
data_2ik = data.frame(sample_size = c(100, 200, 500, 1000),
                      bias_2ik = c(bias_ate_100_2ik, bias_ate_200_2ik, bias_ate_500_2ik, bias_ate_1000_2ik),
                      var_2ik = c(var_ate_100_2ik, var_ate_200_2ik, var_ate_500_2ik, var_ate_1000_2ik),
                      mse_2ik = c(mse_ate_100_2ik, mse_ate_200_2ik, mse_ate_500_2ik, mse_ate_1000_2ik)
)
ggplot(data_2ik, aes(x = sample_size)) +
  geom_line(aes(y = mse_2ik, color = "MSE")) +
  geom_line(aes(y = bias_2ik, color = "Bias")) +
  geom_line(aes(y = var_2ik, color = "Var")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "",
       x = "sample size", y="values", color= ""
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          axis.line = element_line(color = "black"))  # Make axis lines black
ggsave("bvm_2ik.pdf")

# Save in Latex Table
latex_2ik <- stargazer(data_2ik, summary = FALSE)
writeLines(latex_2ik, file.path("latex_2ik.txt"))

####
data_hik = data.frame(sample_size = c(200, 500, 1000),
                      bias_hik = c(bias_ate_200_hik, bias_ate_500_hik, bias_ate_1000_hik),
                      var_hik = c(var_ate_200_hik, var_ate_500_hik, var_ate_1000_hik),
                      mse_hik = c(mse_ate_200_hik, mse_ate_500_hik, mse_ate_1000_hik)
)
ggplot(data_hik, aes(x = sample_size)) +
  geom_line(aes(y = mse_hik, color = "MSE")) +
  geom_line(aes(y = bias_hik, color = "Bias")) +
  geom_line(aes(y = var_hik, color = "Var")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       x = "sample size", y="values"
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          axis.line = element_line(color = "black"))  # Make axis lines black
ggsave("bvm_hik.pdf")

# Save in Latex Table
latex_hik <- stargazer(data_hik, summary = FALSE)
writeLines(latex_hik, file.path("latex_hik.txt"))

###
data_cv = data.frame(sample_size = c(100, 200, 500, 1000),
                     bias_cv = c(bias_ate_100_cv, bias_ate_200_cv, bias_ate_500_cv, bias_ate_1000_cv),
                     var_cv = c(var_ate_100_cv, var_ate_200_cv, var_ate_500_cv, var_ate_1000_cv),
                     mse_cv = c(mse_ate_100_cv, mse_ate_200_cv, mse_ate_500_cv, mse_ate_1000_cv)
)
ggplot(data_cv, aes(x = sample_size)) +
  geom_line(aes(y = mse_cv, color = "MSE")) +
  geom_line(aes(y = bias_cv, color = "Bias")) +
  geom_line(aes(y = var_cv, color = "Var")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       x = "sample size", y="values"
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          axis.line = element_line(color = "black"))  # Make axis lines black
ggsave("bvm_cv.pdf")

# Save in Latex Table
latex_cv <- stargazer(data_cv, summary = FALSE)
writeLines(latex_cv, file.path("latex_cv.txt"))
####




