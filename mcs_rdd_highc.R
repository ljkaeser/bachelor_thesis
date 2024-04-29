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

# define cutoff
c <- 0 

# define functions 
c_m_pos <- function(x)  - 3*(x+0.5)^2 + 8*(x+0.5)^3  + 4*(x+0.5)^5 + 0.625
c_m_neg <- function(x)  - 3*(x-0.5)^2 + 8*(x-0.5)^3  + 4*(x-0.5)^5 + 2.375


# Generate x values for plotting
c_x_pos <- seq(0, 0.5, length.out = 100)  # Generate x values for m_pos (x >= 0)
c_x_neg <- seq(-0.5, 0, length.out = 100) # Generate x values for m_neg (x < 0)

# calculate corresponding y values
c_y_pos <- c_m_pos(c_x_pos)
c_y_neg <- c_m_neg(c_x_neg)

# create a data frame for ggplot
c_data_pos <- data.frame(x = c_x_pos, y = c_y_pos)
c_data_neg <- data.frame(x = c_x_neg, y = c_y_neg)

# plot data generating functions (Figure 2(b))
ggplot() +
  geom_line(data = c_data_pos, aes(x = x, y = y), color = "blue", size =1) +
  geom_line(data = c_data_neg, aes(x = x, y = y), color = "red", size = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey", size = 1) +
  labs(title = "", color = "",
       x = "x",
       y = "y", size = 8) +
  scale_y_continuous(limits = c(-5, 5)) +  # Set the y-axis limits
  theme_minimal() + theme(#panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.title.x = element_text(size = 10),
    #axis.text.y = element_text(size = 10),
    axis.line = element_line(color = "black"))  # Make axis lines black

# Save the plot
ggsave("c_dgf.pdf")

# The real treatment effect we try to estimate in the following 
c_true_ate = c_m_pos(0)-c_m_neg(0)
print(c_true_ate)


###############################################################################

# (2) Run 1000 simulations for sample size 100 #

# define sample size and number of repeats 
n <- 100
n_reps <- 1000

# vector to store values in runs
c_ate_100_ik <- numeric(n_reps)
c_ate_100_cv <- numeric(n_reps)
c_ate_100_2ik <- numeric(n_reps)
c_ate_100_hik <- numeric(n_reps)

c_h_100_ik <- numeric(n_reps)
c_h_100_cv <- numeric(n_reps)
c_h_100_2ik <- numeric(n_reps)
c_h_100_hik <- numeric(n_reps)

c_n_100_h_ik <- numeric(n_reps)
c_n_100_h_2ik <- numeric(n_reps)
c_n_100_h_hik <- numeric(n_reps)
c_n_100_h_cv <- numeric(n_reps)

# loop
for (i in 1:n_reps) {
  set.seed(i)  # Set seed for reproducibility
  
  # Generate new data for each sample
  x <- rnorm(n) # generate observations from standard normal distribution (mean 0, sd 1)
  y <- ifelse(x < c, c_m_neg(x), c_m_pos(x) ) + rnorm(n, sd = 0.5) # add normal error term
  w <- ifelse(x<c,0,1) 
  
  # compute and store IK bandwidth 
  c_h_IK_rdrobust_2014 <- rdbwselect_2014(y,x,c=0,p=1, kernel = "tri", bwselect = "IK")
  c_h_IK <- c_h_IK_rdrobust_2014$bws[1]
  c_h_100_ik[i] <- c_h_IK 
  
  # compute and store 2IK bandwidth 
  c_h_2IK <- 2*c_h_IK
  c_h_100_2ik[i] <- c_h_2IK 
  
  # compute and store hIK bandwidth 
  c_h_hIK <- 0.5*c_h_IK
  c_h_100_hik[i] <- c_h_hIK
  
  # compute and store CV bandwidth 
  c_h_CV_rdrobust_2014 <- rdbwselect_2014(y=y, x = x, bwselect = "CV")
  c_h_CV <- c_h_CV_rdrobust_2014$bws[1]
  c_h_100_cv[i] <- c_h_CV
  
  # compute and store effectively used sample sizes
  c_n_100_h_ik[i] <- sum(x >= (c - c_h_IK) & x <= (c + c_h_IK))
  c_n_100_h_2ik[i] <- sum(x >= (c - c_h_2IK) & x <= (c + c_h_2IK))
  c_n_100_h_hik[i] <- sum(x >= (c - c_h_hIK) & x <= (c + c_h_hIK))
  c_n_100_h_cv[i] <- sum(x >= (c - c_h_CV) & x <= (c + c_h_CV))
  
  # estimate and store ATE with IK bandwidth 
  c_rdd_result_IK <- rdrobust(y, x, c = 0, h = c_h_IK)
  c_ate_100_ik[i] <- c_rdd_result_IK$coef[1]
  
  # estimate and store ATE with 2IK bandwidth 
  c_rdd_result_2IK <- rdrobust(y, x, c = 0, h = c_h_2IK)
  c_ate_100_2ik[i] <- c_rdd_result_2IK$coef[1]
  
  # estimate and store ATE with hIK bandwidth 
  #c_rdd_result_hIK <- rdrobust(y, x, c = 0, h = c_h_hIK)
  #c_ate_100_hik[i] <- c_rdd_result_hIK$coef[1]
  # !! too small
  
  # estimate and store ATE with CV bandwidth 
  c_rdd_result_CV <- rdrobust(y, x, c = 0, h = c_h_CV)
  c_ate_100_cv[i] <- c_rdd_result_CV$coef[1]
}

# mean of ate over runs
c_mean_ate_100_ik <- mean(c_ate_100_ik)
c_mean_ate_100_2ik <- mean(c_ate_100_2ik)
c_mean_ate_100_hik <- mean(c_ate_100_hik)
c_mean_ate_100_cv <- mean(c_ate_100_cv)

# bias over runs
c_bias_ate_100_ik <- c_mean_ate_100_ik - c_true_ate
c_bias_ate_100_2ik <- c_mean_ate_100_2ik - c_true_ate
c_bias_ate_100_hik <- c_mean_ate_100_hik - c_true_ate
c_bias_ate_100_cv <- c_mean_ate_100_cv - c_true_ate

# var over runs
c_var_ate_100_ik <- var(c_ate_100_ik)
c_var_ate_100_2ik <- var(c_ate_100_2ik)
c_var_ate_100_hik <- var(c_ate_100_hik)
c_var_ate_100_cv <- var(c_ate_100_cv)

# mse over runs
c_mse_ate_100_ik <- c_bias_ate_100_ik^2 + c_var_ate_100_ik
c_mse_ate_100_2ik <- c_bias_ate_100_2ik^2 + c_var_ate_100_2ik
c_mse_ate_100_hik <- c_bias_ate_100_hik^2 + c_var_ate_100_hik
c_mse_ate_100_cv <- c_bias_ate_100_cv^2 + c_var_ate_100_cv

# mean of h over runs
c_mean_h_100_ik <- mean(c_h_100_ik)
c_mean_h_100_2ik <- mean(c_h_100_2ik)
c_mean_h_100_hik <- mean(c_h_100_hik)
c_mean_h_100_cv <- mean(c_h_100_cv)

# mean of effective sample size over runs
c_mean_n_h_100_ik <- mean(c_n_100_h_ik)
c_mean_n_h_100_2ik <- mean(c_n_100_h_2ik)
c_mean_n_h_100_hik <- mean(c_n_100_h_hik)
c_mean_n_h_100_cv <- mean(c_n_100_h_cv)

# visualization of MSE in dependence of bandwidth 
c_data100 = data.frame(sample_size = c(100, 100, 100),
                       sample_within_h = c(c_mean_n_h_100_ik,c_mean_n_h_100_2ik, c_mean_n_h_100_cv),
                       procedure = c("IK", "2IK", "CV" ),
                       bandwidth = c(c_mean_h_100_ik, c_mean_h_100_2ik, c_mean_h_100_cv),
                       bias = c(c_bias_ate_100_ik, c_bias_ate_100_2ik, c_bias_ate_100_cv),
                       var = c(c_var_ate_100_ik, c_var_ate_100_2ik, c_var_ate_100_cv),
                       mse = c(c_mse_ate_100_ik, c_mse_ate_100_2ik, c_mse_ate_100_cv))

# plot Figure 4(a)
ggplot(c_data100, aes(x = bandwidth)) +
  geom_line(aes(y = mse, color = "MSE"), size = 1) +
  geom_line(aes(y = bias, color = "Bias"), size = 1) +
  geom_line(aes(y = var, color = "Variance"), size = 1) +
  #geom_text(aes(label = procedure, y=0), vjust = -0.5, hjust = 0.5, size = 5) +  # Add labels for 
  #geom_text(data = c_data100[c_data100$procedure == "IK", ], aes(label = procedure), x = c_data100$bandwidth[1] + 10, y = 10, size = 5) + # Adjust position for "IK" label
  geom_text(aes(x = 0.47, y = 2.4, label = "IK"), size = 5) +
  geom_text(aes(x = 0.44, y = 2.4, label = "CV"), size = 5) +
  geom_text(aes(x = 0.922, y = 2.4, label = "2IK"), size = 5) +
  geom_vline(xintercept = c_data100$bandwidth, linetype = "dotted", color= "grey") +  # Add vertical lines at bandwidth points
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       legend.position = "none",  # Turn off legend
       x = "bandwidth", y="values", size = 8
  ) +
  theme_minimal() + theme(#legend.key.size = unit(2, "lines"),  # Increase size of legend key
    #legend.text = element_text(size = 12),  # Increase size of legend text
    legend.position = "none",  # Turn off legend
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"))  # Make axis lines black


# Save the plot
ggsave("c_bvm_h_100.pdf")


# Save data100 in Latex Table (for Table 3(a))
c_latex_100 <- stargazer(c_data100, summary = FALSE)
writeLines(c_latex_100, file.path("c_latex_100.txt"))

################################################################################

# (3) Run 1000 simulations for sample size 200 #

# define sample size and number of repeats 
n <- 200
n_reps <- 1000

# vector to store values 
c_ate_200_ik <- numeric(n_reps)
c_ate_200_cv <- numeric(n_reps)
c_ate_200_2ik <- numeric(n_reps)
c_ate_200_hik <- numeric(n_reps)

c_h_200_ik <- numeric(n_reps)
c_h_200_cv <- numeric(n_reps)
c_h_200_2ik <- numeric(n_reps)
c_h_200_hik <- numeric(n_reps)

c_n_200_h_ik <- numeric(n_reps)
c_n_200_h_2ik <- numeric(n_reps)
c_n_200_h_hik <- numeric(n_reps)
c_n_200_h_cv <- numeric(n_reps)

# loop
for (i in 1:n_reps) {
  set.seed(i)  # Set seed for reproducibility
  
  # Generate new data for each sample
  x <- rnorm(n) # generate observations from standard normal distribution (mean 0, sd 1)
  y <- ifelse(x < c, c_m_neg(x), c_m_pos(x) ) + rnorm(n, sd = 0.5) # add normal error term
  w <- ifelse(x<c,0,1) 
  
  # IK bandwidth 
  c_h_IK_rdrobust_2014 <- rdbwselect_2014(y,x,c=0,p=1, kernel = "tri", bwselect = "IK")
  c_h_IK <- c_h_IK_rdrobust_2014$bws[1]
  c_h_200_ik[i] <- c_h_IK 
  
  # 2IK bandwidth 
  c_h_2IK <- 2*c_h_IK
  c_h_200_2ik[i] <- c_h_2IK 
  
  # hIK bandwidth 
  c_h_hIK <- 0.5*c_h_IK
  c_h_200_hik[i] <- c_h_hIK
  
  # CV bandwidth 
  c_h_CV_rdrobust_2014 <- rdbwselect_2014(y=y, x = x, bwselect = "CV")
  c_h_CV <- c_h_CV_rdrobust_2014$bws[1]
  c_h_200_cv[i] <- c_h_CV
  
  # effective sample sizes
  c_n_200_h_ik[i] <- sum(x >= (c - c_h_IK) & x <= (c + c_h_IK))
  c_n_200_h_2ik[i] <- sum(x >= (c - c_h_2IK) & x <= (c + c_h_2IK))
  c_n_200_h_hik[i] <- sum(x >= (c - c_h_hIK) & x <= (c + c_h_hIK))
  c_n_200_h_cv[i] <- sum(x >= (c - c_h_CV) & x <= (c + c_h_CV))
  
  # estimation of ATE with h_IK
  c_rdd_result_IK <- rdrobust(y, x, c = 0, h = c_h_IK)
  c_ate_200_ik[i] <- c_rdd_result_IK$coef[1]
  
  # estimation of ATE with h_2IK
  c_rdd_result_2IK <- rdrobust(y, x, c = 0, h = c_h_2IK)
  c_ate_200_2ik[i] <- c_rdd_result_2IK$coef[1]
  
  # estimation of ATE with h_hIK
  c_rdd_result_hIK <- rdrobust(y, x, c = 0, h = c_h_hIK)
  c_ate_200_hik[i] <- c_rdd_result_hIK$coef[1]
  
  # estimation of ATE with h_CV
  c_rdd_result_CV <- rdrobust(y, x, c = 0, h = c_h_CV)
  c_ate_200_cv[i] <- c_rdd_result_CV$coef[1]
}

# mean of ate 
c_mean_ate_200_ik <- mean(c_ate_200_ik)
c_mean_ate_200_2ik <- mean(c_ate_200_2ik)
c_mean_ate_200_hik <- mean(c_ate_200_hik)
c_mean_ate_200_cv <- mean(c_ate_200_cv)

# bias 
c_bias_ate_200_ik <- c_mean_ate_200_ik - c_true_ate
c_bias_ate_200_2ik <- c_mean_ate_200_2ik - c_true_ate
c_bias_ate_200_hik <- c_mean_ate_200_hik - c_true_ate
c_bias_ate_200_cv <- c_mean_ate_200_cv - c_true_ate

# variance
c_var_ate_200_ik <- var(c_ate_200_ik)
c_var_ate_200_2ik <- var(c_ate_200_2ik)
c_var_ate_200_hik <- var(c_ate_200_hik)
c_var_ate_200_cv <- var(c_ate_200_cv)

# mse
c_mse_ate_200_ik <- c_bias_ate_200_ik^2 + c_var_ate_200_ik
c_mse_ate_200_2ik <- c_bias_ate_200_2ik^2 + c_var_ate_200_2ik
c_mse_ate_200_hik <- c_bias_ate_200_hik^2 + c_var_ate_200_hik
c_mse_ate_200_cv <- c_bias_ate_200_cv^2 + c_var_ate_200_cv

# mean h 
c_mean_h_200_ik <- mean(c_h_200_ik)
c_mean_h_200_2ik <- mean(c_h_200_2ik)
c_mean_h_200_hik <- mean(c_h_200_hik)
c_mean_h_200_cv <- mean(c_h_200_cv)

# mean effective sample size 
c_mean_n_h_200_ik <- mean(c_n_200_h_ik)
c_mean_n_h_200_2ik <- mean(c_n_200_h_2ik)
c_mean_n_h_200_hik <- mean(c_n_200_h_hik)
c_mean_n_h_200_cv <- mean(c_n_200_h_cv)

# visualization of MSE in dependence of bandwidth 
c_data200 = data.frame(sample_size = c(200, 200, 200, 200),
                       sample_within_h = c(c_mean_n_h_200_ik,c_mean_n_h_200_2ik, c_mean_n_h_200_hik, c_mean_n_h_200_cv),
                       procedure = c("IK", "2IK", "hIK", "CV" ),
                       bandwidth = c(c_mean_h_200_ik, c_mean_h_200_2ik, c_mean_h_200_hik, c_mean_h_200_cv),
                       bias = c(c_bias_ate_200_ik, c_bias_ate_200_2ik, c_bias_ate_200_hik, c_bias_ate_200_cv),
                       var = c(c_var_ate_200_ik, c_var_ate_200_2ik, c_var_ate_200_hik, c_var_ate_200_cv),
                       mse = c(c_mse_ate_200_ik, c_mse_ate_200_2ik, c_mse_ate_200_hik, c_mse_ate_200_cv))

# plot Figure 4(b)
ggplot(c_data200, aes(x = bandwidth)) +
  geom_line(aes(y = mse, color = "MSE"), size = 1) +
  geom_line(aes(y = bias, color = "Bias"), size = 1) +
  geom_line(aes(y = var, color = "Variance"), size = 1) +
  geom_text(aes(label = procedure, y=0), vjust = -0.5, hjust = 0.5, size = 5) +  # Add labels for 
  geom_vline(xintercept = c_data200$bandwidth, linetype = "dotted", color= "grey") +  # Add vertical lines at bandwidth points
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       x = "bandwidth", y="values"
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          legend.position = "none",  # Turn off legend
                          axis.line = element_line(color = "black"))  # Make axis lines black

# Save the plot
ggsave("c_bvm_h_200.pdf")


# Save data100 in Latex Table (Table 3(b))
c_latex_200 <- stargazer(c_data200, summary = FALSE)
writeLines(c_latex_200, file.path("c_latex_200.txt"))

################################################################################

# (4) Run 1000 simulations for sample size 500 #

# define sample size and number of repeats 
n <- 500
n_reps <- 1000

# vector to store values
c_ate_500_ik <- numeric(n_reps)
c_ate_500_cv <- numeric(n_reps)
c_ate_500_2ik <- numeric(n_reps)
c_ate_500_hik <- numeric(n_reps)

c_h_500_ik <- numeric(n_reps)
c_h_500_cv <- numeric(n_reps)
c_h_500_2ik <- numeric(n_reps)
c_h_500_hik <- numeric(n_reps)

c_n_500_h_ik <- numeric(n_reps)
c_n_500_h_2ik <- numeric(n_reps)
c_n_500_h_hik <- numeric(n_reps)
c_n_500_h_cv <- numeric(n_reps)

# loop
for (i in 1:n_reps) {
  set.seed(i)  # Set seed for reproducibility
  
  # Generate new data for each sample
  x <- rnorm(n) # generate observations from standard normal distribution (mean 0, sd 1)
  y <- ifelse(x < c, c_m_neg(x), c_m_pos(x) ) + rnorm(n, sd = 0.5) # add normal error term
  w <- ifelse(x<c,0,1) 
  
  # IK bandwidth
  c_h_IK_rdrobust_2014 <- rdbwselect_2014(y,x,c=0,p=1, kernel = "tri", bwselect = "IK")
  c_h_IK <- c_h_IK_rdrobust_2014$bws[1]
  c_h_500_ik[i] <- c_h_IK 
  
  # 2IK bandwidth
  c_h_2IK <- 2*c_h_IK
  c_h_500_2ik[i] <- c_h_2IK 
  
  # hIK bandwidth
  c_h_hIK <- 0.5*c_h_IK
  c_h_500_hik[i] <- c_h_hIK
  
  # CV bandwidth
  c_h_CV_rdrobust_2014 <- rdbwselect_2014(y=y, x = x, bwselect = "CV")
  c_h_CV <- c_h_CV_rdrobust_2014$bws[1]
  c_h_500_cv[i] <- c_h_CV
  
  # effective sample size
  c_n_500_h_ik[i] <- sum(x >= (c - c_h_IK) & x <= (c + c_h_IK))
  c_n_500_h_2ik[i] <- sum(x >= (c - c_h_2IK) & x <= (c + c_h_2IK))
  c_n_500_h_hik[i] <- sum(x >= (c - c_h_hIK) & x <= (c + c_h_hIK))
  c_n_500_h_cv[i] <- sum(x >= (c - c_h_CV) & x <= (c + c_h_CV))
  
  # estimation of ATE with h_IK
  c_rdd_result_IK <- rdrobust(y, x, c = 0, h = c_h_IK)
  c_ate_500_ik[i] <- c_rdd_result_IK$coef[1]
  
  # estimation of ATE with h_2IK
  c_rdd_result_2IK <- rdrobust(y, x, c = 0, h = c_h_2IK)
  c_ate_500_2ik[i] <- c_rdd_result_2IK$coef[1]
  
  # estimation of ATE with h_hIK
  c_rdd_result_hIK <- rdrobust(y, x, c = 0, h = c_h_hIK)
  c_ate_500_hik[i] <- c_rdd_result_hIK$coef[1]
  
  # estimation of ATE with h_CV
  c_rdd_result_CV <- rdrobust(y, x, c = 0, h = c_h_CV)
  c_ate_500_cv[i] <- c_rdd_result_CV$coef[1]
}

# mean of ate
c_mean_ate_500_ik <- mean(c_ate_500_ik)
c_mean_ate_500_2ik <- mean(c_ate_500_2ik)
c_mean_ate_500_hik <- mean(c_ate_500_hik)
c_mean_ate_500_cv <- mean(c_ate_500_cv)

# bias
c_bias_ate_500_ik <- c_mean_ate_500_ik - c_true_ate
c_bias_ate_500_2ik <- c_mean_ate_500_2ik - c_true_ate
c_bias_ate_500_hik <- c_mean_ate_500_hik - c_true_ate
c_bias_ate_500_cv <- c_mean_ate_500_cv - c_true_ate

# var
c_var_ate_500_ik <- var(c_ate_500_ik)
c_var_ate_500_2ik <- var(c_ate_500_2ik)
c_var_ate_500_hik <- var(c_ate_500_hik)
c_var_ate_500_cv <- var(c_ate_500_cv)

# mse
c_mse_ate_500_ik <- c_bias_ate_500_ik^2 + c_var_ate_500_ik
c_mse_ate_500_2ik <- c_bias_ate_500_2ik^2 + c_var_ate_500_2ik
c_mse_ate_500_hik <- c_bias_ate_500_hik^2 + c_var_ate_500_hik
c_mse_ate_500_cv <- c_bias_ate_500_cv^2 + c_var_ate_500_cv

# mean h 
c_mean_h_500_ik <- mean(c_h_500_ik)
c_mean_h_500_2ik <- mean(c_h_500_2ik)
c_mean_h_500_hik <- mean(c_h_500_hik)
c_mean_h_500_cv <- mean(c_h_500_cv)

# mean effective sample size 
c_mean_n_h_500_ik <- mean(c_n_500_h_ik)
c_mean_n_h_500_2ik <- mean(c_n_500_h_2ik)
c_mean_n_h_500_hik <- mean(c_n_500_h_hik)
c_mean_n_h_500_cv <- mean(c_n_500_h_cv)

# visualization of MSE in dependence of bandwidth 
c_data500 = data.frame(sample_size = c(500, 500, 500, 500),
                       sample_within_h = c(c_mean_n_h_500_ik,c_mean_n_h_500_2ik, c_mean_n_h_500_hik, c_mean_n_h_500_cv),
                       procedure = c("IK", "2IK", "hIK", "CV" ),
                       bandwidth = c(c_mean_h_500_ik, c_mean_h_500_2ik, c_mean_h_500_hik, c_mean_h_500_cv),
                       bias = c(c_bias_ate_500_ik, c_bias_ate_500_2ik, c_bias_ate_500_hik, c_bias_ate_500_cv),
                       var = c(c_var_ate_500_ik, c_var_ate_500_2ik,c_var_ate_500_hik, c_var_ate_500_cv),
                       mse = c(c_mse_ate_500_ik, c_mse_ate_500_2ik, c_mse_ate_500_hik, c_mse_ate_500_cv))

# plot Figure 4(c)
ggplot(c_data500, aes(x = bandwidth)) +
  geom_line(aes(y = mse, color = "MSE"), size =1 ) +
  geom_line(aes(y = bias, color = "Bias"), size = 1) +
  geom_line(aes(y = var, color = "Variance"), size =1) +
  #geom_text(aes(label = procedure, y=0), vjust = -0.5, hjust = 0.5, size = 5) +  # Add labels for 
  geom_text(aes(x = 0.295, y = 0.25, label = "IK", fontface = "plain"), size = 5) +
  geom_text(aes(x = 0.125, y = 0.25, label = "CV", fontface = "plain"), size = 5) +
  geom_text(aes(x = 0.590, y = 0.25, label = "2IK"), size = 5) +
  geom_text(aes(x = 0.15, y = 0.25, label = "hIK"), size = 5) +
  geom_vline(xintercept = c_data500$bandwidth, linetype = "dotted", color= "grey") +  # Add vertical lines at bandwidth points
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       x = "bandwidth", y="values", size = 8
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          legend.position = "none",  # Turn off legend
                          axis.line = element_line(color = "black"))  # Make axis lines black

# Save the plot
ggsave("c_bvm_h_500.pdf")


# Save data100 in Latex Table (Table 3(c))
c_latex_500 <- stargazer(c_data500, summary = FALSE)
writeLines(c_latex_500, file.path("c_latex_500.txt"))

#########################################################################

# (5) Run 1000 simulations for sample size 1000 #

# define sample size and number of repeats 
n <- 1000
n_reps <- 1000

# vector to store values
c_ate_1000_ik <- numeric(n_reps)
c_ate_1000_cv <- numeric(n_reps)
c_ate_1000_2ik <- numeric(n_reps)
c_ate_1000_hik <- numeric(n_reps)

c_h_1000_ik <- numeric(n_reps)
c_h_1000_cv <- numeric(n_reps)
c_h_1000_2ik <- numeric(n_reps)
c_h_1000_hik <- numeric(n_reps)

c_n_1000_h_ik <- numeric(n_reps)
c_n_1000_h_2ik <- numeric(n_reps)
c_n_1000_h_hik <- numeric(n_reps)
c_n_1000_h_cv <- numeric(n_reps)

# loop
for (i in 1:n_reps) {
  set.seed(i)  # Set seed for reproducibility
  
  # Generate new data for each sample
  x <- rnorm(n) # generate observations from standard normal distribution (mean 0, sd 1)
  y <- ifelse(x < c, c_m_neg(x), c_m_pos(x) ) + rnorm(n, sd = 0.5) # add normal error term
  w <- ifelse(x<c,0,1) 
  
  # IK bandwidth 
  c_h_IK_rdrobust_2014 <- rdbwselect_2014(y,x,c=0,p=1, kernel = "tri", bwselect = "IK")
  c_h_IK <- c_h_IK_rdrobust_2014$bws[1]
  c_h_1000_ik[i] <- c_h_IK 
  
  # 2IK bandwidth 
  c_h_2IK <- 2*c_h_IK
  c_h_1000_2ik[i] <- c_h_2IK 
  
  # hIK bandwidth 
  c_h_hIK <- 0.5*c_h_IK
  c_h_1000_hik[i] <- c_h_hIK
  
  # CV bandwidth 
  c_h_CV_rdrobust_2014 <- rdbwselect_2014(y=y, x = x, bwselect = "CV")
  c_h_CV <- c_h_CV_rdrobust_2014$bws[1]
  c_h_1000_cv[i] <- c_h_CV
  
  # effective sample sizes
  c_n_1000_h_ik[i] <- sum(x >= (c - c_h_IK) & x <= (c + c_h_IK))
  c_n_1000_h_2ik[i] <- sum(x >= (c - c_h_2IK) & x <= (c + c_h_2IK))
  c_n_1000_h_hik[i] <- sum(x >= (c - c_h_hIK) & x <= (c + c_h_hIK))
  c_n_1000_h_cv[i] <- sum(x >= (c - c_h_CV) & x <= (c + c_h_CV))
  
  # estimation of ATE using h_IK
  c_rdd_result_IK <- rdrobust(y, x, c = 0, h = c_h_IK)
  c_ate_1000_ik[i] <- c_rdd_result_IK$coef[1]
  
  # estimation of ATE using h_2IK
  c_rdd_result_2IK <- rdrobust(y, x, c = 0, h = c_h_2IK)
  c_ate_1000_2ik[i] <- c_rdd_result_2IK$coef[1]
  
  # estimation of ATE using h_hIK
  c_rdd_result_hIK <- rdrobust(y, x, c = 0, h = c_h_hIK)
  c_ate_1000_hik[i] <- c_rdd_result_hIK$coef[1]
  
  # estimation of ATE using h_CV
  c_rdd_result_CV <- rdrobust(y, x, c = 0, h = c_h_CV)
  c_ate_1000_cv[i] <- c_rdd_result_CV$coef[1]
}

# mean of ate 
c_mean_ate_1000_ik <- mean(c_ate_1000_ik)
c_mean_ate_1000_2ik <- mean(c_ate_1000_2ik)
c_mean_ate_1000_hik <- mean(c_ate_1000_hik)
c_mean_ate_1000_cv <- mean(c_ate_1000_cv)

# bias
c_bias_ate_1000_ik <- c_mean_ate_1000_ik - c_true_ate
c_bias_ate_1000_2ik <- c_mean_ate_1000_2ik - c_true_ate
c_bias_ate_1000_hik <- c_mean_ate_1000_hik - c_true_ate
c_bias_ate_1000_cv <- c_mean_ate_1000_cv - c_true_ate

# variance
c_var_ate_1000_ik <- var(c_ate_1000_ik)
c_var_ate_1000_2ik <- var(c_ate_1000_2ik)
c_var_ate_1000_hik <- var(c_ate_1000_hik)
c_var_ate_1000_cv <- var(c_ate_1000_cv)

# mse 
c_mse_ate_1000_ik <- c_bias_ate_1000_ik^2 + c_var_ate_1000_ik
c_mse_ate_1000_2ik <- c_bias_ate_1000_2ik^2 + c_var_ate_1000_2ik
c_mse_ate_1000_hik <- c_bias_ate_1000_hik^2 + c_var_ate_1000_hik
c_mse_ate_1000_cv <- c_bias_ate_1000_cv^2 + c_var_ate_1000_cv

# mean h 
c_mean_h_1000_ik <- mean(c_h_1000_ik)
c_mean_h_1000_2ik <- mean(c_h_1000_2ik)
c_mean_h_1000_hik <- mean(c_h_1000_hik)
c_mean_h_1000_cv <- mean(c_h_1000_cv)

# mean effective sample size
c_mean_n_h_1000_ik <- mean(c_n_1000_h_ik)
c_mean_n_h_1000_2ik <- mean(c_n_1000_h_2ik)
c_mean_n_h_1000_hik <- mean(c_n_1000_h_hik)
c_mean_n_h_1000_cv <- mean(c_n_1000_h_cv)

# visualization of MSE in dependence of bandwidth 
c_data1000 = data.frame(sample_size = c(1000, 1000, 1000, 1000),
                        sample_within_h = c(c_mean_n_h_1000_ik,c_mean_n_h_1000_2ik, c_mean_n_h_1000_hik, c_mean_n_h_1000_cv),
                        procedure = c("IK", "2IK", "hIK", "CV" ),
                        bandwidth = c(c_mean_h_1000_ik, c_mean_h_1000_2ik, c_mean_h_1000_hik, c_mean_h_1000_cv),
                        bias = c(c_bias_ate_1000_ik, c_bias_ate_1000_2ik, c_bias_ate_1000_hik, c_bias_ate_1000_cv),
                        var = c(c_var_ate_1000_ik, c_var_ate_1000_2ik,c_var_ate_1000_hik, c_var_ate_1000_cv),
                        mse = c(c_mse_ate_1000_ik, c_mse_ate_1000_2ik, c_mse_ate_1000_hik, c_mse_ate_1000_cv))

# plot Figure 4(d)
ggplot(c_data1000, aes(x = bandwidth)) +
  geom_line(aes(y = mse, color = "MSE"), size = 1) +
  geom_line(aes(y = bias, color = "Bias"), size =1) +
  geom_line(aes(y = var, color = "Variance"), size =1) +
  geom_text(aes(label = procedure, y=0), vjust = -0.5, hjust = 0.5, size = 5) +  # Add labels for 
  geom_vline(xintercept = c_data1000$bandwidth, linetype = "dotted", color= "grey") +  # Add vertical lines at bandwidth points
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       x = "bandwidth", y="values", size = 8
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          legend.position = "none",  # Turn off legend
                          axis.line = element_line(color = "black"))  # Make axis lines black

# Save the plot
ggsave("c_bvm_h_1000.pdf")


# Save data100 in Latex Table (Table 3(d))
c_latex_1000 <- stargazer(c_data1000, summary = FALSE)
writeLines(c_latex_1000, file.path("c_latex_1000.txt"))


################################################################################

# (6) Plot MSE, bias and variance for increasing sample size #

# create data set for IK bandwidth with all sample sizes
c_data_ik = data.frame(sample_size = c(100, 200, 500, 1000),
                       bias_ik = c(c_bias_ate_100_ik, c_bias_ate_200_ik, c_bias_ate_500_ik, c_bias_ate_1000_ik),
                       var_ik = c(c_var_ate_100_ik, c_var_ate_200_ik, c_var_ate_500_ik, c_var_ate_1000_ik),
                       mse_ik = c(c_mse_ate_100_ik, c_mse_ate_200_ik, c_mse_ate_500_ik, c_mse_ate_1000_ik)
)

# plot Figure 6
ggplot(c_data_ik, aes(x = sample_size)) +
  geom_line(aes(y = mse_ik, color = "MSE"), size = 0.5) +
  geom_line(aes(y = bias_ik, color = "Bias"), size =0.5) +
  geom_line(aes(y = var_ik, color = "Var"), size = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       x = "sample size", y= "values", color = "", size = 10
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          axis.line = element_line(color = "black"))  # Make axis lines black

ggsave("c_bvm_ik.pdf")

# Save in Latex Table (Table 5)
c_latex_ik <- stargazer(c_data_ik, summary = FALSE)
writeLines(c_latex_ik, file.path("c_latex_ik.txt"))

### Do the same for other bandwidths (not included in thesis) ###
# 2IK #
c_data_2ik = data.frame(sample_size = c(100, 200, 500, 1000),
                        bias_2ik = c(c_bias_ate_100_2ik, c_bias_ate_200_2ik, c_bias_ate_500_2ik, c_bias_ate_1000_2ik),
                        var_2ik = c(c_var_ate_100_2ik, c_var_ate_200_2ik, c_var_ate_500_2ik, c_var_ate_1000_2ik),
                        mse_2ik = c(c_mse_ate_100_2ik, c_mse_ate_200_2ik, c_mse_ate_500_2ik, c_mse_ate_1000_2ik)
)
ggplot(c_data_2ik, aes(x = sample_size)) +
  geom_line(aes(y = mse_2ik, color = "MSE"),size = 0.5) +
  geom_line(aes(y = bias_2ik, color = "Bias"),size = 0.5) +
  geom_line(aes(y = var_2ik, color = "Var"),size = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       x = "sample size", y="values", size = 10
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          axis.line = element_line(color = "black"))  # Make axis lines black
ggsave("c_bvm_2ik.pdf")

# Save in Latex Table
c_latex_2ik <- stargazer(c_data_2ik, summary = FALSE)
writeLines(c_latex_2ik, file.path("c_latex_2ik.txt"))

# hIK #
c_data_hik = data.frame(sample_size = c(200, 500, 1000),
                        bias_hik = c(c_bias_ate_200_hik, c_bias_ate_500_hik, c_bias_ate_1000_hik),
                        var_hik = c(c_var_ate_200_hik, c_var_ate_500_hik, c_var_ate_1000_hik),
                        mse_hik = c(c_mse_ate_200_hik, c_mse_ate_500_hik, c_mse_ate_1000_hik)
)
ggplot(c_data_hik, aes(x = sample_size)) +
  geom_line(aes(y = mse_hik, color = "MSE"), size = 0.5) +
  geom_line(aes(y = bias_hik, color = "Bias"), size = 0.5) +
  geom_line(aes(y = var_hik, color = "Var"), size = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       x = "sample size", y="values", size = 10
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          axis.line = element_line(color = "black"))  # Make axis lines black
ggsave("c_bvm_hik.pdf")

# Save in Latex Table
c_latex_hik <- stargazer(c_data_hik, summary = FALSE)
writeLines(c_latex_hik, file.path("c_latex_hik.txt"))

# CV # 
c_data_cv = data.frame(sample_size = c(100, 200, 500, 1000),
                       bias_cv = c(c_bias_ate_100_cv, c_bias_ate_200_cv, c_bias_ate_500_cv, c_bias_ate_1000_cv),
                       var_cv = c(c_var_ate_100_cv, c_var_ate_200_cv, c_var_ate_500_cv, c_var_ate_1000_cv),
                       mse_cv = c(c_mse_ate_100_cv, c_mse_ate_200_cv, c_mse_ate_500_cv, c_mse_ate_1000_cv)
)
ggplot(c_data_cv, aes(x = sample_size)) +
  geom_line(aes(y = mse_cv, color = "MSE"), size = 0.5) +
  geom_line(aes(y = bias_cv, color = "Bias"), size = 0.5) +
  geom_line(aes(y = var_cv, color = "Var"), size =0.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  labs(title = "", color = "",
       x = "sample size", y="values", size = 10
  ) +
  theme_minimal() + theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank(),  # Remove minor grid lines
                          axis.line = element_line(color = "black"))  # Make axis lines black
ggsave("c_bvm_cv.pdf")

# Save in Latex Table
c_latex_cv <- stargazer(c_data_cv, summary = FALSE)
writeLines(c_latex_cv, file.path("c_latex_cv.txt"))