setwd("~/Dropbox/Time_delay/SciAdv_revision/code-for-revision/Code_resubmit_S1") # set the working directory contains the code file named "Queueing_functions.R"
source('Queueing_functions.R')

set.seed(123) # set the random seed

# load the input data for estimation
input_matrix_pre = as.matrix(read.table("input_data_example.csv", sep = ",", header = TRUE))

num_trj = dim(input_matrix_pre)[2] - 1 # the number of trajectories used for the estimation
input_matrix = pmax(input_matrix_pre[,2:(num_trj+1)] - input_matrix_pre[1,2:(num_trj+1)], 0)
timespan_pre = input_matrix_pre[,1] - input_matrix_pre[1,1]

basal_epsilon = 0.001 # epsilon value to avoid zero variance of the likelihood function


prior_mean = c(1,1,1,1);
prior_var = c(10^6, 10^6, 10^6, 10^6);

theta_init = c(20, 0.1, 3, 1);
theta_hyper_init = c(1, 1, 1, 1, 1, 1, 1, 1)

## Iteration numbers ====================
burn = 0 # the length of the burn-in period, i.e., the number of first samples which are discarded to reduce the initial condition dependency of MCMC method.
jump = 1 # the reciprocal of the thinning rate. If jump = 10, then every 10 samples from after the burn-in period are chosen as posterior samples. 
effnum = 110 # the effective number of iterations of the MCMC method.

## Proposal variances ===================
prop_var_death = 0.0001
prop_var_n = 0.3
prop_var_tr = 0.3


## Input data setting ====================
inputdata1 <- rowMeans(input_matrix)[-1]
var_list1 = inputdata1/num_trj
var_list1 = var_list1 + basal_epsilon;
timespan = timespan_pre[-1]


## Perform estimation and save the output data =====================
mcmc.result1 <- MCMC_function(data = inputdata1, timespan = timespan, var_list = var_list1, prior_mean = prior_mean, prior_var = prior_var,
                              theta_init = theta_init, burn = burn, jump = jump, effnum = effnum, prop_var = c(prop_var_death, prop_var_n, prop_var_tr))

write.table(x = mcmc.result1$results, file = "post_samples_MBI_single_cell.csv", sep = ",", col.names = c("lambda_b", "lambda_d", "n", "t_r"), row.names = FALSE)
sd_list = c(sd(mcmc.result1$results[,1]),sd(mcmc.result1$results[,2]),sd(mcmc.result1$results[,3]),sd(mcmc.result1$results[,4]))
write.table(x = rbind(colMeans(mcmc.result1$results), sd_list), file = "mean_std_post_samples_MBI_single_cell.csv", sep = ",", col.names = c("lambda_b", "lambda_d", "n", "t_r"), row.names = FALSE)
# For detail description of the format of output files. Please see the README file in our package.

AR_death = colMeans(mcmc.result1$acceptance)[2]
AR_n = colMeans(mcmc.result1$acceptance)[3]
AR_tr = colMeans(mcmc.result1$acceptance)[4]



