setwd("~/Dropbox/Time_delay/SciAdv_revision/code-for-revision/Code_resubmit_S1") # set the working directory contains the code file named "Queueing_functions.R"
source('Queueing_functions.R')

set.seed(123) # set the random seed

# load the input data for estimation
input_matrix_pre = as.matrix(read.table("input_data_example.csv", sep = ",", header = TRUE))

num_trj = dim(input_matrix_pre)[2] - 1 # the number of trajectories used for the estimation
input_matrix = pmax(input_matrix_pre[,2:(num_trj+1)] - input_matrix_pre[1,2:(num_trj+1)], 0)
timespan_pre = input_matrix_pre[,1] - input_matrix_pre[1,1]

basal_epsilon = 0.001 # epsilon value to avoid zero variance of the likelihood function

hyper_prior_mean = c(1,1,1,1,1,1,1,1)
hyper_prior_var = c(10^6,10^6,10^6,10^6,10^6,10^6,10^6,10^6)

theta_init = c(20, 0.1, 3, 1);
theta_hyper_init = c(1, 1, 1, 1, 1, 1, 1, 1)


## Iteration numbers ====================
burn = 0 # the length of the burn-in period, i.e., the number of first samples which are discarded to reduce the initial condition dependency of MCMC method.
jump = 1 # the reciprocal of the thinning rate. If jump = 10, then every 10 samples from after the burn-in period are chosen as posterior samples. 
effnum = 110 # the effective number of iterations of the MCMC method.


## Proposal variances ===================
prop_var_birth = 2
prop_var_death = 0.000001
prop_var_n = 0.05
prop_var_tr = 0.02
hyper_prop_var = c(3, 0.005, 3, 0.00001, 1, 0.1, 0.5, 0.003)



## Input data setting ====================
inputdata1 <- as.matrix(input_matrix[-1, 1:num_trj])
var_list1 = inputdata1
var_list1[1:(nrow(inputdata1)-2),] = var_list1[1:(nrow(inputdata1)-2),] + inputdata1[3:nrow(inputdata1)]
var_list1[1:(nrow(inputdata1)-1),] = var_list1[1:(nrow(inputdata1)-1),] + inputdata1[2:nrow(inputdata1)]
var_list1[2:nrow(inputdata1),] = var_list1[2:nrow(inputdata1),] + inputdata1[1:(nrow(inputdata1)-1)]
var_list1[3:nrow(inputdata1),] = var_list1[3:nrow(inputdata1),] + inputdata1[1:(nrow(inputdata1)-2)]
var_list1[1,] = var_list1[1,]/3
var_list1[2,] = var_list1[2,]/4
var_list1[3:(nrow(inputdata1)-2),] = var_list1[3:(nrow(inputdata1)-2),]/5
var_list1[nrow(inputdata1)-1,] = var_list1[nrow(inputdata1)-1,]/4
var_list1[nrow(inputdata1),] = var_list1[nrow(inputdata1),]/3
var_list1 = var_list1 + basal_epsilon;

timespan = timespan_pre[-1]


## Perform estimation and save the output data =====================
mcmc.result1 <- MCMC_function_ME(data = inputdata1, timespan = timespan, var_list = var_list1, hyper_prior_mean = hyper_prior_mean, hyper_prior_var = hyper_prior_var, 
                                 theta_init = theta_init, theta_hyper_init = theta_hyper_init, burn = burn, jump = jump, effnum = effnum, 
                                 prop_var = c(prop_var_birth, prop_var_death, prop_var_n, prop_var_tr), hyper_prop_var = hyper_prop_var)

for(jj in 1:num_trj){
  write.table(x = cbind(mcmc.result1$A_sample[,jj],mcmc.result1$B_sample[,jj],mcmc.result1$alpha_sample[,jj],mcmc.result1$beta_sample[,jj]), 
              file = paste("post_samples_",jj,"th_cell_MBI_multiple_cells.csv", sep = ""), sep = ",", col.names = c("lambda_b", "lambda_d", "n", "t_r"), row.names = FALSE)
  sd_list = c(sd(mcmc.result1$A_sample[,jj]),sd(mcmc.result1$B_sample[,jj]),sd(mcmc.result1$alpha_sample[,jj]),sd(mcmc.result1$beta_sample[,jj]))
  write.table(x = rbind(colMeans(cbind(mcmc.result1$A_sample[,jj],mcmc.result1$B_sample[,jj],mcmc.result1$alpha_sample[,jj],mcmc.result1$beta_sample[,jj])), sd_list), 
              file = paste("mean_std_post_samples_",jj,"th_cell_MBI_multiple_cells.csv", sep = ""), sep = ",", col.names = c("lambda_b", "lambda_d", "n", "t_r"), row.names = FALSE)
}
write.table(x = mcmc.result1$MCMC_sample_hyp, file = "post_samples_hyperparam_MBI.csv", sep = ",", col.names = c("alpha_lambda_b", "beta_lambda_b", "alpha_lambda_d", "beta_lambda_d", "alpha_n", "beta_n", "alpha_t_r", "beta_t_r"), row.names = FALSE)
hyp_sd_list = c(sd(mcmc.result1$MCMC_sample_hyp[,1]),sd(mcmc.result1$MCMC_sample_hyp[,2]),sd(mcmc.result1$MCMC_sample_hyp[,3]),sd(mcmc.result1$MCMC_sample_hyp[,4]),
                sd(mcmc.result1$MCMC_sample_hyp[,5]),sd(mcmc.result1$MCMC_sample_hyp[,6]),sd(mcmc.result1$MCMC_sample_hyp[,7]),sd(mcmc.result1$MCMC_sample_hyp[,8]))
write.table(x = rbind(colMeans(mcmc.result1$MCMC_sample_hyp), hyp_sd_list), file = "mean_std_post_samples_hyperparam_MBI.csv", sep = ",", col.names = c("alpha_lambda_b", "beta_lambda_b", "alpha_lambda_d", "beta_lambda_d", "alpha_n", "beta_n", "alpha_t_r", "beta_t_r"), row.names = FALSE)

AR_param = matrix(NA, nrow = num_trj, ncol  =4)
AR_param[,1] = colMeans(mcmc.result1$A_update)
AR_param[,2] = colMeans(mcmc.result1$B_update)
AR_param[,3] = colMeans(mcmc.result1$alpha_update)
AR_param[,4] = colMeans(mcmc.result1$beta_update)
AR_hyper = colMeans(mcmc.result1$update_matrix_hyp)

# For detail description of the format of output files. Please see the README file in our package.


