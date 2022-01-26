setwd("~/Dropbox/Time_delay/SciAdv_revision/code-for-revision/Code_resubmit_S1") # set the working directory contains the code file named "Queueing_functions.R"
source('Queueing_functions.R')

set.seed(123) # set the random seed

# load the input data for estimation
input_matrix_pre = as.matrix(read.table("input_data_feedback_example.csv", sep = ",", header = TRUE))

num_trj = dim(input_matrix_pre)[2] - 1 # the number of trajectories used for the estimation
input_matrix = pmax(input_matrix_pre[,2:(num_trj+1)] - input_matrix_pre[1,2:(num_trj+1)], 0)
timespan_pre = input_matrix_pre[,1] - input_matrix_pre[1,1]

basal_epsilon = 0.001 # epsilon value to avoid zero variance of the likelihood function

hyper_prior_mean = rep(1, 12)
hyper_prior_var = rep(10^6, 12)

theta_init = c(20, 0.1, 3, 1, 20, 1);
theta_hyper_init = rep(1, 12)

## Iteration numbers ====================
burn = 0 # the length of the burn-in period, i.e., the number of first samples which are discarded to reduce the initial condition dependency of MCMC method.
jump = 1 # the reciprocal of the thinning rate. If jump = 10, then every 10 samples from after the burn-in period are chosen as posterior samples. 
effnum = 110 # the effective number of iterations of the MCMC method.
selrow0 = seq(from = burn + jump, by = jump, to = burn + jump * effnum)


## Proposal variances ===================
prop_var_lambda_b = 100
prop_var_lambda_d = 0.0001
prop_var_n = 0.3
prop_var_tr = 0.3
prop_var_R = 100
prop_var_tau2 = 0.3
hyper_prop_var = c(3, 0.005, 3, 0.00001, 1, 0.1, 0.5, 0.003, 3, 0.005, 3, 0.01)


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
mcmc.result1 <- MCMC_function_NFL_ME_delay_feed_share_alpha_fine(data = inputdata1, timespan = timespan, var_list = var_list1, hyper_prior_mean = hyper_prior_mean, hyper_prior_var = hyper_prior_var, 
                                                 theta_init = theta_init, theta_hyper_init = theta_hyper_init, burn = burn, jump = jump, effnum = effnum, 
                                                 prop_var = c(prop_var_lambda_b, prop_var_lambda_d, prop_var_n, prop_var_tr, prop_var_R, prop_var_tau2), 
                                                 hyper_prop_var = hyper_prop_var, save.file.name = "feedback_multiple_cells_same_n_test", summary.file.name = "feedback_multiple_cells_same_n_test")

for(jj in 1:num_trj){
  write.table(x = cbind(mcmc.result1$A_sample[selrow0,jj],mcmc.result1$B_sample[selrow0,jj],mcmc.result1$alpha_sample[selrow0,jj],mcmc.result1$beta_sample[selrow0,jj],mcmc.result1$AR_sample[selrow0,jj],mcmc.result1$tau2_sample[selrow0,jj]), 
              file = paste("post_samples_",jj,"th_cell_MBI_multiple_cells_same_n_feedback.csv", sep = ""), sep = ",", col.names = c("lambda_b", "lambda_d", "n", "t_r", "R", "tau2"), row.names = FALSE)
  sd_list = c(sd(mcmc.result1$A_sample[selrow0,jj]),sd(mcmc.result1$B_sample[selrow0,jj]),sd(mcmc.result1$alpha_sample[selrow0,jj]),sd(mcmc.result1$beta_sample[selrow0,jj]),
              sd(mcmc.result1$AR_sample[selrow0,jj]),sd(mcmc.result1$tau2_sample[selrow0,jj]))
  write.table(x = rbind(colMeans(cbind(mcmc.result1$A_sample[selrow0,jj],mcmc.result1$B_sample[selrow0,jj],mcmc.result1$alpha_sample[selrow0,jj],mcmc.result1$beta_sample[selrow0,jj],mcmc.result1$AR_sample[selrow0,jj],mcmc.result1$tau2_sample[selrow0,jj])), sd_list), 
              file = paste("mean_std_post_samples_",jj,"th_cell_MBI_multiple_cells_same_n_feedback.csv", sep = ""), sep = ",", col.names = c("lambda_b", "lambda_d", "n", "t_r", "R", "tau2"), row.names = FALSE)
}

post_sample_hyp = mcmc.result1$MCMC_sample_hyp[selrow0,]



write.table(x = mcmc.result1$MCMC_sample_hyp, file = "post_samples_hyperparam_MBI_same_n_feedback.csv", sep = ",", 
            col.names = c("alpha_lambda_b", "beta_lambda_b", "alpha_lambda_d", "beta_lambda_d", 
                          "alpha_n", "beta_n", "alpha_t_r", "beta_t_r", "alpha_R", "beta_R", "alpha_tau2", "beta_tau2"), row.names = FALSE)


hyp_sd_list = c(sd(post_sample_hyp[,1]),sd(post_sample_hyp[,2]),sd(post_sample_hyp[,3]),sd(post_sample_hyp[,4]),
                sd(post_sample_hyp[,5]),sd(post_sample_hyp[,6]),sd(post_sample_hyp[,7]),sd(post_sample_hyp[,8]),
                sd(post_sample_hyp[,9]),sd(post_sample_hyp[,10]),sd(post_sample_hyp[,11]),sd(post_sample_hyp[,12]))

write.table(x = rbind(colMeans(post_sample_hyp), hyp_sd_list), file = "mean_std_post_samples_hyperparam_MBI_same_n_feedback.csv", sep = ",", 
            col.names = c("alpha_lambda_b", "beta_lambda_b", "alpha_lambda_d", "beta_lambda_d", "alpha_n", "beta_n", "alpha_t_r",
                          "beta_t_r", "alpha_R", "beta_R", "alpha_tau2", "beta_tau2"), 
            row.names = FALSE)

AR_lambda_b = colMeans(mcmc.result1$A_update)
AR_lambda_d = colMeans(mcmc.result1$B_update)
AR_n = colMeans(mcmc.result1$alpha_update)
AR_tr = colMeans(mcmc.result1$beta_update)
AR_R = colMeans(mcmc.result1$AR_update)
AR_tau2 = colMeans(mcmc.result1$tau2_update)
AR_hyper = colMeans(mcmc.result1$update_matrix_hyp)

# For detail description of the format of output files. Please see the README file in our package.

