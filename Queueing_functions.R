# library(ramcmc) 
# library(mvtnorm)
library(MASS)

mean_trajectory <- function(timespan, theta){
  length_of_T = length(timespan)
  mean_trj = rep(0, length_of_T)
  integral_val_partition = rep(0, length_of_T)
  integral_val = rep(0, length_of_T)
  
  integrand2 <- function(tau, theta){
    return(pgamma(tau, shape = theta[3], scale = theta[4]) * exp(theta[2]* tau))
  }
  for (ii in 1:length_of_T){
    if (ii == 1){
      integral_val_partition[ii] = integrate(integrand2, lower=0, upper=timespan[ii], theta = theta, abs.tol = 1e-15)$value
    }
    else{
      integral_val_partition[ii] = integrate(integrand2, lower=timespan[ii-1], upper=timespan[ii], theta = theta, abs.tol = 1e-15)$value
    }
  }
  integral_val <- cumsum(integral_val_partition)
  
  for (ii in 1:length_of_T){
    t = timespan[ii];
    # cat(ii); cat("\n");
    mean_trj[ii] = theta[1] * exp(-theta[2] * t) * integral_val[ii]
  }
  return(mean_trj)
}

log_of_gampdf <- function(x,shape,scale){
  logval = -log(gamma(shape-ceiling(shape)+1)) - shape * log(scale) + (shape-1) * log(x) - 1/scale * x
  for (d in 1:length(logval)){
    if ((ceiling(shape[d]) - 1) > 0){
      for (jj in 1:(ceiling(shape[d]) - 1)){
        logval[d] = logval[d] - log(shape[d] - jj)
      }
    }
  }
  return(logval)
}

log_of_normpdf <- function(x, mu, sigma){
  logval = -1/2 * log(2*pi) - log(sigma) - (x-mu)^2 / (2*sigma^2);
  return(logval)
}

MCMC_function <- function(data, timespan, var_list, prior_mean, prior_var, theta_init, burn, jump, effnum, prop_var = c(0.000001, 0.1, 0.003)){
  
  total_num = burn + jump * effnum;
  selrow = seq(from = burn + jump, by = jump, to = total_num)
  MCMC_sample = matrix(0, nrow = total_num, ncol = length(theta_init))
  update_matrix =  matrix(0, nrow = total_num, ncol = length(theta_init))
  
  
  
  MCMC_sample[1,] = theta_init
  update_matrix[2:total_num,1] = 1
  
  for (rep in 2:total_num){
    if (rep %% 1000 == 0){
      cat(rep); cat("\n");
      cat(summary.file.name); cat("\n");
    }
    # 
    mean_trj1 = mean_trajectory(timespan, c(1, MCMC_sample[rep-1, 2], MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4]))
    
    post_var = 1 / (sum(mean_trj1^2 / var_list) + prior_var[1]^(-1))
    post_mean = post_var * (sum(data * mean_trj1 / var_list) + prior_mean[1]/prior_var[1])
    
    MCMC_sample[rep, 1] = rnorm(1, mean = post_mean, sd = sqrt(post_var))
    
    # MCMC for theta(2): decay rate of X(t);
    # Random walk Metropolis algorithm (use normal proposal, symmetric)
    
    prop_mean2 = MCMC_sample[rep-1,2]
    prop_var2 = prop_var[1] # tune this parameter
    
    theta2_cand = rnorm(1, mean = prop_mean2, sd = sqrt(prop_var2))
    
    if (theta2_cand > 0){
      gamma_pri_alpha_2 = prior_mean[2]^2 / prior_var[2]
      gamma_pri_beta_2 = prior_var[2] / prior_mean[2]
      
      log_pri = log_of_gampdf(MCMC_sample[rep-1, 2], gamma_pri_alpha_2, gamma_pri_beta_2)
      log_pri_st = log_of_gampdf(theta2_cand, gamma_pri_alpha_2, gamma_pri_beta_2)
      
      mean_trj2 = mean_trajectory(timespan, c(MCMC_sample[rep,1], MCMC_sample[rep-1, 2], MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4]))
      mean_trj2_st = mean_trajectory(timespan, c(MCMC_sample[rep,1], theta2_cand, MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4]))
      
      log_lik = sum(log_of_normpdf(data, mean_trj2, sqrt(var_list)))
      log_lik_st = sum(log_of_normpdf(data, mean_trj2_st, sqrt(var_list)))
      
      acceptance_ratio2 = exp(log_pri_st - log_pri + log_lik_st - log_lik)
      
      if (runif(1) < acceptance_ratio2){
        MCMC_sample[rep, 2] = theta2_cand
        update_matrix[rep,2] = 1
      }else{
        MCMC_sample[rep, 2] = MCMC_sample[rep-1,2]
      }
    }else{
      MCMC_sample[rep, 2] = MCMC_sample[rep-1, 2]
    }
    
    # MCMC for theta(3): the shape parameter of delay distribution
    
    # MCMC_sample(rep,3) = MCMC_sample(rep-1,3);
    
    prop_mean3 = MCMC_sample[rep-1,3];
    prop_var3 = prop_var[2]; # tune this parameter
    
    theta3_cand = rnorm(1, mean = prop_mean3, sd = sqrt(prop_var3));
    if (theta3_cand > 0){
      gamma_pri_alpha_3 = prior_mean[3]^2 / prior_var[3]
      gamma_pri_beta_3 = prior_var[3]/prior_mean[3] 
      
      log_pri = log_of_gampdf(MCMC_sample[rep-1, 3], gamma_pri_alpha_3, gamma_pri_beta_3)
      log_pri_st = log_of_gampdf(theta3_cand, gamma_pri_alpha_3, gamma_pri_beta_3)
      
      mean_trj3 = mean_trajectory(timespan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4]))
      mean_trj3_st = mean_trajectory(timespan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], theta3_cand, MCMC_sample[rep-1, 4]))
      
      log_lik = sum(log_of_normpdf(data, mean_trj3, sqrt(var_list)))
      log_lik_st = sum(log_of_normpdf(data, mean_trj3_st, sqrt(var_list)))
      
      acceptance_ratio3 = exp(log_pri_st - log_pri + log_lik_st - log_lik)
      
      if (runif(1) < acceptance_ratio3){
        MCMC_sample[rep, 3] = theta3_cand
        update_matrix[rep, 3] = 1
      }else{
        MCMC_sample[rep, 3] = MCMC_sample[rep-1,3]
      }
    }else{
      MCMC_sample[rep, 3] = MCMC_sample[rep-1, 3]
    }
    
    # === MCMC for theta(4): the rate parameter of delay distribution ===
    
    #     MCMC_sample(rep,4) = MCMC_sample(rep-1,4);
    
    
    prop_mean4 = MCMC_sample[rep-1,4];
    prop_var4 = prop_var[3]; # tune this parameter
    
    theta4_cand = rnorm(1, mean = prop_mean4, sd = sqrt(prop_var4))
    if (theta4_cand > 0){
      gamma_pri_alpha_4 = prior_mean[4]^2 / prior_var[4]
      gamma_pri_beta_4 = prior_var[4]/prior_mean[4]
      
      log_pri = log_of_gampdf(MCMC_sample[rep-1, 4], gamma_pri_alpha_4, gamma_pri_beta_4)
      log_pri_st = log_of_gampdf(theta4_cand, gamma_pri_alpha_4, gamma_pri_beta_4)
      
      mean_trj4 = mean_trajectory(timespan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep, 3], MCMC_sample[rep-1, 4]))
      mean_trj4_st = mean_trajectory(timespan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep, 3], theta4_cand))
      
      log_lik = sum(log_of_normpdf(data, mean_trj4, sqrt(var_list)))
      log_lik_st = sum(log_of_normpdf(data, mean_trj4_st, sqrt(var_list)))
      
      acceptance_ratio4 = exp(log_pri_st - log_pri + log_lik_st - log_lik)
      
      if (runif(1) < acceptance_ratio4){
        MCMC_sample[rep, 4] = theta4_cand
        update_matrix[rep, 4] = 1
      }else{
        MCMC_sample[rep, 4] = MCMC_sample[rep-1, 4]
      }
    }else{
      MCMC_sample[rep, 4] = MCMC_sample[rep-1, 4]
    }
  }
  results = MCMC_sample[selrow,]
  acceptance = update_matrix
  my_list <- list("results" = results, "acceptance" = acceptance)
  return(my_list)
}


MCMC_function_fixB <- function(data, timespan, var_list, prior_mean, prior_var, theta_init, burn, jump, effnum, prop_var = c(0.000001, 0.1, 0.003), B){
  
  total_num = burn + jump * effnum;
  MCMC_sample = matrix(0, nrow = total_num, ncol = length(theta_init))
  update_matrix =  matrix(0, nrow = total_num, ncol = length(theta_init))
  
  
  
  MCMC_sample[1,] = theta_init
  MCMC_sample[1,2] = B
  update_matrix[2:total_num,1:2] = 1
  
  for (rep in 2:total_num){
    #  if (rep %% 10 == 0){
    # cat(rep); cat("\n");
    # }
    
    mean_trj1 = mean_trajectory(timespan, c(1, MCMC_sample[rep-1, 2], MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4]))
    
    post_var = 1 / (sum(mean_trj1^2 / var_list) + prior_var[1]^(-1))
    post_mean = post_var * (sum(data * mean_trj1 / var_list) + prior_mean[1]/prior_var[1])
    
    MCMC_sample[rep, 1] = rnorm(1, mean = post_mean, sd = sqrt(post_var))
    
    # MCMC for theta(2): decay rate of X(t);
    # Random walk Metropolis algorithm (use normal proposal, symmetric)
    
    prop_mean2 = MCMC_sample[rep-1,2]
    prop_var2 = prop_var[1] # tune this parameter
    
    # theta2_cand = rnorm(1, mean = prop_mean2, sd = sqrt(prop_var2))
    
    # if (theta2_cand > 0){
    #   gamma_pri_alpha_2 = prior_mean[2]^2 / prior_var[2]
    #   gamma_pri_beta_2 = prior_var[2] / prior_mean[2]
    #   
    #   log_pri = log_of_gampdf(MCMC_sample[rep-1, 2], gamma_pri_alpha_2, gamma_pri_beta_2)
    #   log_pri_st = log_of_gampdf(theta2_cand, gamma_pri_alpha_2, gamma_pri_beta_2)
    #   
    #   mean_trj2 = mean_trajectory(timespan, c(MCMC_sample[rep,1], MCMC_sample[rep-1, 2], MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4]))
    #   mean_trj2_st = mean_trajectory(timespan, c(MCMC_sample[rep,1], theta2_cand, MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4]))
    #   
    #   log_lik = sum(log_of_normpdf(data, mean_trj2, sqrt(var_list)))
    #   log_lik_st = sum(log_of_normpdf(data, mean_trj2_st, sqrt(var_list)))
    #   
    #   acceptance_ratio2 = exp(log_pri_st - log_pri + log_lik_st - log_lik)
    #   
    #   if (runif(1) < acceptance_ratio2){
    #     MCMC_sample[rep, 2] = theta2_cand
    #     update_matrix[rep,2] = 1
    #   }else{
    #     MCMC_sample[rep, 2] = MCMC_sample[rep-1,2]
    #   }
    # }else{
    #   MCMC_sample[rep, 2] = MCMC_sample[rep-1, 2]
    # }
    MCMC_sample[rep, 2] = MCMC_sample[rep-1, 2]
    # MCMC for theta(3): the shape parameter of delay distribution
    
    # MCMC_sample(rep,3) = MCMC_sample(rep-1,3);
    
    prop_mean3 = MCMC_sample[rep-1,3];
    prop_var3 = prop_var[2]; # tune this parameter
    
    theta3_cand = rnorm(1, mean = prop_mean3, sd = sqrt(prop_var3));
    if (theta3_cand > 0){
      gamma_pri_alpha_3 = prior_mean[3]^2 / prior_var[3]
      gamma_pri_beta_3 = prior_var[3]/prior_mean[3] 
      
      log_pri = log_of_gampdf(MCMC_sample[rep-1, 3], gamma_pri_alpha_3, gamma_pri_beta_3)
      log_pri_st = log_of_gampdf(theta3_cand, gamma_pri_alpha_3, gamma_pri_beta_3)
      
      mean_trj3 = mean_trajectory(timespan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4]))
      mean_trj3_st = mean_trajectory(timespan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], theta3_cand, MCMC_sample[rep-1, 4]))
      
      log_lik = sum(log_of_normpdf(data, mean_trj3, sqrt(var_list)))
      log_lik_st = sum(log_of_normpdf(data, mean_trj3_st, sqrt(var_list)))
      
      acceptance_ratio3 = exp(log_pri_st - log_pri + log_lik_st - log_lik)
      
      if (runif(1) < acceptance_ratio3){
        MCMC_sample[rep, 3] = theta3_cand
        update_matrix[rep, 3] = 1
      }else{
        MCMC_sample[rep, 3] = MCMC_sample[rep-1,3]
      }
    }else{
      MCMC_sample[rep, 3] = MCMC_sample[rep-1, 3]
    }
    
    # === MCMC for theta(4): the rate parameter of delay distribution ===
    
    #     MCMC_sample(rep,4) = MCMC_sample(rep-1,4);
    
    
    prop_mean4 = MCMC_sample[rep-1,4];
    prop_var4 = prop_var[3]; # tune this parameter
    
    theta4_cand = rnorm(1, mean = prop_mean4, sd = sqrt(prop_var4))
    if (theta4_cand > 0){
      gamma_pri_alpha_4 = prior_mean[4]^2 / prior_var[4]
      gamma_pri_beta_4 = prior_var[4]/prior_mean[4]
      
      log_pri = log_of_gampdf(MCMC_sample[rep-1, 4], gamma_pri_alpha_4, gamma_pri_beta_4)
      log_pri_st = log_of_gampdf(theta4_cand, gamma_pri_alpha_4, gamma_pri_beta_4)
      
      mean_trj4 = mean_trajectory(timespan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep, 3], MCMC_sample[rep-1, 4]))
      mean_trj4_st = mean_trajectory(timespan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep, 3], theta4_cand))
      
      log_lik = sum(log_of_normpdf(data, mean_trj4, sqrt(var_list)))
      log_lik_st = sum(log_of_normpdf(data, mean_trj4_st, sqrt(var_list)))
      
      acceptance_ratio4 = exp(log_pri_st - log_pri + log_lik_st - log_lik)
      
      if (runif(1) < acceptance_ratio4){
        MCMC_sample[rep, 4] = theta4_cand
        update_matrix[rep, 4] = 1
      }else{
        MCMC_sample[rep, 4] = MCMC_sample[rep-1, 4]
      }
    }else{
      MCMC_sample[rep, 4] = MCMC_sample[rep-1, 4]
    }
  }
  results = MCMC_sample;
  acceptance = update_matrix;
  my_list <- list("results" = results, "acceptance" = acceptance)
  return(my_list)
}


MCMC_function_ME <- function(data, timespan, var_list, hyper_prior_mean = rep(1,8), hyper_prior_var = rep(1000,8), 
                             theta_init, theta_hyper_init, burn = 0, jump = 1, effnum = 100, 
                             prop_var = c(50, 0.000001, 0.1, 0.003), hyper_prop_var = c(50,1,0.001,1,1,1,0.01,1), flatpri = FALSE, save.file.name = "dummy1", summary.file.name = "dummy2"){
  
  # In this function, data shold be a matrix with more than one column. Each column represents a timeseries data.
  # var_list is also a matrix with the size identical to data.
  
  num_repeat = burn + jump * effnum;
  selrow = seq(from = burn + jump, by = jump, to = num_repeat)
  num_data = dim(data)[2]
  num_param = length(theta_init)
  
  MSE_abs_cumul = c()
  MSE_rel_cumul = c()
  MSE_first_rel_cumul = c()
  
  Prod_posterior_mean_cumul = c()
  death_posterior_mean_cumul = c()
  alpha_posterior_mean_cumul = c()    
  beta_posterior_mean_cumul = c()
  timedelay_posterior_mean_cumul = c()
  
  Prod_cv_list_cumul = c()
  death_cv_list_cumul = c()
  alpha_cv_list_cumul = c()
  beta_cv_list_cumul = c()
  Meandelay_cv_list_cumul = c()
  
  birth_acf_lag5_cumul = c()
  death_acf_lag5_cumul = c()
  alpha_acf_lag5_cumul = c()
  beta_acf_lag5_cumul = c()
  meandelay_acf_lag5_cumul = c()
  birth_acf_lag10_cumul = c()
  death_acf_lag10_cumul = c()
  alpha_acf_lag10_cumul = c()
  beta_acf_lag10_cumul = c()
  meandelay_acf_lag10_cumul = c()
  
  
  A_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  B_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  beta_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  
  A_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  B_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  beta_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  
  MCMC_sample_hyp = matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  update_matrix_hyp =  matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  
  A_sample[1,] = theta_init[1]
  B_sample[1,] = theta_init[2]
  alpha_sample[1,] = theta_init[3]
  beta_sample[1,] = theta_init[4]
  
  MCMC_sample_hyp[1,] = theta_hyper_init
  
  thetaCurrent <- matrix(NA, nrow = num_data, ncol = 4)
  thetaCurrent[,1] <- theta_init[1]
  thetaCurrent[,2] <- theta_init[2]
  thetaCurrent[,3] <- theta_init[3]
  thetaCurrent[,4] <- theta_init[4]
  
  thetaCand <- thetaCurrent
  
  thetaHypCurrent <- theta_hyper_init
  thetaHypCand <- thetaHypCurrent
  
  for (rep in 2:num_repeat){
    if (rep %% 1000 == 0){
      cat(rep); cat("\n");
      cat(summary.file.name); cat("\n");
    }
    ## Start to sample the parameters for each individual
    for (dataid in 1:num_data){
      for (param_id in 1:num_param){
        param_cand <- rnorm(1, mean = thetaCurrent[dataid,param_id], sd = sqrt(prop_var[param_id]))
        if (param_cand > 0){
          thetaCand[dataid, param_id] = param_cand
          mean_trj_current <- mean_trajectory(timespan, thetaCurrent[dataid,])
          mean_trj_cand <- mean_trajectory(timespan, thetaCand[dataid,])
          log_pri = log_of_gampdf(thetaCurrent[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          log_pri_cand = log_of_gampdf(thetaCand[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          
          log_lik = sum(log_of_normpdf(data[,dataid], mean_trj_current, sqrt(var_list[,dataid])))
          log_lik_cand = sum(log_of_normpdf(data[,dataid], mean_trj_cand, sqrt(var_list[,dataid])))
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaCurrent = thetaCand
            switch(param_id, 
                   A_update[rep,dataid] <- 1, 
                   B_update[rep,dataid] <- 1, 
                   alpha_update[rep,dataid] <- 1, 
                   beta_update[rep,dataid] <- 1)
          }else{
            thetaCand = thetaCurrent
          }
        }
      }
    }
    
    ## Start to sample the hyper parameters for all populations.
    for (hyp_param_id in 1:(2*num_param)){
      hyp_param_cand <- rnorm(1, mean = thetaHypCurrent[hyp_param_id], sd = sqrt(hyper_prop_var[hyp_param_id]))
      
      if (hyp_param_cand > 0){
        thetaHypCand[hyp_param_id] = hyp_param_cand
        
        gamma_pri_alpha = hyper_prior_mean[hyp_param_id]^2 / hyper_prior_var[hyp_param_id]
        gamma_pri_beta = hyper_prior_var[hyp_param_id]/hyper_prior_mean[hyp_param_id]
        
        log_pri = log_of_gampdf(thetaHypCurrent[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
        log_pri_cand = log_of_gampdf(thetaHypCand[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
        
        if (hyp_param_id %% 2 == 1){
          log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id], num_data), rep(thetaHypCurrent[hyp_param_id+1], num_data)))
          log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id], num_data), rep(thetaHypCand[hyp_param_id+1], num_data)))
        }else{
          log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id-1], num_data), rep(thetaHypCurrent[hyp_param_id], num_data)))
          log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id-1], num_data), rep(thetaHypCand[hyp_param_id], num_data)))
        }
        
        acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
        
        if (runif(1) < acceptance_ratio){
          thetaHypCurrent = thetaHypCand
          update_matrix_hyp[rep, hyp_param_id] = 1
        }else{
          thetaHypCand = thetaHypCurrent
        }
      }
    }
    A_sample[rep,] = thetaCurrent[,1]
    B_sample[rep,] = thetaCurrent[,2]
    alpha_sample[rep,] = thetaCurrent[,3]
    beta_sample[rep,] = thetaCurrent[,4]
    MCMC_sample_hyp[rep, ] = thetaHypCurrent
    
    
    if (rep %% 55000 == 0 | rep == num_repeat){
      # mcmc.result1 <- list("A_sample" = A_sample, "B_sample" = B_sample, "alpha_sample" = alpha_sample, "beta_sample" = beta_sample,
      #                      "A_update" = A_update, "B_update" = B_update, "alpha_update" = alpha_update, "beta_update" = beta_update,
      #                      "MCMC_sample_hyp" = MCMC_sample_hyp, "update_matrix_hyp" = update_matrix_hyp, "rep" = rep)
      # inputdata1 = data
      # save(mcmc.result1, inputdata1, file = paste(save.file.name, ".Rdata", sep = ""))
    }
    
    
    # if (rep %% 200 == 0  | rep == num_repeat){
    #   cat(paste("Current iteration: ",rep, sep =""));
    #   cat("\n");
    #   selrow = seq(rep/2, by = 100, to = rep)
    #   MSE_abs = rep(NA, num_data)
    #   MSE_rel = rep(NA, num_data)
    #   MSE_first_rel = rep(NA, num_data)
    #   timedelay_posterior_mean = colMeans(alpha_sample[selrow,]*beta_sample[selrow,])
    #   Prod_cv_list = rep(NA, num_data)
    #   death_cv_list = rep(NA, num_data)
    #   alpha_cv_list = rep(NA, num_data)
    #   beta_cv_list = rep(NA, num_data)
    #   Meandelay_cv_list = rep(NA, num_data)
    #   
    #   birth_acf_lag5 = rep(NA, num_data)
    #   death_acf_lag5 = rep(NA, num_data)
    #   alpha_acf_lag5 = rep(NA, num_data)
    #   beta_acf_lag5 = rep(NA, num_data)
    #   meandelay_acf_lag5 = rep(NA, num_data)
    #   birth_acf_lag10 = rep(NA, num_data)
    #   death_acf_lag10 = rep(NA, num_data)
    #   alpha_acf_lag10 = rep(NA, num_data)
    #   beta_acf_lag10 = rep(NA, num_data)
    #   meandelay_acf_lag10 = rep(NA, num_data)
    #   first_interval = which(timespan < (2.5*mean(timedelay_posterior_mean)))
    #   
    #   mean_trj_matrix = matrix(0, nrow = num_data, ncol = length(timespan) + 1)
    #   
    #   for (dataid in 1:num_data){
    #     mean_trj_compare <- mean_trajectory(timespan, c(mean(A_sample[selrow,dataid]), mean(B_sample[selrow,dataid]), mean(alpha_sample[selrow,dataid]),mean(beta_sample[selrow,dataid])))
    #     mean_trj_matrix[dataid, 2:(length(timespan)+1)] = mean_trj_compare
    #     
    #     MSE_abs[dataid] = (sum((data[,dataid] - mean_trj_compare)^2))^0.5/length(timespan)
    #     MSE_rel[dataid] = MSE_abs[dataid] / data[length(timespan), dataid]
    #     if (length(timespan[first_interval]) == 0){
    #       MSE_first_rel[dataid] = 0  
    #     }else{
    #       MSE_first_rel[dataid] = (sum((data[first_interval,dataid] - mean_trj_compare[first_interval])^2))^0.5 / length(timespan[first_interval])
    #       MSE_first_rel[dataid] = MSE_first_rel[dataid] / data[length(timespan[first_interval]), dataid]
    #     }
    #     Prod_cv_list[dataid] = sd(A_sample[selrow, dataid]) / mean(A_sample[selrow, dataid])
    #     death_cv_list[dataid] = sd(B_sample[selrow, dataid]) / mean(B_sample[selrow, dataid])
    #     alpha_cv_list[dataid] = sd(alpha_sample[selrow, dataid]) / mean(alpha_sample[selrow, dataid])
    #     beta_cv_list[dataid] = sd(beta_sample[selrow, dataid]) / mean(beta_sample[selrow, dataid])
    #     Meandelay_cv_list[dataid] = sd(alpha_sample[selrow, dataid]*beta_sample[selrow, dataid]) / mean(alpha_sample[selrow, dataid]*beta_sample[selrow, dataid])
    #     
    #     alpha_acf_lag5[dataid] = acf(alpha_sample[selrow, dataid], plot = FALSE)$acf[6]
    #     alpha_acf_lag10[dataid] = acf(alpha_sample[selrow, dataid], plot = FALSE)$acf[11]
    #     beta_acf_lag5[dataid] = acf(beta_sample[selrow, dataid], plot = FALSE)$acf[6]
    #     beta_acf_lag10[dataid] = acf(beta_sample[selrow, dataid], plot = FALSE)$acf[11]
    #     birth_acf_lag5[dataid] = acf(A_sample[selrow, dataid], plot = FALSE)$acf[6]
    #     birth_acf_lag10[dataid] = acf(A_sample[selrow, dataid], plot = FALSE)$acf[11]
    #     death_acf_lag5[dataid] = acf(B_sample[selrow, dataid], plot = FALSE)$acf[6]
    #     death_acf_lag10[dataid] = acf(B_sample[selrow, dataid], plot = FALSE)$acf[11]
    #     meandelay_acf_lag5[dataid] = acf(alpha_sample[selrow, dataid]*beta_sample[selrow, dataid], plot = FALSE)$acf[6]
    #     meandelay_acf_lag10[dataid] = acf(alpha_sample[selrow, dataid]*beta_sample[selrow, dataid], plot = FALSE)$acf[11]
    #   }
    #   
    #   cat(paste("== ", summary.file.name ," = \n", sep = ""));
    #   cat(paste("== constant thinning rate == \n", sep = ""));
    #   # cat(paste("MSE abs mean:", round(mean(MSE_abs),3), ",  max:", round(max(MSE_abs),3), "\n", sep = ""))
    #   cat(paste("MSE rel mean:", round(mean(MSE_rel), 5), ",  max:", round(max(MSE_rel), 5), "\n", sep = ""))
    #   cat(paste("MSE first rel mean:", round(mean(MSE_first_rel), 5), ",  max:", round(max(MSE_first_rel), 5), "\n", sep = ""))
    #   cat(paste("   Timedelay \n", sep = ""))
    #   cat(paste("Mean:", round(mean(timedelay_posterior_mean),3), ",  CV:", round(sd(timedelay_posterior_mean)/mean(timedelay_posterior_mean),3), 
    #             "  CV's mean:", round(mean(Meandelay_cv_list),3), ",  CV's max:", round(max(Meandelay_cv_list),3), 
    #             "\n", sep = ""))
    #   cat(paste("   Production \n", sep = ""))
    #   cat(paste("Mean:", round(mean(colMeans(A_sample[selrow,])),3), ",  CV:", round(sd(colMeans(A_sample[selrow,]))/mean(colMeans(A_sample[selrow,])),3), 
    #             "  CV's mean:", round(mean(Prod_cv_list),3), ",  CV's max:", round(max(Prod_cv_list),3), 
    #             "\n", sep = ""))
    #   cat(paste("   Alpha \n", sep = ""))
    #   cat(paste("Mean:", round(mean(colMeans(alpha_sample[selrow,])),3), ",  CV:", round(sd(colMeans(alpha_sample[selrow,]))/mean(colMeans(alpha_sample[selrow,])),3), "\n", sep = ""))
    #   if(length(acf(alpha_sample[selrow,1], plot = FALSE)$acf) >= 11){
    #     cat(paste("ACF lag 5:", round(acf(alpha_sample[selrow,1], plot = FALSE)$acf[6],3), ",  ACF lag 10:", round(acf(alpha_sample[selrow,1], plot = FALSE)$acf[11],3), "\n", sep = ""))
    #   }else if(length(acf(alpha_sample[selrow,1], plot = FALSE)$acf) >= 6){
    #     cat(paste("ACF lag 5:", round(acf(alpha_sample[selrow,1], plot = FALSE)$acf[6],3), "\n", sep = ""))
    #   }
    #   cat("\n");
    #   
    #   
    #   ## ======= save all posterior samples ============== ## 
    #   MSE_abs_cumul = rbind(MSE_abs_cumul, MSE_abs)
    #   MSE_rel_cumul = rbind(MSE_rel_cumul, MSE_rel)
    #   MSE_first_rel_cumul = rbind(MSE_first_rel_cumul, MSE_first_rel)
    #   
    #   Prod_posterior_mean_cumul = rbind(Prod_posterior_mean_cumul, colMeans(A_sample[selrow,]))
    #   death_posterior_mean_cumul = rbind(death_posterior_mean_cumul, colMeans(B_sample[selrow,]))
    #   alpha_posterior_mean_cumul = rbind(alpha_posterior_mean_cumul, colMeans(alpha_sample[selrow,]))
    #   beta_posterior_mean_cumul = rbind(beta_posterior_mean_cumul, colMeans(beta_sample[selrow,]))
    #   timedelay_posterior_mean_cumul = rbind(timedelay_posterior_mean_cumul,timedelay_posterior_mean)
    #   
    #   
    #   Prod_cv_list_cumul = rbind(Prod_cv_list_cumul, Prod_cv_list)
    #   death_cv_list_cumul = rbind(death_cv_list_cumul, death_cv_list)
    #   alpha_cv_list_cumul = rbind(alpha_cv_list_cumul, alpha_cv_list)
    #   beta_cv_list_cumul = rbind(beta_cv_list_cumul, beta_cv_list)
    #   Meandelay_cv_list_cumul = rbind(Meandelay_cv_list_cumul, Meandelay_cv_list)
    #   
    #   
    #   birth_acf_lag5_cumul = rbind(birth_acf_lag5_cumul, birth_acf_lag5)
    #   birth_acf_lag10_cumul = rbind(birth_acf_lag10_cumul, birth_acf_lag10)
    #   death_acf_lag5_cumul = rbind(death_acf_lag5_cumul, death_acf_lag5)
    #   death_acf_lag10_cumul = rbind(death_acf_lag10_cumul, death_acf_lag10)
    #   alpha_acf_lag5_cumul = rbind(alpha_acf_lag5_cumul, alpha_acf_lag5)
    #   alpha_acf_lag10_cumul = rbind(alpha_acf_lag10_cumul, alpha_acf_lag10)
    #   beta_acf_lag5_cumul = rbind(beta_acf_lag5_cumul, beta_acf_lag5)
    #   beta_acf_lag10_cumul = rbind(beta_acf_lag10_cumul, beta_acf_lag10)
    #   meandelay_acf_lag5_cumul = rbind(meandelay_acf_lag5_cumul, meandelay_acf_lag5)
    #   meandelay_acf_lag10_cumul = rbind(meandelay_acf_lag10_cumul, meandelay_acf_lag10)
    #   
    #   # write.table(A_sample[1:rep,], file = paste(save.file.name, "_post_sample_birth.csv", sep =""), sep = ",", col.names = FALSE, row.names = FALSE)
    #   # write.table(B_sample[1:rep,], file = paste(save.file.name, "_post_sample_death.csv", sep =""), sep = ",", col.names = FALSE, row.names = FALSE)
    #   # write.table(alpha_sample[1:rep,], file = paste(save.file.name, "_post_sample_alpha.csv", sep =""), sep = ",", col.names = FALSE, row.names = FALSE)
    #   # write.table(beta_sample[1:rep,], file = paste(save.file.name, "_post_sample_beta.csv", sep =""), sep = ",", col.names = FALSE, row.names = FALSE)
    #   # write.table(MCMC_sample_hyp[1:rep,], file = paste(save.file.name, "_post_sample_hyper.csv", sep =""), sep = ",", col.names = FALSE, row.names = FALSE)
    # }
  } 
  my_list <- list("A_sample" = A_sample[selrow,], "B_sample" = B_sample[selrow,], "alpha_sample" = alpha_sample[selrow,], "beta_sample" = beta_sample[selrow,],
                  "A_update" = A_update, "B_update" = B_update, "alpha_update" = alpha_update, "beta_update" = beta_update,
                  "MCMC_sample_hyp" = MCMC_sample_hyp[selrow,], "update_matrix_hyp" = update_matrix_hyp)
  return(my_list)
}

MCMC_function_ME_share_alpha <- function(data, timespan, var_list, hyper_prior_mean = rep(1,8), hyper_prior_var = rep(1000,8), 
                                         theta_init, theta_hyper_init, burn = 0, jump = 1, effnum = 100, 
                                         prop_var = c(50, 0.000001, 0.1, 0.003), hyper_prop_var = c(50,1,0.001,1,1,1,0.01,1), save.file.name = "dummy1", summary.file.name = "dummy2"){
  
  # In this function, data shold be a matrix with more than one column. Each column represents a timeseries data.
  # var_list is also a matrix with the size identical to data.
  
  num_repeat = burn + jump * effnum;
  selrow = seq(from = burn + jump, by = jump, to = num_repeat)
  
  num_data = dim(data)[2]
  num_param = length(theta_init)
  
  MSE_abs_cumul = c()
  MSE_rel_cumul = c()
  MSE_first_rel_cumul = c()
  
  Prod_posterior_mean_cumul = c()
  death_posterior_mean_cumul = c()
  alpha_posterior_mean_cumul = c()    
  beta_posterior_mean_cumul = c()
  timedelay_posterior_mean_cumul = c()
  
  Prod_cv_list_cumul = c()
  death_cv_list_cumul = c()
  alpha_cv_list_cumul = c()
  beta_cv_list_cumul = c()
  Meandelay_cv_list_cumul = c()
  
  birth_acf_lag5_cumul = c()
  death_acf_lag5_cumul = c()
  alpha_acf_lag5_cumul = c()
  beta_acf_lag5_cumul = c()
  meandelay_acf_lag5_cumul = c()
  birth_acf_lag10_cumul = c()
  death_acf_lag10_cumul = c()
  alpha_acf_lag10_cumul = c()
  beta_acf_lag10_cumul = c()
  meandelay_acf_lag10_cumul = c()
  
  A_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  B_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  beta_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  
  A_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  B_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  beta_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  
  MCMC_sample_hyp = matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  update_matrix_hyp =  matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  
  A_sample[1,] = theta_init[1]
  B_sample[1,] = theta_init[2]
  alpha_sample[1,] = theta_init[3]
  beta_sample[1,] = theta_init[4]
  
  MCMC_sample_hyp[1,] = theta_hyper_init
  
  thetaCurrent <- matrix(NA, nrow = num_data, ncol = 4)
  thetaCurrent[,1] <- theta_init[1]
  thetaCurrent[,2] <- theta_init[2]
  thetaCurrent[,3] <- theta_init[3]
  thetaCurrent[,4] <- theta_init[4]
  
  thetaCand <- thetaCurrent
  
  thetaHypCurrent <- theta_hyper_init
  thetaHypCand <- thetaHypCurrent
  
  for (rep in 2:num_repeat){
    if (rep %% 1000 == 0){
      cat(rep); cat("\n");
      cat(summary.file.name); cat("\n");
    }
    
    ## Start to sample the parameters for each individual
    for (dataid in 1:num_data){
      for (param_id in 1:num_param){
        param_cand <- rnorm(1, mean = thetaCurrent[dataid,param_id], sd = sqrt(prop_var[param_id]))
        if (param_cand > 0){
          if(param_id == 3){
            # thetaCand[, param_id] = param_cand
            next
          }else{
            thetaCand[dataid, param_id] = param_cand
          }  
          mean_trj_current <- mean_trajectory(timespan, thetaCurrent[dataid,])
          mean_trj_cand <- mean_trajectory(timespan, thetaCand[dataid,])
          log_pri = log_of_gampdf(thetaCurrent[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          log_pri_cand = log_of_gampdf(thetaCand[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          
          log_lik = sum(log_of_normpdf(data[,dataid], mean_trj_current, sqrt(var_list[,dataid])))
          log_lik_cand = sum(log_of_normpdf(data[,dataid], mean_trj_cand, sqrt(var_list[,dataid])))
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaCurrent = thetaCand
            switch(param_id, 
                   A_update[rep,dataid] <- 1, 
                   B_update[rep,dataid] <- 1, 
                   alpha_update[rep,dataid] <- 1, 
                   beta_update[rep,dataid] <- 1)
          }else{
            thetaCand = thetaCurrent
          }
        }
      }
    }
    
    ## Sample the shared alpha values for all individual
    
    param_cand <- rnorm(1, mean = thetaCurrent[1,3], sd = sqrt(prop_var[3]))
    
    if (param_cand > 0){
      thetaCand[, 3] = param_cand
      log_pri = log_of_gampdf(thetaCurrent[1,3], thetaHypCurrent[5], thetaHypCurrent[6])
      log_pri_cand = log_of_gampdf(thetaCand[1,3], thetaHypCurrent[5], thetaHypCurrent[6])
      log_lik = 0
      log_lik_cand = 0
      for (dataid in 1:num_data){
        mean_trj_current <- mean_trajectory(timespan, thetaCurrent[dataid,])
        mean_trj_cand <- mean_trajectory(timespan, thetaCand[dataid,])
        
        log_lik = log_lik + sum(log_of_normpdf(data[,dataid], mean_trj_current, sqrt(var_list[,dataid])))
        log_lik_cand = log_lik_cand + sum(log_of_normpdf(data[,dataid], mean_trj_cand, sqrt(var_list[,dataid])))
      }
      acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
      
      if (runif(1) < acceptance_ratio){
        thetaCurrent = thetaCand
        switch(3, 
               A_update[rep,dataid] <- 1, 
               B_update[rep,dataid] <- 1, 
               alpha_update[rep,dataid] <- 1, 
               beta_update[rep,dataid] <- 1)
      }else{
        thetaCand = thetaCurrent
      }
    }
    
    
    
    ## Start to sample the hyper parameters for all populations.
    for (hyp_param_id in 1:(2*num_param)){
      if (hyp_param_id == 5 | hyp_param_id == 6){
        thetaHypCand = thetaHypCurrent
      }else{
        
        hyp_param_cand <- rnorm(1, mean = thetaHypCurrent[hyp_param_id], sd = sqrt(hyper_prop_var[hyp_param_id]))
        
        if (hyp_param_cand > 0){
          thetaHypCand[hyp_param_id] = hyp_param_cand
          
          gamma_pri_alpha = hyper_prior_mean[hyp_param_id]^2 / hyper_prior_var[hyp_param_id]
          gamma_pri_beta = hyper_prior_var[hyp_param_id]/hyper_prior_mean[hyp_param_id]
          
          log_pri = log_of_gampdf(thetaHypCurrent[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
          log_pri_cand = log_of_gampdf(thetaHypCand[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
          
          if (hyp_param_id %% 2 == 1){
            log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id], num_data), rep(thetaHypCurrent[hyp_param_id+1], num_data)))
            log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id], num_data), rep(thetaHypCand[hyp_param_id+1], num_data)))
          }else{
            log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id-1], num_data), rep(thetaHypCurrent[hyp_param_id], num_data)))
            log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id-1], num_data), rep(thetaHypCand[hyp_param_id], num_data)))
          }
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaHypCurrent = thetaHypCand
            update_matrix_hyp[rep, hyp_param_id] = 1
          }else{
            thetaHypCand = thetaHypCurrent
          }
        }
      }
    }
    A_sample[rep,] = thetaCurrent[,1]
    B_sample[rep,] = thetaCurrent[,2]
    alpha_sample[rep,] = thetaCurrent[,3]
    beta_sample[rep,] = thetaCurrent[,4]
    MCMC_sample_hyp[rep, ] = thetaHypCurrent
    
    
  }
  my_list <- list("A_sample" = A_sample[selrow,], "B_sample" = B_sample[selrow,], "alpha_sample" = alpha_sample[selrow,], "beta_sample" = beta_sample[selrow,],
                  "A_update" = A_update, "B_update" = B_update, "alpha_update" = alpha_update, "beta_update" = beta_update,
                  "MCMC_sample_hyp" = MCMC_sample_hyp[selrow,], "update_matrix_hyp" = update_matrix_hyp)
  return(my_list)
}




MCMC_function_ME_lambdaDpri <- function(data, timespan, var_list, hyper_prior_mean = rep(1,8), hyper_prior_var = rep(1000,8), 
                                        theta_init, theta_hyper_init, burn = 0, jump = 1, effnum = 100, 
                                        prop_var = c(50, 0.000001, 0.1, 0.003), hyper_prop_var = c(50,1,0.001,1,1,1,0.01,1),
                                        lambdaDalpha = 10^(-3), lambdaDbeta = 10^3, save.file.name, summary.file.name){
  
  # In this function, data shold be a matrix with more than one column. Each column represents a timeseries data.
  # var_list is also a matrix with the size identical to data.
  
  num_repeat = burn + jump * effnum;
  num_data = dim(data)[2]
  num_param = length(theta_init)
  
  
  MSE_abs_cumul = c()
  MSE_rel_cumul = c()
  MSE_first_rel_cumul = c()
  
  Prod_posterior_mean_cumul = c()
  death_posterior_mean_cumul = c()
  alpha_posterior_mean_cumul = c()    
  beta_posterior_mean_cumul = c()
  timedelay_posterior_mean_cumul = c()
  
  Prod_cv_list_cumul = c()
  death_cv_list_cumul = c()
  alpha_cv_list_cumul = c()
  beta_cv_list_cumul = c()
  Meandelay_cv_list_cumul = c()
  
  birth_acf_lag5_cumul = c()
  death_acf_lag5_cumul = c()
  alpha_acf_lag5_cumul = c()
  beta_acf_lag5_cumul = c()
  meandelay_acf_lag5_cumul = c()
  birth_acf_lag10_cumul = c()
  death_acf_lag10_cumul = c()
  alpha_acf_lag10_cumul = c()
  beta_acf_lag10_cumul = c()
  meandelay_acf_lag10_cumul = c()
  
  MSE_abs_cumul2 = c()
  MSE_rel_cumul2 = c()
  MSE_first_rel_cumul2 = c()
  
  Prod_posterior_mean_cumul2 = c()
  death_posterior_mean_cumul2 = c()
  alpha_posterior_mean_cumul2 = c()    
  beta_posterior_mean_cumul2 = c()
  timedelay_posterior_mean_cumul2 = c()
  
  Prod_cv_list_cumul2 = c()
  death_cv_list_cumul2 = c()
  alpha_cv_list_cumul2 = c()
  beta_cv_list_cumul2 = c()
  Meandelay_cv_list_cumul2 = c()
  
  birth_acf_lag5_cumul2 = c()
  death_acf_lag5_cumul2 = c()
  alpha_acf_lag5_cumul2 = c()
  beta_acf_lag5_cumul2 = c()
  meandelay_acf_lag5_cumul2 = c()
  birth_acf_lag10_cumul2 = c()
  death_acf_lag10_cumul2 = c()
  alpha_acf_lag10_cumul2 = c()
  beta_acf_lag10_cumul2 = c()
  meandelay_acf_lag10_cumul2 = c()
  
  
  A_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  B_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  beta_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  
  A_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  B_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  beta_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  
  MCMC_sample_hyp = matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  update_matrix_hyp =  matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  
  A_sample[1,] = theta_init[1]
  B_sample[1,] = theta_init[2]
  alpha_sample[1,] = theta_init[3]
  beta_sample[1,] = theta_init[4]
  
  MCMC_sample_hyp[1,] = theta_hyper_init
  MCMC_sample_hyp[1,3] = lambdaDalpha
  MCMC_sample_hyp[1,4] = lambdaDbeta
  
  thetaCurrent <- matrix(NA, nrow = num_data, ncol = 4)
  thetaCurrent[,1] <- theta_init[1]
  thetaCurrent[,2] <- theta_init[2]
  thetaCurrent[,3] <- theta_init[3]
  thetaCurrent[,4] <- theta_init[4]
  
  thetaCand <- thetaCurrent
  
  thetaHypCurrent <- theta_hyper_init
  thetaHypCurrent[3] = lambdaDalpha
  thetaHypCurrent[4] = lambdaDbeta
  thetaHypCand <- thetaHypCurrent
  
  for (rep in 2:num_repeat){
    if (rep %% 1000 == 0){
      cat(rep); cat("\n");
      cat(summary.file.name); cat("\n");
    }
    
    ## Start to sample the parameters for each individual
    for (dataid in 1:num_data){
      for (param_id in 1:num_param){
        param_cand <- rnorm(1, mean = thetaCurrent[dataid,param_id], sd = sqrt(prop_var[param_id]))
        if (param_cand > 0){
          thetaCand[dataid, param_id] = param_cand
          mean_trj_current <- mean_trajectory(timespan, thetaCurrent[dataid,])
          mean_trj_cand <- mean_trajectory(timespan, thetaCand[dataid,])
          log_pri = log_of_gampdf(thetaCurrent[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          log_pri_cand = log_of_gampdf(thetaCand[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          
          log_lik = sum(log_of_normpdf(data[,dataid], mean_trj_current, sqrt(var_list[,dataid])))
          log_lik_cand = sum(log_of_normpdf(data[,dataid], mean_trj_cand, sqrt(var_list[,dataid])))
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaCurrent = thetaCand
            switch(param_id, 
                   A_update[rep,dataid] <- 1, 
                   B_update[rep,dataid] <- 1, 
                   alpha_update[rep,dataid] <- 1, 
                   beta_update[rep,dataid] <- 1)
          }else{
            thetaCand = thetaCurrent
          }
        }
      }
    }
    
    ## Start to sample the hyper parameters for all populations.
    for (hyp_param_id in 1:(2*num_param)){
      if (hyp_param_id == 3 | hyp_param_id == 4){
        thetaHypCand = thetaHypCurrent
      }else{
        hyp_param_cand <- rnorm(1, mean = thetaHypCurrent[hyp_param_id], sd = sqrt(hyper_prop_var[hyp_param_id]))
        
        if (hyp_param_cand > 0){
          thetaHypCand[hyp_param_id] = hyp_param_cand
          
          gamma_pri_alpha = hyper_prior_mean[hyp_param_id]^2 / hyper_prior_var[hyp_param_id]
          gamma_pri_beta = hyper_prior_var[hyp_param_id]/hyper_prior_mean[hyp_param_id]
          
          log_pri = log_of_gampdf(thetaHypCurrent[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
          log_pri_cand = log_of_gampdf(thetaHypCand[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
          
          if (hyp_param_id %% 2 == 1){
            log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id], num_data), rep(thetaHypCurrent[hyp_param_id+1], num_data)))
            log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id], num_data), rep(thetaHypCand[hyp_param_id+1], num_data)))
          }else{
            log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id-1], num_data), rep(thetaHypCurrent[hyp_param_id], num_data)))
            log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id-1], num_data), rep(thetaHypCand[hyp_param_id], num_data)))
          }
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaHypCurrent = thetaHypCand
            update_matrix_hyp[rep, hyp_param_id] = 1
          }else{
            thetaHypCand = thetaHypCurrent
          }
        }
      }
    }
    A_sample[rep,] = thetaCurrent[,1]
    B_sample[rep,] = thetaCurrent[,2]
    alpha_sample[rep,] = thetaCurrent[,3]
    beta_sample[rep,] = thetaCurrent[,4]
    MCMC_sample_hyp[rep, ] = thetaHypCurrent
    if (rep %% 55000 == 0 | rep == num_repeat){
      mcmc.result1 <- list("A_sample" = A_sample, "B_sample" = B_sample, "alpha_sample" = alpha_sample, "beta_sample" = beta_sample,
                           "A_update" = A_update, "B_update" = B_update, "alpha_update" = alpha_update, "beta_update" = beta_update,
                           "MCMC_sample_hyp" = MCMC_sample_hyp, "update_matrix_hyp" = update_matrix_hyp, "rep" = rep)
      inputdata1 = data
      # save(mcmc.result1, inputdata1, file = paste(save.file.name, ".Rdata", sep = ""))
    }
    
  }
  return(mcmc.result1)
}

MCMC_function_ME_lambdaDpri_merge_share_alpha <- function(data, timespan, var_list, hyper_prior_mean = rep(1,8), hyper_prior_var = rep(1000,8), 
                                                          theta_init, theta_hyper_init, burn = 0, jump = 1, effnum = 100, 
                                                          prop_var = c(50, 0.000001, 0.1, 0.003), hyper_prop_var = c(50,1,0.001,1,1,1,0.01,1),
                                                          lambdaDalpha = 10^(-3), lambdaDbeta = 10^3, save.file.name, summary.file.name){
  
  # In this function, data shold be a matrix with more than one column. Each column represents a timeseries data.
  # var_list is also a matrix with the size identical to data.
  
  num_repeat = burn + jump * effnum;
  num_data = dim(data)[2]
  num_param = length(theta_init)
  
  
  MSE_abs_cumul = c()
  MSE_rel_cumul = c()
  MSE_first_rel_cumul = c()
  
  Prod_posterior_mean_cumul = c()
  death_posterior_mean_cumul = c()
  alpha_posterior_mean_cumul = c()    
  beta_posterior_mean_cumul = c()
  timedelay_posterior_mean_cumul = c()
  
  Prod_cv_list_cumul = c()
  death_cv_list_cumul = c()
  alpha_cv_list_cumul = c()
  beta_cv_list_cumul = c()
  Meandelay_cv_list_cumul = c()
  
  birth_acf_lag5_cumul = c()
  death_acf_lag5_cumul = c()
  alpha_acf_lag5_cumul = c()
  beta_acf_lag5_cumul = c()
  meandelay_acf_lag5_cumul = c()
  birth_acf_lag10_cumul = c()
  death_acf_lag10_cumul = c()
  alpha_acf_lag10_cumul = c()
  beta_acf_lag10_cumul = c()
  meandelay_acf_lag10_cumul = c()
  
  MSE_abs_cumul2 = c()
  MSE_rel_cumul2 = c()
  MSE_first_rel_cumul2 = c()
  
  Prod_posterior_mean_cumul2 = c()
  death_posterior_mean_cumul2 = c()
  alpha_posterior_mean_cumul2 = c()    
  beta_posterior_mean_cumul2 = c()
  timedelay_posterior_mean_cumul2 = c()
  
  Prod_cv_list_cumul2 = c()
  death_cv_list_cumul2 = c()
  alpha_cv_list_cumul2 = c()
  beta_cv_list_cumul2 = c()
  Meandelay_cv_list_cumul2 = c()
  
  birth_acf_lag5_cumul2 = c()
  death_acf_lag5_cumul2 = c()
  alpha_acf_lag5_cumul2 = c()
  beta_acf_lag5_cumul2 = c()
  meandelay_acf_lag5_cumul2 = c()
  birth_acf_lag10_cumul2 = c()
  death_acf_lag10_cumul2 = c()
  alpha_acf_lag10_cumul2 = c()
  beta_acf_lag10_cumul2 = c()
  meandelay_acf_lag10_cumul2 = c()
  
  
  A_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  B_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  beta_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  
  A_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  B_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  beta_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  
  MCMC_sample_hyp = matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  update_matrix_hyp =  matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  
  A_sample[1,] = theta_init[1]
  B_sample[1,] = theta_init[2]
  alpha_sample[1,] = theta_init[3]
  beta_sample[1,] = theta_init[4]
  
  MCMC_sample_hyp[1,] = theta_hyper_init
  MCMC_sample_hyp[1,3] = lambdaDalpha
  MCMC_sample_hyp[1,4] = lambdaDbeta
  
  thetaCurrent <- matrix(NA, nrow = num_data, ncol = 4)
  thetaCurrent[,1] <- theta_init[1]
  thetaCurrent[,2] <- theta_init[2]
  thetaCurrent[,3] <- theta_init[3]
  thetaCurrent[,4] <- theta_init[4]
  
  thetaCand <- thetaCurrent
  
  thetaHypCurrent <- theta_hyper_init
  thetaHypCurrent[3] = lambdaDalpha
  thetaHypCurrent[4] = lambdaDbeta
  thetaHypCand <- thetaHypCurrent
  for (rep in 2:num_repeat){
    # if (rep %% 100 == 0){
    #   cat(rep); cat("\n");
    # }
    
    ## Start to sample the parameters for each individual
    for (dataid in 1:num_data){
      for (param_id in 1:num_param){
        param_cand <- rnorm(1, mean = thetaCurrent[dataid,param_id], sd = sqrt(prop_var[param_id]))
        if (param_cand > 0){
          if(param_id == 3){
            thetaCand[, param_id] = param_cand
            
          }else{
            thetaCand[dataid, param_id] = param_cand
          }  
          mean_trj_current <- mean_trajectory(timespan[!is.na(timespan[,dataid]),dataid], thetaCurrent[dataid,])
          mean_trj_cand <- mean_trajectory(timespan[!is.na(timespan[,dataid]),dataid], thetaCand[dataid,])
          log_pri = log_of_gampdf(thetaCurrent[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          log_pri_cand = log_of_gampdf(thetaCand[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          
          log_lik = sum(log_of_normpdf(data[!is.na(data[,dataid]),dataid], mean_trj_current, sqrt(var_list[!is.na(var_list[,dataid]),dataid])))
          log_lik_cand = sum(log_of_normpdf(data[!is.na(data[,dataid]),dataid], mean_trj_cand, sqrt(var_list[!is.na(var_list[,dataid]),dataid])))
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaCurrent = thetaCand
            switch(param_id, 
                   A_update[rep,dataid] <- 1, 
                   B_update[rep,dataid] <- 1, 
                   alpha_update[rep,dataid] <- 1, 
                   beta_update[rep,dataid] <- 1)
          }else{
            thetaCand = thetaCurrent
          }
        }
      }
    }
    
    ## Start to sample the hyper parameters for all populations.
    for (hyp_param_id in 1:(2*num_param)){
      if (hyp_param_id == 3 | hyp_param_id == 4 |hyp_param_id == 5 | hyp_param_id == 6){
        thetaHypCand = thetaHypCurrent
      }else{
        hyp_param_cand <- rnorm(1, mean = thetaHypCurrent[hyp_param_id], sd = sqrt(hyper_prop_var[hyp_param_id]))
        
        if (hyp_param_cand > 0){
          thetaHypCand[hyp_param_id] = hyp_param_cand
          
          gamma_pri_alpha = hyper_prior_mean[hyp_param_id]^2 / hyper_prior_var[hyp_param_id]
          gamma_pri_beta = hyper_prior_var[hyp_param_id]/hyper_prior_mean[hyp_param_id]
          
          log_pri = log_of_gampdf(thetaHypCurrent[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
          log_pri_cand = log_of_gampdf(thetaHypCand[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
          
          if (hyp_param_id %% 2 == 1){
            log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id], num_data), rep(thetaHypCurrent[hyp_param_id+1], num_data)))
            log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id], num_data), rep(thetaHypCand[hyp_param_id+1], num_data)))
          }else{
            log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id-1], num_data), rep(thetaHypCurrent[hyp_param_id], num_data)))
            log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id-1], num_data), rep(thetaHypCand[hyp_param_id], num_data)))
          }
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaHypCurrent = thetaHypCand
            update_matrix_hyp[rep, hyp_param_id] = 1
          }else{
            thetaHypCand = thetaHypCurrent
          }
        }
      }
    }
    A_sample[rep,] = thetaCurrent[,1]
    B_sample[rep,] = thetaCurrent[,2]
    alpha_sample[rep,] = thetaCurrent[,3]
    beta_sample[rep,] = thetaCurrent[,4]
    MCMC_sample_hyp[rep, ] = thetaHypCurrent
    
    
    
    if (rep %% 55000 == 0 | rep == num_repeat){
      mcmc.result1 <- list("A_sample" = A_sample, "B_sample" = B_sample, "alpha_sample" = alpha_sample, "beta_sample" = beta_sample,
                           "A_update" = A_update, "B_update" = B_update, "alpha_update" = alpha_update, "beta_update" = beta_update,
                           "MCMC_sample_hyp" = MCMC_sample_hyp, "update_matrix_hyp" = update_matrix_hyp, "rep" = rep)
      inputdata1 = data
      # save(mcmc.result1, inputdata1, file = paste(save.file.name, ".Rdata", sep = ""))
    }
  }
  return(mcmc.result1)
}


MCMC_function_ME_lambdaDpri_merge <- function(data, timespan, var_list, hyper_prior_mean = rep(1,8), hyper_prior_var = rep(1000,8), 
                                              theta_init, theta_hyper_init, burn = 0, jump = 1, effnum = 100, 
                                              prop_var = c(50, 0.000001, 0.1, 0.003), hyper_prop_var = c(50,1,0.001,1,1,1,0.01,1),
                                              lambdaDalpha = 10^(-3), lambdaDbeta = 10^3, save.file.name, summary.file.name){
  
  # In this function, data shold be a matrix with more than one column. Each column represents a timeseries data.
  # var_list is also a matrix with the size identical to data.
  
  num_repeat = burn + jump * effnum;
  num_data = dim(data)[2]
  num_param = length(theta_init)
  
  
  MSE_abs_cumul = c()
  MSE_rel_cumul = c()
  MSE_first_rel_cumul = c()
  
  Prod_posterior_mean_cumul = c()
  death_posterior_mean_cumul = c()
  alpha_posterior_mean_cumul = c()    
  beta_posterior_mean_cumul = c()
  timedelay_posterior_mean_cumul = c()
  
  Prod_cv_list_cumul = c()
  death_cv_list_cumul = c()
  alpha_cv_list_cumul = c()
  beta_cv_list_cumul = c()
  Meandelay_cv_list_cumul = c()
  
  birth_acf_lag5_cumul = c()
  death_acf_lag5_cumul = c()
  alpha_acf_lag5_cumul = c()
  beta_acf_lag5_cumul = c()
  meandelay_acf_lag5_cumul = c()
  birth_acf_lag10_cumul = c()
  death_acf_lag10_cumul = c()
  alpha_acf_lag10_cumul = c()
  beta_acf_lag10_cumul = c()
  meandelay_acf_lag10_cumul = c()
  
  MSE_abs_cumul2 = c()
  MSE_rel_cumul2 = c()
  MSE_first_rel_cumul2 = c()
  
  Prod_posterior_mean_cumul2 = c()
  death_posterior_mean_cumul2 = c()
  alpha_posterior_mean_cumul2 = c()    
  beta_posterior_mean_cumul2 = c()
  timedelay_posterior_mean_cumul2 = c()
  
  Prod_cv_list_cumul2 = c()
  death_cv_list_cumul2 = c()
  alpha_cv_list_cumul2 = c()
  beta_cv_list_cumul2 = c()
  Meandelay_cv_list_cumul2 = c()
  
  birth_acf_lag5_cumul2 = c()
  death_acf_lag5_cumul2 = c()
  alpha_acf_lag5_cumul2 = c()
  beta_acf_lag5_cumul2 = c()
  meandelay_acf_lag5_cumul2 = c()
  birth_acf_lag10_cumul2 = c()
  death_acf_lag10_cumul2 = c()
  alpha_acf_lag10_cumul2 = c()
  beta_acf_lag10_cumul2 = c()
  meandelay_acf_lag10_cumul2 = c()
  
  
  A_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  B_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  beta_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  
  A_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  B_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  beta_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  
  MCMC_sample_hyp = matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  update_matrix_hyp =  matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  
  A_sample[1,] = theta_init[1]
  B_sample[1,] = theta_init[2]
  alpha_sample[1,] = theta_init[3]
  beta_sample[1,] = theta_init[4]
  
  MCMC_sample_hyp[1,] = theta_hyper_init
  MCMC_sample_hyp[1,3] = lambdaDalpha
  MCMC_sample_hyp[1,4] = lambdaDbeta
  
  thetaCurrent <- matrix(NA, nrow = num_data, ncol = 4)
  thetaCurrent[,1] <- theta_init[1]
  thetaCurrent[,2] <- theta_init[2]
  thetaCurrent[,3] <- theta_init[3]
  thetaCurrent[,4] <- theta_init[4]
  
  thetaCand <- thetaCurrent
  
  thetaHypCurrent <- theta_hyper_init
  thetaHypCurrent[3] = lambdaDalpha
  thetaHypCurrent[4] = lambdaDbeta
  thetaHypCand <- thetaHypCurrent
  for (rep in 2:num_repeat){
    # if (rep %% 100 == 0){
    #   cat(rep); cat("\n");
    # }
    
    ## Start to sample the parameters for each individual
    for (dataid in 1:num_data){
      for (param_id in 1:num_param){
        param_cand <- rnorm(1, mean = thetaCurrent[dataid,param_id], sd = sqrt(prop_var[param_id]))
        if (param_cand > 0){
          thetaCand[dataid, param_id] = param_cand
          mean_trj_current <- mean_trajectory(timespan[!is.na(timespan[,dataid]),dataid], thetaCurrent[dataid,])
          mean_trj_cand <- mean_trajectory(timespan[!is.na(timespan[,dataid]),dataid], thetaCand[dataid,])
          log_pri = log_of_gampdf(thetaCurrent[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          log_pri_cand = log_of_gampdf(thetaCand[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          
          log_lik = sum(log_of_normpdf(data[!is.na(data[,dataid]),dataid], mean_trj_current, sqrt(var_list[!is.na(var_list[,dataid]),dataid])))
          log_lik_cand = sum(log_of_normpdf(data[!is.na(data[,dataid]),dataid], mean_trj_cand, sqrt(var_list[!is.na(var_list[,dataid]),dataid])))
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaCurrent = thetaCand
            switch(param_id, 
                   A_update[rep,dataid] <- 1, 
                   B_update[rep,dataid] <- 1, 
                   alpha_update[rep,dataid] <- 1, 
                   beta_update[rep,dataid] <- 1)
          }else{
            thetaCand = thetaCurrent
          }
        }
      }
    }
    
    ## Start to sample the hyper parameters for all populations.
    for (hyp_param_id in 1:(2*num_param)){
      if (hyp_param_id == 3 | hyp_param_id == 4){
        thetaHypCand = thetaHypCurrent
      }else{
        hyp_param_cand <- rnorm(1, mean = thetaHypCurrent[hyp_param_id], sd = sqrt(hyper_prop_var[hyp_param_id]))
        
        if (hyp_param_cand > 0){
          thetaHypCand[hyp_param_id] = hyp_param_cand
          
          gamma_pri_alpha = hyper_prior_mean[hyp_param_id]^2 / hyper_prior_var[hyp_param_id]
          gamma_pri_beta = hyper_prior_var[hyp_param_id]/hyper_prior_mean[hyp_param_id]
          
          log_pri = log_of_gampdf(thetaHypCurrent[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
          log_pri_cand = log_of_gampdf(thetaHypCand[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
          
          if (hyp_param_id %% 2 == 1){
            log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id], num_data), rep(thetaHypCurrent[hyp_param_id+1], num_data)))
            log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id], num_data), rep(thetaHypCand[hyp_param_id+1], num_data)))
          }else{
            log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id-1], num_data), rep(thetaHypCurrent[hyp_param_id], num_data)))
            log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id-1], num_data), rep(thetaHypCand[hyp_param_id], num_data)))
          }
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaHypCurrent = thetaHypCand
            update_matrix_hyp[rep, hyp_param_id] = 1
          }else{
            thetaHypCand = thetaHypCurrent
          }
        }
      }
    }
    A_sample[rep,] = thetaCurrent[,1]
    B_sample[rep,] = thetaCurrent[,2]
    alpha_sample[rep,] = thetaCurrent[,3]
    beta_sample[rep,] = thetaCurrent[,4]
    MCMC_sample_hyp[rep, ] = thetaHypCurrent
    if (rep %% 55000 == 0 | rep == num_repeat){
      mcmc.result1 <- list("A_sample" = A_sample, "B_sample" = B_sample, "alpha_sample" = alpha_sample, "beta_sample" = beta_sample,
                           "A_update" = A_update, "B_update" = B_update, "alpha_update" = alpha_update, "beta_update" = beta_update,
                           "MCMC_sample_hyp" = MCMC_sample_hyp, "update_matrix_hyp" = update_matrix_hyp, "rep" = rep)
      inputdata1 = data
      # save(mcmc.result1, inputdata1, file = paste(save.file.name, ".Rdata", sep = ""))
    }
    
  }
  return(mcmc.result1)
}

MCMC_function_ME_lambdaDpri_share_alpha <- function(data, timespan, var_list, hyper_prior_mean = rep(1,8), hyper_prior_var = rep(1000,8), 
                                                    theta_init, theta_hyper_init, burn = 0, jump = 1, effnum = 100, 
                                                    prop_var = c(50, 0.000001, 0.1, 0.003), hyper_prop_var = c(50,1,0.001,1,1,1,0.01,1),
                                                    lambdaDalpha = 10^(-3), lambdaDbeta = 10^3, save.file.name, summary.file.name){
  
  # In this function, data should be a matrix with more than one column. Each column represents a timeseries data.
  # var_list is also a matrix with the size identical to data.
  
  num_repeat = burn + jump * effnum;
  num_data = dim(data)[2]
  num_param = length(theta_init)
  
  
  MSE_abs_cumul = c()
  MSE_rel_cumul = c()
  MSE_first_rel_cumul = c()
  
  Prod_posterior_mean_cumul = c()
  death_posterior_mean_cumul = c()
  alpha_posterior_mean_cumul = c()    
  beta_posterior_mean_cumul = c()
  timedelay_posterior_mean_cumul = c()
  
  Prod_cv_list_cumul = c()
  death_cv_list_cumul = c()
  alpha_cv_list_cumul = c()
  beta_cv_list_cumul = c()
  Meandelay_cv_list_cumul = c()
  
  birth_acf_lag5_cumul = c()
  death_acf_lag5_cumul = c()
  alpha_acf_lag5_cumul = c()
  beta_acf_lag5_cumul = c()
  meandelay_acf_lag5_cumul = c()
  birth_acf_lag10_cumul = c()
  death_acf_lag10_cumul = c()
  alpha_acf_lag10_cumul = c()
  beta_acf_lag10_cumul = c()
  meandelay_acf_lag10_cumul = c()
  
  MSE_abs_cumul2 = c()
  MSE_rel_cumul2 = c()
  MSE_first_rel_cumul2 = c()
  
  Prod_posterior_mean_cumul2 = c()
  death_posterior_mean_cumul2 = c()
  alpha_posterior_mean_cumul2 = c()    
  beta_posterior_mean_cumul2 = c()
  timedelay_posterior_mean_cumul2 = c()
  
  Prod_cv_list_cumul2 = c()
  death_cv_list_cumul2 = c()
  alpha_cv_list_cumul2 = c()
  beta_cv_list_cumul2 = c()
  Meandelay_cv_list_cumul2 = c()
  
  birth_acf_lag5_cumul2 = c()
  death_acf_lag5_cumul2 = c()
  alpha_acf_lag5_cumul2 = c()
  beta_acf_lag5_cumul2 = c()
  meandelay_acf_lag5_cumul2 = c()
  birth_acf_lag10_cumul2 = c()
  death_acf_lag10_cumul2 = c()
  alpha_acf_lag10_cumul2 = c()
  beta_acf_lag10_cumul2 = c()
  meandelay_acf_lag10_cumul2 = c()
  
  
  
  A_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  B_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  beta_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  
  A_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  B_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  beta_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  
  MCMC_sample_hyp = matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  update_matrix_hyp =  matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  
  A_sample[1,] = theta_init[1]
  B_sample[1,] = theta_init[2]
  alpha_sample[1,] = theta_init[3]
  beta_sample[1,] = theta_init[4]
  
  MCMC_sample_hyp[1,] = theta_hyper_init
  MCMC_sample_hyp[1,3] = lambdaDalpha
  MCMC_sample_hyp[1,4] = lambdaDbeta
  
  thetaCurrent <- matrix(NA, nrow = num_data, ncol = 4)
  thetaCurrent[,1] <- theta_init[1]
  thetaCurrent[,2] <- theta_init[2]
  thetaCurrent[,3] <- theta_init[3]
  thetaCurrent[,4] <- theta_init[4]
  
  thetaCand <- thetaCurrent
  
  thetaHypCurrent <- theta_hyper_init
  thetaHypCurrent[3] = lambdaDalpha
  thetaHypCurrent[4] = lambdaDbeta
  thetaHypCand <- thetaHypCurrent
  for (rep in 2:num_repeat){
    if (rep %% 100 == 0 & rep < 2000){
      cat(rep); cat("\n");
    }
    
    ## Start to sample the parameters for each individual
    for (dataid in 1:num_data){
      for (param_id in 1:num_param){
        param_cand <- rnorm(1, mean = thetaCurrent[dataid,param_id], sd = sqrt(prop_var[param_id]))
        if (param_cand > 0){
          if(param_id == 3){
            # thetaCand[, param_id] = param_cand
            next
          }else{
            thetaCand[dataid, param_id] = param_cand
          }  
          mean_trj_current <- mean_trajectory(timespan, thetaCurrent[dataid,])
          mean_trj_cand <- mean_trajectory(timespan, thetaCand[dataid,])
          log_pri = log_of_gampdf(thetaCurrent[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          log_pri_cand = log_of_gampdf(thetaCand[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          
          log_lik = sum(log_of_normpdf(data[,dataid], mean_trj_current, sqrt(var_list[,dataid])))
          log_lik_cand = sum(log_of_normpdf(data[,dataid], mean_trj_cand, sqrt(var_list[,dataid])))
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaCurrent = thetaCand
            switch(param_id, 
                   A_update[rep,dataid] <- 1, 
                   B_update[rep,dataid] <- 1, 
                   alpha_update[rep,dataid] <- 1, 
                   beta_update[rep,dataid] <- 1)
          }else{
            thetaCand = thetaCurrent
          }
        }
      }
    }
    ## Sample the shared alpha values for all individual
    
    param_cand <- rnorm(1, mean = thetaCurrent[1,3], sd = sqrt(prop_var[3]))
    
    if (param_cand > 0){
      thetaCand[, 3] = param_cand
      log_pri = log_of_gampdf(thetaCurrent[1,3], thetaHypCurrent[5], thetaHypCurrent[6])
      log_pri_cand = log_of_gampdf(thetaCand[1,3], thetaHypCurrent[5], thetaHypCurrent[6])
      log_lik = 0
      log_lik_cand = 0
      for (dataid in 1:num_data){
        mean_trj_current <- mean_trajectory(timespan, thetaCurrent[dataid,])
        mean_trj_cand <- mean_trajectory(timespan, thetaCand[dataid,])
        
        log_lik = log_lik + sum(log_of_normpdf(data[,dataid], mean_trj_current, sqrt(var_list[,dataid])))
        log_lik_cand = log_lik_cand + sum(log_of_normpdf(data[,dataid], mean_trj_cand, sqrt(var_list[,dataid])))
      }
      acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
      
      if (runif(1) < acceptance_ratio){
        thetaCurrent = thetaCand
        switch(3, 
               A_update[rep,dataid] <- 1, 
               B_update[rep,dataid] <- 1, 
               alpha_update[rep,dataid] <- 1, 
               beta_update[rep,dataid] <- 1)
      }else{
        thetaCand = thetaCurrent
      }
    }
    
    ## Start to sample the hyper parameters for all populations.
    for (hyp_param_id in 1:(2*num_param)){
      if (hyp_param_id == 3 | hyp_param_id == 4 |hyp_param_id == 5 | hyp_param_id == 6){
        thetaHypCand = thetaHypCurrent
      }else{
        hyp_param_cand <- rnorm(1, mean = thetaHypCurrent[hyp_param_id], sd = sqrt(hyper_prop_var[hyp_param_id]))
        
        if (hyp_param_cand > 0){
          thetaHypCand[hyp_param_id] = hyp_param_cand
          
          gamma_pri_alpha = hyper_prior_mean[hyp_param_id]^2 / hyper_prior_var[hyp_param_id]
          gamma_pri_beta = hyper_prior_var[hyp_param_id]/hyper_prior_mean[hyp_param_id]
          
          log_pri = log_of_gampdf(thetaHypCurrent[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
          log_pri_cand = log_of_gampdf(thetaHypCand[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
          
          if (hyp_param_id %% 2 == 1){
            log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id], num_data), rep(thetaHypCurrent[hyp_param_id+1], num_data)))
            log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id], num_data), rep(thetaHypCand[hyp_param_id+1], num_data)))
          }else{
            log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id-1], num_data), rep(thetaHypCurrent[hyp_param_id], num_data)))
            log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id-1], num_data), rep(thetaHypCand[hyp_param_id], num_data)))
          }
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaHypCurrent = thetaHypCand
            update_matrix_hyp[rep, hyp_param_id] = 1
          }else{
            thetaHypCand = thetaHypCurrent
          }
        }
      }
    }
    A_sample[rep,] = thetaCurrent[,1]
    B_sample[rep,] = thetaCurrent[,2]
    alpha_sample[rep,] = thetaCurrent[,3]
    beta_sample[rep,] = thetaCurrent[,4]
    MCMC_sample_hyp[rep, ] = thetaHypCurrent
    
    if (rep %% 55000 == 0 | rep == num_repeat){
      mcmc.result1 <- list("A_sample" = A_sample, "B_sample" = B_sample, "alpha_sample" = alpha_sample, "beta_sample" = beta_sample,
                           "A_update" = A_update, "B_update" = B_update, "alpha_update" = alpha_update, "beta_update" = beta_update,
                           "MCMC_sample_hyp" = MCMC_sample_hyp, "update_matrix_hyp" = update_matrix_hyp, "rep" = rep)
      inputdata1 = data
      # save(mcmc.result1, inputdata1, file = paste(save.file.name, ".Rdata", sep = ""))
    }
  }
  return(mcmc.result1)
}



TimeDelayGillespieforX <- function(A.X, B.X, alpha.X, beta.X, repnum = 3000, maxT = 100, Volume = 1, time.interval = 1){
  X <- 0
  XList <- rep(NA, repnum)
  Xbirth <- rep(0, maxT/time.interval)
  Xdeath <- rep(0, maxT/time.interval)
  currentTime <- 0
  TList <- rep(NA, repnum)
  n <- 1
  k <- 1
  stackTimeX <- c()
  
  for (i in 1:repnum){
    a1 <- Volume * A.X
    a2 <- B.X * X
    a0 <- sum(a1,a2)
    # r2 <- runif(1)
    currentTime <- currentTime + rexp(1, rate = a0)
    
    stackTimeX <- sort(stackTimeX)  
    if(!(is.null(stackTimeX))){
      minStack <- min(stackTimeX)
    } else {
      minStack <- Inf
    }
    if (currentTime < minStack){                                    
      r1 <- runif(1)
      if (r1 < a1/a0){
        XList[i] <- X
        TList[i] <- currentTime
        stackTimeX <- c(stackTimeX, currentTime + k*rgamma(n=1, shape = alpha.X, scale = beta.X))
        #stackTimeX <- c(stackTimeX, currentTime)
      } else{
        # death of X occurs
        X <- X-1;
        XList[i] <- X
        TList[i] <- currentTime
        Xdeath[ceiling(currentTime/time.interval)] = Xdeath[ceiling(currentTime/time.interval)] + 1
      } 
    } else{
      X <- X+1;
      XList[i] <- X
      TList[i] <- minStack
      currentTime <- minStack
      stackTimeX <- stackTimeX[-1]
      Xbirth[ceiling(currentTime/time.interval)] = Xbirth[ceiling(currentTime/time.interval)] + 1
    }
    if (currentTime > maxT){
      break
    }
  }
  
  XList <- XList/Volume
  Xbirth <- Xbirth/Volume
  Xdeath <- Xdeath/Volume
  
  my_list <- list("XList" = XList, "TList" = TList, "Xbirth" = Xbirth, "Xdeath" = Xdeath)
  return(my_list)
}



TimeDelayGillespieforX_paramChange <- function(A.X, B.X, alpha.X, beta.X, repnum = 3000, maxT = 100, Volume = 1, changeTime = 3, A.X.new, B.X.new, alpha.X.new, beta.X.new){
  X <- 0
  XList <- rep(NA, repnum)
  Xbirth <- rep(0, maxT)
  Xdeath <- rep(0, maxT)
  currentTime <- 0
  TList <- rep(NA, repnum)
  n <- 1
  k <- 1
  stackTimeX <- c()
  
  for (i in 1:repnum){
    if (currentTime > changeTime){
      A.X = A.X.new
      B.X = B.X.new
      alpha.X = alpha.X.new
      beta.X = beta.X.new
    }
    a1 <- Volume * A.X
    a2 <- B.X * X
    a0 <- sum(a1,a2)
    # r2 <- runif(1)
    currentTime <- currentTime + rexp(1, rate = a0)
    
    stackTimeX <- sort(stackTimeX)  
    if(!(is.null(stackTimeX))){
      minStack <- min(stackTimeX)
    } else {
      minStack <- Inf
    }
    if (currentTime < minStack){                                    
      r1 <- runif(1)
      if (r1 < a1/a0){
        XList[i] <- X
        TList[i] <- currentTime
        stackTimeX <- c(stackTimeX, currentTime + k*rgamma(n=1, shape = alpha.X, scale = beta.X))
        #stackTimeX <- c(stackTimeX, currentTime)
      } else{
        # death of X occurs
        X <- X-1;
        XList[i] <- X
        TList[i] <- currentTime
        Xdeath[ceiling(currentTime)] = Xdeath[ceiling(currentTime)] + 1
      } 
    } else{
      X <- X+1;
      XList[i] <- X
      TList[i] <- minStack
      currentTime <- minStack
      stackTimeX <- stackTimeX[-1]
      Xbirth[ceiling(currentTime)] = Xbirth[ceiling(currentTime)] + 1
    }
    if (currentTime > maxT){
      break
    }
  }
  
  XList <- XList/Volume
  Xbirth <- Xbirth/Volume
  Xdeath <- Xdeath/Volume
  
  my_list <- list("XList" = XList, "TList" = TList, "Xbirth" = Xbirth, "Xdeath" = Xdeath)
  return(my_list)
}


TimeDelayGillespieforX_NFL <- function(A.X, B.X, alpha.X, beta.X, A_nfl, repnum = 3000, maxT = 100, Volume = 1, time.interval = 1){
  X <- 0
  XList <- rep(NA, repnum)
  Xbirth <- rep(0, maxT/time.interval)
  Xdeath <- rep(0, maxT/time.interval)
  currentTime <- 0
  TList <- rep(NA, repnum)
  n <- 1
  k <- 1
  stackTimeX <- c()
  
  for (i in 1:repnum){
    a1 <- Volume * A.X * (max(0, 1-X/A_nfl))
    a2 <- B.X * X
    a0 <- sum(a1,a2)
    
    # r2 <- runif(1)
    if(a0 == 0 & is.null(stackTimeX) & B.X == 0){
      break
    }else if(a0 == 0 & B.X == 0){
      minStack <- min(stackTimeX)
      X <- X+1;
      XList[i] <- X
      TList[i] <- minStack
      currentTime <- minStack
      stackTimeX <- stackTimeX[-1]
      Xbirth[ceiling(currentTime/time.interval)] = Xbirth[ceiling(currentTime/time.interval)] + 1
      if (currentTime > maxT){
        break
      }
    }else{
      #   
      currentTime <- currentTime + rexp(1, rate = a0)
      
      stackTimeX <- sort(stackTimeX)  
      if(!(is.null(stackTimeX))){
        minStack <- min(stackTimeX)
      } else {
        minStack <- Inf
      }
      if (currentTime < minStack){                                    
        r1 <- runif(1)
        if (r1 < a1/a0){
          XList[i] <- X
          TList[i] <- currentTime
          stackTimeX <- c(stackTimeX, currentTime + k*rgamma(n=1, shape = alpha.X, scale = beta.X))
          #stackTimeX <- c(stackTimeX, currentTime)
        } else{
          # death of X occurs
          X <- X-1;
          XList[i] <- X
          TList[i] <- currentTime
          Xdeath[ceiling(currentTime/time.interval)] = Xdeath[ceiling(currentTime/time.interval)] + 1
        } 
      } else{
        X <- X+1;
        XList[i] <- X
        TList[i] <- minStack
        currentTime <- minStack
        stackTimeX <- stackTimeX[-1]
        Xbirth[ceiling(currentTime/time.interval)] = Xbirth[ceiling(currentTime/time.interval)] + 1
      }
      if (currentTime > maxT){
        break
      }
    }
  }
  
  XList <- XList/Volume
  Xbirth <- Xbirth/Volume
  Xdeath <- Xdeath/Volume
  
  my_list <- list("XList" = XList, "TList" = TList, "Xbirth" = Xbirth, "Xdeath" = Xdeath)
  return(my_list)
}



numer_int <- function(xspan, yval){
  # xspan has to be evenly spaced.
  if(length(xspan) == 1){
    intval0 = 0
  }else{
    h = xspan[2] - xspan[1]
    # if the number of data points is even then it uses the trapezoidal rule.
    if(length(xspan)%%2 == 0){
      intval0 = 0
      for(jj in 2:length(xspan)){
        intval0 = intval0 + (xspan[jj] - xspan[jj-1]) * (yval[jj]+yval[jj-1])/2
      }
    }else{ # if the number of data points is odd then it uses Simpson's rule.
      intval0 = 0
      for(jj in seq(from=2,by=2,to=length(xspan))){
        intval0 = intval0 + h * (yval[jj-1] + 4*yval[jj] + yval[jj+1])/3
      }
    }
  }
  return(intval0)
}


mean_trajectory_NFL <- function(timespan, theta){
  lambda_b = theta[1]
  lambda_d = theta[2]
  alpha = theta[3]
  beta = theta[4]
  A = theta[5]
  
  h = timespan[2] - timespan[1]
  xspan = c(0, timespan)
  y = rep(NA, length(xspan))
  hs = rep(NA, length(xspan))
  hs_tilde = rep(NA, length(xspan))
  y[1] = 0
  hs[1] = 0
  hs_tilde[1] = 0
  
  for(ii in 1:(length(xspan)-1)){
    hs[ii] = numer_int(xspan[1:ii], lambda_b * pmax(0, 1-1/A * y[1:ii]) * dgamma(max(xspan[1:ii]) - xspan[1:ii], shape = alpha, scale = beta))
    f0 = hs[ii] + numer_int(xspan[1:ii], hs[1:ii] * (-lambda_d) * exp(-lambda_d*(max(xspan[1:ii]) - xspan[1:ii])))
    y_tilde = y[ii] + h * f0
    y_tilde_ext = c(y[1:ii], y_tilde)
    hs_tilde[ii+1] = numer_int(xspan[1:(ii+1)],lambda_b * pmax(0, 1-1/A * y_tilde_ext) * dgamma(max(xspan[1:(ii+1)]) - xspan[1:(ii+1)], shape = alpha, scale = beta))
    # f0_tilde = hs_tilde[ii+1] + numer_int(xspan[1:(ii+1)], hs_tilde[1:(ii+1)] * (-lambda_d) * exp(-lambda_d*(max(xspan[1:(ii+1)]) - xspan[1:(ii+1)])))
    f0_tilde = hs_tilde[ii+1] + numer_int(xspan[1:(ii+1)], c(hs[1:ii], hs_tilde[ii+1]) * (-lambda_d) * exp(-lambda_d*(max(xspan[1:(ii+1)]) - xspan[1:(ii+1)])))
    y[ii+1] = y[ii] + h/2 * (f0 + f0_tilde) 
  }
  return(y[-1])
}

mean_trajectory_NFL_delay_feed <- function(timespan, theta){
  lambda_b = theta[1]
  lambda_d = theta[2]
  alpha = theta[3]
  beta = theta[4]
  A = theta[5]
  tau2 = theta[6]
  
  h = timespan[2] - timespan[1]
  xspan = c(0, timespan)
  
  delay_id = sum(xspan <= tau2) # if xspan[1] !=0 then it might be adjusted.
  p0 = (tau2 - xspan[delay_id])/h
  
  y = rep(NA, length(xspan))
  hs = rep(NA, length(xspan))
  hs_tilde = rep(NA, length(xspan))
  y_delay = rep(NA, length(xspan))
  y[1] = 0
  y_delay[1] = 0
  
  hs[1] = 0
  hs_tilde[1] = 0
  
  
  for(ii in 1:(length(xspan)-1)){
    hs[ii] = numer_int(xspan[1:ii], lambda_b * pmax(0, 1-1/A * y_delay[1:ii]) * dgamma(max(xspan[1:ii]) - xspan[1:ii], shape = alpha, scale = beta))
    f0 = hs[ii] + numer_int(xspan[1:ii], hs[1:ii] * (-lambda_d) * exp(-lambda_d*(max(xspan[1:ii]) - xspan[1:ii])))
    y_tilde = y[ii] + h * f0
    if(ii+1 <= delay_id){
      y_delay_tilde = 0  
    }else{
      if(is.na(y[ii+2-delay_id])){
        y_delay_tilde = p0 * y[ii+1-delay_id] + (1-p0)*y_tilde  
      }else{
        y_delay_tilde = p0 * y[ii+1-delay_id] + (1-p0)*y[ii+2-delay_id]
      }
    }
    
    y_delay_tilde_ext = c(y_delay[1:ii], y_delay_tilde)
    hs_tilde[ii+1] = numer_int(xspan[1:(ii+1)],lambda_b * pmax(0, 1-1/A * y_delay_tilde_ext) * dgamma(max(xspan[1:(ii+1)]) - xspan[1:(ii+1)], shape = alpha, scale = beta))
    # f0_tilde = hs_tilde[ii+1] + numer_int(xspan[1:(ii+1)], hs_tilde[1:(ii+1)] * (-lambda_d) * exp(-lambda_d*(max(xspan[1:(ii+1)]) - xspan[1:(ii+1)])))
    f0_tilde = hs_tilde[ii+1] + numer_int(xspan[1:(ii+1)], c(hs[1:ii], hs_tilde[ii+1]) * (-lambda_d) * exp(-lambda_d*(max(xspan[1:(ii+1)]) - xspan[1:(ii+1)])))
    y[ii+1] = y[ii] + h/2 * (f0 + f0_tilde)
    if(ii+1 <= delay_id){
      y_delay[ii+1] = 0  
    }else{
      y_delay[ii+1] = p0 * y[ii+1-delay_id] + (1-p0)*y[ii+2-delay_id]
    }
  }
  return(y[-1])
}



moving_avergae_column <- function(inputdata, windowsize = 5){
  # windowsize must be an odd number.
  # inputdata must be a column vector.
  
  data_length = length(inputdata)
  centeridx = (windowsize+1)/2
  data_shadow_matrix = matrix(NA, nrow = data_length, ncol = windowsize)
  wing_length = (windowsize-1)/2
  
  for (jj in 1:wing_length){
    data_shadow_matrix[(1+jj):data_length, jj] = inputdata[1:(data_length-jj)]
  }
  for (jj in 1:wing_length){
    data_shadow_matrix[1:(data_length-jj), wing_length+jj] = inputdata[(1+jj):data_length]
  } 
  data_shadow_matrix[, windowsize] = inputdata
  outputdata = rowMeans(data_shadow_matrix, na.rm = TRUE)
  return(outputdata)
}



MCMC_function_NFL_delay_feed_fine <- function(data, timespan, var_list, prior_mean, prior_var, theta_init, burn, jump, effnum, prop_var = c(1,0.001, 0.1, 0.003, 1, 1),
                                              save.file.name, summary.file.name){
  # In this function, data should be a matrix with more than one column. Each column represents a timeseries data.
  # var_list is also a matrix with the size identical to data.
  
  total_num = burn + jump * effnum;
  MCMC_sample = matrix(0, nrow = total_num, ncol = length(theta_init))
  update_matrix =  matrix(0, nrow = total_num, ncol = length(theta_init))
  
  MCMC_sample[1,] = theta_init
  
  for (rep in 2:total_num){
    if (rep %% 1000 == 0){
      cat(rep); cat("\n");
      cat(summary.file.name); cat("\n");
    }
    # 
    
    prop_mean1 = MCMC_sample[rep-1,1]
    prop_var1 = prop_var[1] # tune this parameter
    
    theta1_cand = rnorm(1, mean = prop_mean1, sd = sqrt(prop_var1))
    
    if (theta1_cand > 0){
      gamma_pri_alpha_1 = prior_mean[1]^2 / prior_var[1]
      gamma_pri_beta_1 = prior_var[1] / prior_mean[1]
      
      log_pri = log_of_gampdf(MCMC_sample[rep-1, 1], gamma_pri_alpha_1, gamma_pri_beta_1)
      log_pri_st = log_of_gampdf(theta1_cand, gamma_pri_alpha_1, gamma_pri_beta_1)
      
      alpha_tmp = MCMC_sample[rep-1, 3]
      beta_tmp = MCMC_sample[rep-1, 4]
      h = timespan[2] - timespan[1]
      if(alpha_tmp == 1){
        divisor0 = beta_tmp
      }else{
        divisor0 = (alpha_tmp-1)*beta_tmp
      }
      dd = ceiling(5*h / divisor0)
      ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
      subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
      
      mean_trj1 = mean_trajectory_NFL_delay_feed(ext_tspan, c(MCMC_sample[rep-1,1], MCMC_sample[rep-1, 2], MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4], MCMC_sample[rep-1, 5], MCMC_sample[rep-1, 6]))[subseq0]
      mean_trj1_st = mean_trajectory_NFL_delay_feed(ext_tspan, c(theta1_cand, MCMC_sample[rep-1, 2], MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4], MCMC_sample[rep-1, 5], MCMC_sample[rep-1, 6]))[subseq0]
      
      log_lik = sum(log_of_normpdf(data, mean_trj1, sqrt(var_list)))
      log_lik_st = sum(log_of_normpdf(data, mean_trj1_st, sqrt(var_list)))
      
      acceptance_ratio1 = exp(log_pri_st - log_pri + log_lik_st - log_lik)
      
      if (runif(1) < acceptance_ratio1){
        MCMC_sample[rep, 1] = theta1_cand
        update_matrix[rep,1] = 1
      }else{
        MCMC_sample[rep, 1] = MCMC_sample[rep-1,1]
      }
    }else{
      MCMC_sample[rep, 1] = MCMC_sample[rep-1, 1]
    } 
    
    # MCMC for theta(2): decay rate of X(t);
    # Random walk Metropolis algorithm (use normal proposal, symmetric)
    
    prop_mean2 = MCMC_sample[rep-1,2]
    prop_var2 = prop_var[2] # tune this parameter
    
    theta2_cand = rnorm(1, mean = prop_mean2, sd = sqrt(prop_var2))
    
    if (theta2_cand > 0){
      gamma_pri_alpha_2 = prior_mean[2]^2 / prior_var[2]
      gamma_pri_beta_2 = prior_var[2] / prior_mean[2]
      
      log_pri = log_of_gampdf(MCMC_sample[rep-1, 2], gamma_pri_alpha_2, gamma_pri_beta_2)
      log_pri_st = log_of_gampdf(theta2_cand, gamma_pri_alpha_2, gamma_pri_beta_2)
      
      alpha_tmp = MCMC_sample[rep-1, 3]
      beta_tmp = MCMC_sample[rep-1, 4]
      h = timespan[2] - timespan[1]
      if(alpha_tmp == 1){
        divisor0 = beta_tmp
      }else{
        divisor0 = (alpha_tmp-1)*beta_tmp
      }
      dd = ceiling(5*h / divisor0)
      ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
      subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
      
      mean_trj2 = mean_trajectory_NFL_delay_feed(ext_tspan, c(MCMC_sample[rep,1], MCMC_sample[rep-1, 2], MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4], MCMC_sample[rep-1, 5], MCMC_sample[rep-1, 6]))[subseq0]
      mean_trj2_st = mean_trajectory_NFL_delay_feed(ext_tspan, c(MCMC_sample[rep,1], theta2_cand, MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4], MCMC_sample[rep-1, 5], MCMC_sample[rep-1, 6]))[subseq0]
      
      log_lik = sum(log_of_normpdf(data, mean_trj2, sqrt(var_list)))
      log_lik_st = sum(log_of_normpdf(data, mean_trj2_st, sqrt(var_list)))
      
      acceptance_ratio2 = exp(log_pri_st - log_pri + log_lik_st - log_lik)
      
      if (runif(1) < acceptance_ratio2){
        MCMC_sample[rep, 2] = theta2_cand
        update_matrix[rep,2] = 1
      }else{
        MCMC_sample[rep, 2] = MCMC_sample[rep-1,2]
      }
    }else{
      MCMC_sample[rep, 2] = MCMC_sample[rep-1, 2]
    }
    
    # MCMC for theta(3): the shape parameter of delay distribution
    
    # MCMC_sample(rep,3) = MCMC_sample(rep-1,3);
    
    prop_mean3 = MCMC_sample[rep-1,3];
    prop_var3 = prop_var[3]; # tune this parameter
    
    theta3_cand = rnorm(1, mean = prop_mean3, sd = sqrt(prop_var3));
    if (theta3_cand >= 1){
      gamma_pri_alpha_3 = prior_mean[3]^2 / prior_var[3]
      gamma_pri_beta_3 = prior_var[3]/prior_mean[3] 
      
      log_pri = log_of_gampdf(MCMC_sample[rep-1, 3], gamma_pri_alpha_3, gamma_pri_beta_3)
      log_pri_st = log_of_gampdf(theta3_cand, gamma_pri_alpha_3, gamma_pri_beta_3)
      
      
      
      alpha_tmp = MCMC_sample[rep-1, 3]
      beta_tmp = MCMC_sample[rep-1, 4]
      h = timespan[2] - timespan[1]
      if(alpha_tmp == 1){
        divisor0 = beta_tmp
      }else{
        divisor0 = (alpha_tmp-1)*beta_tmp
      }
      dd = ceiling(5*h / divisor0)
      ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
      subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
      
      mean_trj3 = mean_trajectory_NFL_delay_feed(ext_tspan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep-1, 3], MCMC_sample[rep-1, 4], MCMC_sample[rep-1, 5], MCMC_sample[rep-1, 6]))[subseq0]
      
      
      alpha_tmp = theta3_cand
      beta_tmp = MCMC_sample[rep-1, 4]
      h = timespan[2] - timespan[1]
      if(alpha_tmp == 1){
        divisor0 = beta_tmp
      }else{
        divisor0 = (alpha_tmp-1)*beta_tmp
      }
      dd = ceiling(5*h / divisor0)
      ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
      subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
      mean_trj3_st = mean_trajectory_NFL_delay_feed(ext_tspan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], theta3_cand, MCMC_sample[rep-1, 4], MCMC_sample[rep-1, 5], MCMC_sample[rep-1, 6]))[subseq0]
      
      log_lik = sum(log_of_normpdf(data, mean_trj3, sqrt(var_list)))
      log_lik_st = sum(log_of_normpdf(data, mean_trj3_st, sqrt(var_list)))
      
      acceptance_ratio3 = exp(log_pri_st - log_pri + log_lik_st - log_lik)
      
      if (runif(1) < acceptance_ratio3){
        MCMC_sample[rep, 3] = theta3_cand
        update_matrix[rep, 3] = 1
      }else{
        MCMC_sample[rep, 3] = MCMC_sample[rep-1,3]
      }
    }else{
      MCMC_sample[rep, 3] = MCMC_sample[rep-1, 3]
    }
    
    # === MCMC for theta(4): the rate parameter of delay distribution ===
    
    #     MCMC_sample(rep,4) = MCMC_sample(rep-1,4);
    
    
    prop_mean4 = MCMC_sample[rep-1,4];
    prop_var4 = prop_var[4]; # tune this parameter
    
    theta4_cand = rnorm(1, mean = prop_mean4, sd = sqrt(prop_var4))
    if (theta4_cand > 0){
      gamma_pri_alpha_4 = prior_mean[4]^2 / prior_var[4]
      gamma_pri_beta_4 = prior_var[4]/prior_mean[4]
      
      log_pri = log_of_gampdf(MCMC_sample[rep-1, 4], gamma_pri_alpha_4, gamma_pri_beta_4)
      log_pri_st = log_of_gampdf(theta4_cand, gamma_pri_alpha_4, gamma_pri_beta_4)
      
      
      alpha_tmp = MCMC_sample[rep, 3]
      beta_tmp = MCMC_sample[rep-1, 4]
      h = timespan[2] - timespan[1]
      if(alpha_tmp == 1){
        divisor0 = beta_tmp
      }else{
        divisor0 = (alpha_tmp-1)*beta_tmp
      }
      dd = ceiling(5*h / divisor0)
      ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
      subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
      
      mean_trj4 = mean_trajectory_NFL_delay_feed(ext_tspan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep, 3], MCMC_sample[rep-1, 4], MCMC_sample[rep-1, 5], MCMC_sample[rep-1, 6]))[subseq0]
      
      alpha_tmp = MCMC_sample[rep, 3]
      beta_tmp = theta4_cand
      h = timespan[2] - timespan[1]
      if(alpha_tmp == 1){
        divisor0 = beta_tmp
      }else{
        divisor0 = (alpha_tmp-1)*beta_tmp
      }
      dd = ceiling(5*h / divisor0)
      ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
      subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
      
      mean_trj4_st = mean_trajectory_NFL_delay_feed(ext_tspan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep, 3], theta4_cand, MCMC_sample[rep-1, 5], MCMC_sample[rep-1, 6]))[subseq0]
      
      
      log_lik = sum(log_of_normpdf(data, mean_trj4, sqrt(var_list)))
      log_lik_st = sum(log_of_normpdf(data, mean_trj4_st, sqrt(var_list)))
      
      acceptance_ratio4 = exp(log_pri_st - log_pri + log_lik_st - log_lik)
      
      if (runif(1) < acceptance_ratio4){
        MCMC_sample[rep, 4] = theta4_cand
        update_matrix[rep, 4] = 1
      }else{
        MCMC_sample[rep, 4] = MCMC_sample[rep-1, 4]
      }
    }else{
      MCMC_sample[rep, 4] = MCMC_sample[rep-1, 4]
    }
    
    
    # === MCMC for theta(5): the feedback sensitivity ===
    
    #     MCMC_sample(rep,5) = MCMC_sample(rep-1,5);
    
    
    prop_mean5 = MCMC_sample[rep-1,5];
    prop_var5 = prop_var[5]; # tune this parameter
    
    theta5_cand = rnorm(1, mean = prop_mean5, sd = sqrt(prop_var5))
    if (theta5_cand > 0){
      gamma_pri_alpha_5 = prior_mean[5]^2 / prior_var[5]
      gamma_pri_beta_5 = prior_var[5]/prior_mean[5]
      
      log_pri = log_of_gampdf(MCMC_sample[rep-1, 5], gamma_pri_alpha_5, gamma_pri_beta_5)
      log_pri_st = log_of_gampdf(theta5_cand, gamma_pri_alpha_5, gamma_pri_beta_5)
      
      
      alpha_tmp = MCMC_sample[rep, 3]
      beta_tmp = MCMC_sample[rep, 4]
      h = timespan[2] - timespan[1]
      if(alpha_tmp == 1){
        divisor0 = beta_tmp
      }else{
        divisor0 = (alpha_tmp-1)*beta_tmp
      }
      dd = ceiling(5*h / divisor0)
      ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
      subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
      
      mean_trj5 = mean_trajectory_NFL_delay_feed(ext_tspan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep, 3], MCMC_sample[rep, 4], MCMC_sample[rep-1, 5], MCMC_sample[rep-1, 6]))[subseq0]
      mean_trj5_st = mean_trajectory_NFL_delay_feed(ext_tspan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep, 3], MCMC_sample[rep, 4], theta5_cand, MCMC_sample[rep-1, 6]))[subseq0]
      
      log_lik = sum(log_of_normpdf(data, mean_trj5, sqrt(var_list)))
      log_lik_st = sum(log_of_normpdf(data, mean_trj5_st, sqrt(var_list)))
      
      acceptance_ratio5 = exp(log_pri_st - log_pri + log_lik_st - log_lik)
      
      if (runif(1) < acceptance_ratio5){
        MCMC_sample[rep, 5] = theta5_cand
        update_matrix[rep, 5] = 1
      }else{
        MCMC_sample[rep, 5] = MCMC_sample[rep-1, 5]
      }
    }else{
      MCMC_sample[rep, 5] = MCMC_sample[rep-1, 5]
    }
    
    
    # === MCMC for theta(6): the feedback sensitivity ===
    
    #     MCMC_sample(rep,6) = MCMC_sample(rep-1,6);
    
    
    prop_mean6 = MCMC_sample[rep-1,6];
    prop_var6 = prop_var[6]; # tune this parameter
    
    theta6_cand = rnorm(1, mean = prop_mean6, sd = sqrt(prop_var6))
    if (theta6_cand > 0){
      gamma_pri_alpha_6 = prior_mean[6]^2 / prior_var[6]
      gamma_pri_beta_6 = prior_var[6]/prior_mean[6]
      
      log_pri = log_of_gampdf(MCMC_sample[rep-1, 6], gamma_pri_alpha_6, gamma_pri_beta_6)
      log_pri_st = log_of_gampdf(theta6_cand, gamma_pri_alpha_6, gamma_pri_beta_6)
      
      
      alpha_tmp = MCMC_sample[rep, 3]
      beta_tmp = MCMC_sample[rep, 4]
      h = timespan[2] - timespan[1]
      if(alpha_tmp == 1){
        divisor0 = beta_tmp
      }else{
        divisor0 = (alpha_tmp-1)*beta_tmp
      }
      dd = ceiling(5*h / divisor0)
      ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
      subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
      
      mean_trj6 = mean_trajectory_NFL_delay_feed(ext_tspan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep, 3], MCMC_sample[rep, 4], MCMC_sample[rep, 5], MCMC_sample[rep-1, 6]))[subseq0]
      mean_trj6_st = mean_trajectory_NFL_delay_feed(ext_tspan, c(MCMC_sample[rep,1], MCMC_sample[rep, 2], MCMC_sample[rep, 3], MCMC_sample[rep, 4], MCMC_sample[rep, 5],theta6_cand))[subseq0]
      
      log_lik = sum(log_of_normpdf(data, mean_trj6, sqrt(var_list)))
      log_lik_st = sum(log_of_normpdf(data, mean_trj6_st, sqrt(var_list)))
      
      acceptance_ratio6 = exp(log_pri_st - log_pri + log_lik_st - log_lik)
      
      if (runif(1) < acceptance_ratio6){
        MCMC_sample[rep, 6] = theta6_cand
        update_matrix[rep, 6] = 1
      }else{
        MCMC_sample[rep, 6] = MCMC_sample[rep-1, 6]
      }
    }else{
      MCMC_sample[rep, 6] = MCMC_sample[rep-1, 6]
    }
    
  }
  results = MCMC_sample;
  acceptance = update_matrix;
  my_list <- list("results" = results, "acceptance" = acceptance)
  return(my_list)
}



MCMC_function_NFL_ME_delay_feed_fine <- function(data, timespan, var_list, hyper_prior_mean = rep(1,12), hyper_prior_var = rep(1000,12), 
                                                 theta_init, theta_hyper_init, burn = 0, jump = 1, effnum = 100, 
                                                 prop_var = c(50, 0.000001, 0.1, 0.003, 50, 1), hyper_prop_var = c(50,1,0.001,1,1,1,0.01,1,1,1,1,1),
                                                 save.file.name, summary.file.name){
  
  # In this function, data should be a matrix with more than one column. Each column represents a timeseries data.
  # var_list is also a matrix with the size identical to data.
  
  num_repeat = burn + jump * effnum;
  num_data = dim(data)[2]
  num_param = length(theta_init)
  
  A_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  B_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  beta_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  AR_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  tau2_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  
  A_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  B_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  beta_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  AR_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  tau2_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  
  MCMC_sample_hyp = matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  update_matrix_hyp =  matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  
  A_sample[1,] = theta_init[1]
  B_sample[1,] = theta_init[2]
  alpha_sample[1,] = theta_init[3]
  beta_sample[1,] = theta_init[4]
  AR_sample[1,] = theta_init[5]
  tau2_sample[1,] = theta_init[6]
  
  MCMC_sample_hyp[1,] = theta_hyper_init
  
  thetaCurrent <- matrix(NA, nrow = num_data, ncol = num_param)
  thetaCurrent[,1] <- theta_init[1]
  thetaCurrent[,2] <- theta_init[2]
  thetaCurrent[,3] <- theta_init[3]
  thetaCurrent[,4] <- theta_init[4]
  thetaCurrent[,5] <- theta_init[5]
  thetaCurrent[,6] <- theta_init[6]
  
  thetaCand <- thetaCurrent
  
  thetaHypCurrent <- theta_hyper_init
  
  thetaHypCand <- thetaHypCurrent
  
  for (rep in 2:num_repeat){
    if (rep %% 1000 == 0){
      cat(rep); cat("\n");
      cat(summary.file.name); cat("\n");
    }
    
    ## Start to sample the parameters for each individual
    for (dataid in 1:num_data){
      for (param_id in 1:num_param){
        param_cand <- rnorm(1, mean = thetaCurrent[dataid,param_id], sd = sqrt(prop_var[param_id]))
        if (param_cand > 0){
          if(param_id == 3 & param_cand < 1){
            # thetaCand[, param_id] = param_cand
            next
          }else{
            thetaCand[dataid, param_id] = param_cand
          }  
          
          
          
          alpha_tmp = thetaCurrent[dataid,3]
          beta_tmp = thetaCurrent[dataid,4]
          h = timespan[2] - timespan[1]
          if(alpha_tmp == 1){
            divisor0 = beta_tmp
          }else{
            divisor0 = (alpha_tmp-1)*beta_tmp
          }
          dd = ceiling(5*h / divisor0)
          ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
          subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
          
          mean_trj_current <- mean_trajectory_NFL_delay_feed(ext_tspan, thetaCurrent[dataid,])[subseq0]
          
          alpha_tmp = thetaCand[dataid,3]
          beta_tmp = thetaCand[dataid,4]
          h = timespan[2] - timespan[1]
          if(alpha_tmp == 1){
            divisor0 = beta_tmp
          }else{
            divisor0 = (alpha_tmp-1)*beta_tmp
          }
          dd = ceiling(5*h / divisor0)
          ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
          subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
          
          mean_trj_cand <- mean_trajectory_NFL_delay_feed(ext_tspan, thetaCand[dataid,])[subseq0]
          
          
          log_pri = log_of_gampdf(thetaCurrent[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          log_pri_cand = log_of_gampdf(thetaCand[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          
          log_lik = sum(log_of_normpdf(data[,dataid], mean_trj_current, sqrt(var_list[,dataid])))
          log_lik_cand = sum(log_of_normpdf(data[,dataid], mean_trj_cand, sqrt(var_list[,dataid])))
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaCurrent = thetaCand
            switch(param_id, 
                   A_update[rep,dataid] <- 1, 
                   B_update[rep,dataid] <- 1, 
                   alpha_update[rep,dataid] <- 1, 
                   beta_update[rep,dataid] <- 1,
                   AR_update[rep,dataid] <- 1,
                   tau2_update[rep,dataid] <- 1)
          }else{
            thetaCand = thetaCurrent
          }
        }
      }
    }
    
    ## Sample the shared alpha values for all individual
    
    ## Start to sample the hyper parameters for all populations.
    for (hyp_param_id in 1:(2*num_param)){
      
      hyp_param_cand <- rnorm(1, mean = thetaHypCurrent[hyp_param_id], sd = sqrt(hyper_prop_var[hyp_param_id]))
      
      if (hyp_param_cand > 0){
        thetaHypCand[hyp_param_id] = hyp_param_cand
        
        gamma_pri_alpha = hyper_prior_mean[hyp_param_id]^2 / hyper_prior_var[hyp_param_id]
        gamma_pri_beta = hyper_prior_var[hyp_param_id]/hyper_prior_mean[hyp_param_id]
        
        log_pri = log_of_gampdf(thetaHypCurrent[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
        log_pri_cand = log_of_gampdf(thetaHypCand[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
        
        if (hyp_param_id %% 2 == 1){
          log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id], num_data), rep(thetaHypCurrent[hyp_param_id+1], num_data)))
          log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id], num_data), rep(thetaHypCand[hyp_param_id+1], num_data)))
        }else{
          log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id-1], num_data), rep(thetaHypCurrent[hyp_param_id], num_data)))
          log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id-1], num_data), rep(thetaHypCand[hyp_param_id], num_data)))
        }
        
        acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
        
        if (runif(1) < acceptance_ratio){
          thetaHypCurrent = thetaHypCand
          update_matrix_hyp[rep, hyp_param_id] = 1
        }else{
          thetaHypCand = thetaHypCurrent
        }
      }
      
    }
    A_sample[rep,] = thetaCurrent[,1]
    B_sample[rep,] = thetaCurrent[,2]
    alpha_sample[rep,] = thetaCurrent[,3]
    beta_sample[rep,] = thetaCurrent[,4]
    AR_sample[rep,] = thetaCurrent[,5]
    tau2_sample[rep,] = thetaCurrent[,6]
    MCMC_sample_hyp[rep, ] = thetaHypCurrent
    
    
    if (rep %% 5500 == 0 | rep == num_repeat){
      mcmc.result1 <- list("A_sample" = A_sample, "B_sample" = B_sample, "alpha_sample" = alpha_sample, "beta_sample" = beta_sample, "AR_sample" = AR_sample, "tau2_sample" = tau2_sample,
                           "A_update" = A_update, "B_update" = B_update, "alpha_update" = alpha_update, "beta_update" = beta_update, "AR_update" = AR_update, "tau2_update" = tau2_update,
                           "MCMC_sample_hyp" = MCMC_sample_hyp, "update_matrix_hyp" = update_matrix_hyp, "rep" = rep)
      inputdata1 = data
      # save(mcmc.result1, inputdata1, file = paste(save.file.name, ".Rdata", sep = ""))
    }
  }
  return(mcmc.result1)
}


MCMC_function_NFL_ME_delay_feed_share_alpha_fine <- function(data, timespan, var_list, hyper_prior_mean = rep(1,12), hyper_prior_var = rep(1000,12), 
                                                             theta_init, theta_hyper_init, burn = 0, jump = 1, effnum = 100, 
                                                             prop_var = c(50, 0.000001, 0.1, 0.003, 50, 1), hyper_prop_var = c(50,1,0.001,1,1,1,0.01,1,1,1,1,1),
                                                             save.file.name, summary.file.name){
  
  # In this function, data should be a matrix with more than one column. Each column represents a timeseries data.
  # var_list is also a matrix with the size identical to data.
  
  num_repeat = burn + jump * effnum;
  num_data = dim(data)[2]
  num_param = length(theta_init)
  
  A_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  B_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  beta_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  AR_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  tau2_sample = matrix(0, nrow = num_repeat, ncol = num_data)
  
  A_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  B_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  alpha_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  beta_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  AR_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  tau2_update =  matrix(0, nrow = num_repeat, ncol = num_data)
  
  MCMC_sample_hyp = matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  update_matrix_hyp =  matrix(0, nrow = num_repeat, ncol = length(theta_hyper_init))
  
  A_sample[1,] = theta_init[1]
  B_sample[1,] = theta_init[2]
  alpha_sample[1,] = theta_init[3]
  beta_sample[1,] = theta_init[4]
  AR_sample[1,] = theta_init[5]
  tau2_sample[1,] = theta_init[6]
  
  MCMC_sample_hyp[1,] = theta_hyper_init
  
  thetaCurrent <- matrix(NA, nrow = num_data, ncol = num_param)
  thetaCurrent[,1] <- theta_init[1]
  thetaCurrent[,2] <- theta_init[2]
  thetaCurrent[,3] <- theta_init[3]
  thetaCurrent[,4] <- theta_init[4]
  thetaCurrent[,5] <- theta_init[5]
  thetaCurrent[,6] <- theta_init[6]
  
  thetaCand <- thetaCurrent
  
  thetaHypCurrent <- theta_hyper_init
  
  thetaHypCand <- thetaHypCurrent
  
  for (rep in 2:num_repeat){
    if (rep %% 1000 == 0){
      cat(rep); cat("\n");
      cat(summary.file.name); cat("\n");
    }
    
    ## Start to sample the parameters for each individual
    for (dataid in 1:num_data){
      for (param_id in 1:num_param){
        param_cand <- rnorm(1, mean = thetaCurrent[dataid,param_id], sd = sqrt(prop_var[param_id]))
        if (param_cand > 0){
          if(param_id == 3){
            # thetaCand[, param_id] = param_cand
            next
          }else{
            thetaCand[dataid, param_id] = param_cand
          }  
          
          
          
          alpha_tmp = thetaCurrent[dataid,3]
          beta_tmp = thetaCurrent[dataid,4]
          h = timespan[2] - timespan[1]
          if(alpha_tmp == 1){
            divisor0 = beta_tmp
          }else{
            divisor0 = (alpha_tmp-1)*beta_tmp
          }
          dd = ceiling(5*h / divisor0)
          ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
          subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
          
          mean_trj_current <- mean_trajectory_NFL_delay_feed(ext_tspan, thetaCurrent[dataid,])[subseq0]
          
          alpha_tmp = thetaCand[dataid,3]
          beta_tmp = thetaCand[dataid,4]
          h = timespan[2] - timespan[1]
          if(alpha_tmp == 1){
            divisor0 = beta_tmp
          }else{
            divisor0 = (alpha_tmp-1)*beta_tmp
          }
          dd = ceiling(5*h / divisor0)
          ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
          subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
          
          mean_trj_cand <- mean_trajectory_NFL_delay_feed(ext_tspan, thetaCand[dataid,])[subseq0]
          
          
          log_pri = log_of_gampdf(thetaCurrent[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          log_pri_cand = log_of_gampdf(thetaCand[dataid,param_id], thetaHypCurrent[2*param_id-1], thetaHypCurrent[2*param_id])
          
          log_lik = sum(log_of_normpdf(data[,dataid], mean_trj_current, sqrt(var_list[,dataid])))
          log_lik_cand = sum(log_of_normpdf(data[,dataid], mean_trj_cand, sqrt(var_list[,dataid])))
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaCurrent = thetaCand
            switch(param_id, 
                   A_update[rep,dataid] <- 1, 
                   B_update[rep,dataid] <- 1, 
                   alpha_update[rep,dataid] <- 1, 
                   beta_update[rep,dataid] <- 1,
                   AR_update[rep,dataid] <- 1,
                   tau2_update[rep,dataid] <- 1)
          }else{
            thetaCand = thetaCurrent
          }
        }
      }
    }
    
    ## Sample the shared alpha values for all individual
    
    param_cand <- rnorm(1, mean = thetaCurrent[1,3], sd = sqrt(prop_var[3]))
    
    if (param_cand >= 1){
      thetaCand[, 3] = param_cand
      log_pri = log_of_gampdf(thetaCurrent[1,3], thetaHypCurrent[5], thetaHypCurrent[6])
      log_pri_cand = log_of_gampdf(thetaCand[1,3], thetaHypCurrent[5], thetaHypCurrent[6])
      log_lik = 0
      log_lik_cand = 0
      for (dataid in 1:num_data){
        
        
        alpha_tmp = thetaCurrent[dataid,3]
        beta_tmp = thetaCurrent[dataid,4]
        h = timespan[2] - timespan[1]
        if(alpha_tmp == 1){
          divisor0 = beta_tmp
        }else{
          divisor0 = (alpha_tmp-1)*beta_tmp
        }
        dd = ceiling(5*h / divisor0)
        ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
        subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
        
        mean_trj_current <- mean_trajectory_NFL_delay_feed(ext_tspan, thetaCurrent[dataid,])[subseq0]
        
        alpha_tmp = thetaCand[dataid,3]
        beta_tmp = thetaCand[dataid,4]
        h = timespan[2] - timespan[1]
        if(alpha_tmp == 1){
          divisor0 = beta_tmp
        }else{
          divisor0 = (alpha_tmp-1)*beta_tmp
        }
        dd = ceiling(5*h / divisor0)
        ext_tspan = seq(from = h/dd, to = max(timespan), by = h/dd)
        subseq0 = round(seq(from = dd, to = dd * max(timespan)/h, by = dd))
        
        mean_trj_cand <- mean_trajectory_NFL_delay_feed(ext_tspan, thetaCand[dataid,])[subseq0]
        
        log_lik = log_lik + sum(log_of_normpdf(data[,dataid], mean_trj_current, sqrt(var_list[,dataid])))
        log_lik_cand = log_lik_cand + sum(log_of_normpdf(data[,dataid], mean_trj_cand, sqrt(var_list[,dataid])))
      }
      acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
      
      if (runif(1) < acceptance_ratio){
        thetaCurrent = thetaCand
        switch(3, 
               A_update[rep,dataid] <- 1, 
               B_update[rep,dataid] <- 1, 
               alpha_update[rep,dataid] <- 1, 
               beta_update[rep,dataid] <- 1,
               AR_update[rep,dataid] <- 1,
               tau2_update[rep,dataid] <- 1)
      }else{
        thetaCand = thetaCurrent
      }
    }
    
    ## Start to sample the hyper parameters for all populations.
    for (hyp_param_id in 1:(2*num_param)){
      if (hyp_param_id == 5 | hyp_param_id == 6){
        thetaHypCand = thetaHypCurrent
      }else{
        hyp_param_cand <- rnorm(1, mean = thetaHypCurrent[hyp_param_id], sd = sqrt(hyper_prop_var[hyp_param_id]))
        
        if (hyp_param_cand > 0){
          thetaHypCand[hyp_param_id] = hyp_param_cand
          
          gamma_pri_alpha = hyper_prior_mean[hyp_param_id]^2 / hyper_prior_var[hyp_param_id]
          gamma_pri_beta = hyper_prior_var[hyp_param_id]/hyper_prior_mean[hyp_param_id]
          
          log_pri = log_of_gampdf(thetaHypCurrent[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
          log_pri_cand = log_of_gampdf(thetaHypCand[hyp_param_id], gamma_pri_alpha, gamma_pri_beta)
          
          if (hyp_param_id %% 2 == 1){
            log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id], num_data), rep(thetaHypCurrent[hyp_param_id+1], num_data)))
            log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id], num_data), rep(thetaHypCand[hyp_param_id+1], num_data)))
          }else{
            log_lik = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCurrent[hyp_param_id-1], num_data), rep(thetaHypCurrent[hyp_param_id], num_data)))
            log_lik_cand = sum(log_of_gampdf(thetaCurrent[,ceiling(hyp_param_id/2)], rep(thetaHypCand[hyp_param_id-1], num_data), rep(thetaHypCand[hyp_param_id], num_data)))
          }
          
          acceptance_ratio = exp(log_pri_cand - log_pri + log_lik_cand - log_lik)
          
          if (runif(1) < acceptance_ratio){
            thetaHypCurrent = thetaHypCand
            update_matrix_hyp[rep, hyp_param_id] = 1
          }else{
            thetaHypCand = thetaHypCurrent
          }
        }
      }
    }
    A_sample[rep,] = thetaCurrent[,1]
    B_sample[rep,] = thetaCurrent[,2]
    alpha_sample[rep,] = thetaCurrent[,3]
    beta_sample[rep,] = thetaCurrent[,4]
    AR_sample[rep,] = thetaCurrent[,5]
    tau2_sample[rep,] = thetaCurrent[,6]
    MCMC_sample_hyp[rep, ] = thetaHypCurrent
    
    
    if (rep %% 5500 == 0 | rep == num_repeat){
      mcmc.result1 <- list("A_sample" = A_sample, "B_sample" = B_sample, "alpha_sample" = alpha_sample, "beta_sample" = beta_sample, "AR_sample" = AR_sample, "tau2_sample" = tau2_sample,
                           "A_update" = A_update, "B_update" = B_update, "alpha_update" = alpha_update, "beta_update" = beta_update, "AR_update" = AR_update, "tau2_update" = tau2_update,
                           "MCMC_sample_hyp" = MCMC_sample_hyp, "update_matrix_hyp" = update_matrix_hyp, "rep" = rep)
      inputdata1 = data
      # save(mcmc.result1, inputdata1, file = paste(save.file.name, ".Rdata", sep = ""))
    }
    
  }
  return(mcmc.result1)
}



