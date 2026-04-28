#### NEW PROJECT
### LATE PHASE DOSE FINDING
# Update

library(Rcpp) # Load package 'Rcpp'
library(MyPackageV2) # contains RCPP code for parallilisation 
library(rBeta2009)
library(microbenchmark)
library(gridExtra)

library(ggplot2)

library(doParallel)  # Parallelisation
library(foreach)  # Parallelisation

library(survival)
library(survminer)


library(dplyr)
library(knitr)
library(kableExtra)

library(gt)
library(simtrial)

library(cmprsk)





##### Code ##### 

# Collapse for cleanliness
if (TRUE){

  qq_inv_function = function(q, n = 160, method = 1){
    
    if (method == 1){
      
      count = 0
      
      temp_prob = 1 - rcpp_exact_beta(1 + count, 1+ n - count, 1 + round(0.96*n), 1 + round(0.04*n)) # 1 + because Beta(1,1) prior
      
      temp_prop_plus_1 = temp_prob
      
      while((temp_prob - q)*(temp_prop_plus_1 - q) > 0){
        count = count + 1
        if (count == n){
          return(n)
        }
        
        temp_prob = temp_prop_plus_1
        
        temp_prop_plus_1 = 1 - rcpp_exact_beta(1 + count, 1+ n - count, 1 + round(0.96*n), 1 + round(0.04*n))
      }
      
      return(count)
    }
    
    if (method == 2){
      
      count = 0
      
      temp_prob = norm_diff_cdf(1 + count, 1+ n - count, 1 + round(0.96*n), 1 + round(0.04*n) , c = 0 )
      
      temp_prop_plus_1 = temp_prob
      
      while((temp_prob - q)*(temp_prop_plus_1 - q) > 0){
        count = count + 1
        if (count == n){
          return(n)
        }
        
        temp_prob = temp_prop_plus_1
        
        temp_prop_plus_1 = norm_diff_cdf(1 + count, 1+ n - count, 1 + round(0.96*n), 1 + round(0.04*n) , c = 0 )
      }
      
      return(count)
    }
    
    if (method == 3){
      
      count = 0
      
      temp_prob = beta_sim(1 + count, 1+ n - count, 1 + round(0.96*n), 1 + round(0.04*n) , c = 0 ) # 1 + because Beta(1,1) prior
      
      temp_prop_plus_1 = temp_prob
      
      while((temp_prob - q)*(temp_prop_plus_1 - q) > 0){
        count = count + 1
        if (count == n){
          return(n)
        }
        
        temp_prob = temp_prop_plus_1
        
        temp_prop_plus_1 = beta_sim(1 + count, 1+ n - count, 1 + round(0.96*n), 1 + round(0.04*n) , c = 0 ) # 1 + because Beta(1,1) prior
      }
      
      return(count)
    }
    
    
  }
  
  safer_allo = function(data_treat, data_con, safe_prob, eta = 1){
    # Initially built for surivial endpoints but the CLT still makes it hold
    # safe prob is just the ratio of safe treat over safe treat + safe con
    theta_treat = mean(data_treat)
    theta_con = mean(data_con)
    
    var_treat = var(data_treat)
    var_con = var(data_con)
    var_total = var_treat + var_con
    if(var_total == 0){
      return(safe_prob) # error catching in case
    }
    
    Z = (theta_treat - theta_con)/sqrt(var_total)
    phi = pnorm(abs(Z))
    
    if (phi <= 0.5){
      return(0.5)
    }
    else if (phi < 1){
      return(0.5 + (safe_prob - 0.5)*(1 - (1 - ((phi - 0.5)/0.5)^eta)))
    }
    else{
      return(safe_prob)
    }
    
    
  }
  
  log_exact_beta = function(A1,B1,A2,B2){
    # Does Miller's exact beta calculation in log space
    #
    
    prob = 0
    
    for (i in 0:(A2-1)){
      to_add_log = lbeta(A1 + i, B1 + B2) - log(B2 + i) - lbeta(i+1,B2) - lbeta(A1,B1)
      prob = prob + exp(to_add_log)
    }
    
    return(prob)
    
  }
  
  beta_sim = function(a_1,b_1,a_2,b_2,c = 0.05, n_sim = 10000){
    prob = 0
    
    return(sum(rbeta(n_sim,a_1,b_1) > rbeta(n_sim,a_2,b_2) - c)/n_sim)
    
  }
  
  var_beta = function(a,b){
    # variance of a beta(a,b)
    
    return( a*b/(((a+b)^2)*(a+b+1)) )
  }
  
  mean_beta = function(a,b){
    # mean of a beta(a,b)
    if (a + b == 2){
      if (a == b){
        return(0.5) # then it's just a uniform
      }
    }
    return( (a - 1)/(a + b - 2))
    
  }
  
  norm_diff_cdf = function(a_1,b_1,a_2,b_2,c = 0.05){
    # does a normal approximation of the Prob(theta_1 > theta_0 - c | D)
    
    mu_1 = mean_beta(a_1,b_1)
    mu_2 = mean_beta(a_2,b_2)
    
    var_1 = var_beta(a_1,b_1)
    var_2 = var_beta(a_2,b_2)
    
    return(pnorm(c, mu_2 - mu_1, sqrt(var_1 + var_2))) # Probability that the treatment is at least c less than the control
    
    
  }
  
  
  NI_trial = function(eff_treat, tox_treat, eff_con, tox_con, N, IAn = N, NI_margin = 0.06, conf = 0.99, test = 1, alpha = 0.05){
    ## Will compare the treatment and control in an NI trial context
    # eff_treat and tox_treat are the treatments efficacy and toxicity
    # eff_con and tox_con are the control's treatment and efficacy
    # N is maximum sample size
    # IAn is a vector of interims placements (choose divisible by 2 for ease)
    # NI margin is the level of clinically acceptable reduction in efficacy
    # conf is a vector (matching the length of IAn) of thresholds to compare the posterior probabilities against
    # beta 1,1 prior
    
    # test = 1 is standard NI posterior test, adjusting for the NI margin
    # test = 2 is posterior that doesn't adjust for the NI
    # test = 3 is the Bayes factor
    # test = 4 is normal approximation, alpha is only used here to calculate the width of the CI
    
    eff_treat_vec = numeric(0)
    tox_treat_vec = numeric(0)
    eff_con_vec = numeric(0)
    tox_con_vec = numeric(0)
    
    for (i in 1:length(IAn)){  
      if (i == 1){ # how many to simulate
        
        to_sim = IAn[1]/2
  
      }
      else{
        to_sim = (IAn[i] - IAn[i-1])/2
      }
      
      # simulating
      
      eff_treat_vec = append(eff_treat_vec, rbinom(to_sim,1,eff_treat))
      tox_treat_vec = append(tox_treat_vec, rbinom(to_sim,1,tox_treat))
      eff_con_vec = append(eff_con_vec, rbinom(to_sim,1,eff_con))
      tox_con_vec = append(tox_con_vec, rbinom(to_sim,1,tox_con))
          
      
      eff_treat_suc = sum(eff_treat_vec) # successes in the treatment
      eff_treat_fail = IAn[i]/2 - eff_treat_suc # failures in the treatment
      
      eff_con_suc = sum(eff_con_vec) # successes in the treatment
      eff_con_fail = IAn[i]/2 - eff_con_suc # failures in the treatment
      
      
      # analysis
      if (test == 1){
        if (NI_margin == 0){
          # do rcpp beta
          
          post_prob = rcpp_exact_beta(1 + eff_treat_suc, 1+eff_treat_fail, 1 + eff_con_suc, 1 + eff_con_fail) # 1 + because Beta(1,1) prior
          
        }
        else{
          # do the sim beta
          
          post_prob = beta_sim(1 + eff_treat_suc, 1+eff_treat_fail, 1 + eff_con_suc, 1 + eff_con_fail , c = NI_margin ) # 1 + because Beta(1,1) prior
        }
        
        if (post_prob > conf[i]){ # if we pass the decision boundary
          return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
        } 
      }
      else if (test == 3){
        
        post_prob = bayes_factor_comp(1 + eff_treat_suc, 1+eff_treat_fail, 1 + eff_con_suc, 1 + eff_con_fail) # 1 + because Beta(1,1) prior
  
        
        if (post_prob < conf[i]){ # if we pass the decision boundary
          return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
        } 
      }
      
      else if (test == 4){
        
        # assuming normal then finding 95% approximate CI 
        num_IA = length(IAn)
        Pocock_bounds = numeric(num_IA)
        for (j in 1:num_IA){
          Pocock_bounds[j] = sqrt(num_IA/j)
        }
        
        
        n_treat = eff_treat_suc + eff_treat_fail
        mle_treat = eff_treat_suc/n_treat
        var_treat = (eff_treat_suc*eff_treat_fail) / n_treat^3
        
        n_con = eff_con_suc + eff_con_fail
        mle_con = eff_con_suc/n_con
        var_con = (eff_con_suc*eff_con_fail) / n_con^3
        
        var_diff = var_con + var_treat
        
        # Now, mle_treat - mle_con ~ N(delta, var_diff) under normal assumptions, where delta is the NI margin
        # So, mle_treat - mle_con +- 1.96*sqrt(var_diff) should be a 95% CI
        # If mle_treat - mle_con - 1.96*var_diff > delta, then we have NI 
        
        CI_conf = qnorm(1 - alpha)*Pocock_bounds[i]
        
        if (mle_con - mle_treat + CI_conf*sqrt(var_diff) < NI_margin){ # if we pass the decision boundary # NOTE not currently tuned for multiple IA 
          return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
        } 
      }
      
    }
    
    return(c(0,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,N)) # else if we never pass we fail
    
  }
  
  NI_trial_RAR = function(eff_treat, tox_treat, eff_con, tox_con, N, IAn = N, NI_margin = 0.06, conf = c(0.99,0.98,0.97,0.96), test = 1, alpha = 0.05, RAR = 0, w = 0.5, sigma = 0.1, fut_stop = FALSE, safe_stop = FALSE, safe_conf = rep(0.2,length(IAn)), norm_approx = TRUE){
    ## Will compare the treatment and control in an NI trial context
    # eff_treat and tox_treat are the treatments efficacy and toxicity
    # eff_con and tox_con are the control's treatment and efficacy
    # N is maximum sample size
    # IAn is a vector of interims placements (choose divisible by 2 for ease)
    # NI margin is the level of clinically acceptable reduction in efficacy
    # conf is a vector (matching the length of IAn) of thresholds to compare the posterior probabilities against
    # beta 1,1 prior
    # w is the weighting for RITS. w close to 1 faovurs the efficacy, 0 favours the safety
    # sigma tunes gaussian RITS
    # fut_stop = TRUE means there is futility stopping (stochastic Curtailment approach)
    
    # test = 1 is standard NI posterior test, adjusting for the NI margin
    # test = 2 is posterior that doesn't adjust for the NI
    # test = 3 is the Bayes factor
    # test = 4 is normal approximation, alpha is only used here to calculate the width of the CI
    
    # RAR = 0 is just equal randomisation
    # RAR = 1 is RITS
    # RAR = 2 is Gaussian RITS
    # RAR = 3 is safety
    # RAR = 4 is just standard efficacy based BRAR
    
    allo_prob = 0.5 # starts at fixed, just defining here to avoid errors later
    
    eff_treat_vec = numeric(0)
    tox_treat_vec = numeric(0)
    eff_con_vec = numeric(0)
    tox_con_vec = numeric(0)
    
    for (i in 1:length(IAn)){  
      if (i == 1){ # how many to simulate
        
        to_sim = IAn[1]/2 # burn in, no RAR
        to_sim_con = to_sim
        to_sim_treat = to_sim
        
      }
      else{
        
        if (RAR == 0){
          # no RAR
          allo_prob = 0.5
        }
        else if (RAR == 1){
          # RITS
          tox_treat_1 = sum(tox_treat_vec)
          tox_treat_0 = length(tox_treat_vec) - tox_treat_1 
          
          tox_con_1 = sum(tox_con_vec)
          tox_con_0 = length(tox_con_vec) - tox_con_1
        
          
          allo_prob_eff = 1 - rcpp_exact_beta(1 + eff_treat_suc,1 + eff_treat_fail,1 + eff_con_suc,1 + eff_con_fail)
          allo_prob_tox = rcpp_exact_beta(1 + tox_treat_1, 1 + tox_treat_0, 1 + tox_con_1, 1 + tox_con_0)
          
          allo_prob = w*allo_prob_eff + (1 - w)*allo_prob_tox
          
          allo_prob = allo_prob^(IAn[i-1]/(2*N))/(allo_prob^(IAn[i-1]/(2*N)) + (1-allo_prob)^(IAn[i-1]/(2*N))) # tuning with n/2N
          
          
        }
        else if (RAR == 2){
          # Gaussian RITS
          
          tox_treat_1 = sum(tox_treat_vec)
          tox_treat_0 = length(tox_treat_vec) - tox_treat_1 
          
          tox_con_1 = sum(tox_con_vec)
          tox_con_0 = length(tox_con_vec) - tox_con_1
          
          
          allo_prob_eff = 1 - rcpp_exact_beta(1 + eff_treat_suc,1 + eff_treat_fail,1 + eff_con_suc,1 + eff_con_fail)
          allo_prob_tox = rcpp_exact_beta(1 + tox_treat_1, 1 + tox_treat_0, 1 + tox_con_1, 1 + tox_con_0)
          
          
          g_x = exp( - ((allo_prob_eff - 0.5)^2) / sigma )
          
          allo_prob = (1 - g_x)*allo_prob_eff + g_x*allo_prob_tox
          
          allo_prob = allo_prob^(IAn[i-1]/(2*N))/(allo_prob^(IAn[i-1]/(2*N)) + (1-allo_prob)^(IAn[i-1]/(2*N))) # tuning with n/2N
          
          
        }
        else if (RAR == 3){
          # SAFER
          tox_treat_1 = sum(tox_treat_vec)
          tox_treat_0 = length(tox_treat_vec) - tox_treat_1 
          
          tox_con_1 = sum(tox_con_vec)
          tox_con_0 = length(tox_con_vec) - tox_con_1
          
          allo_prob_tox = rcpp_exact_beta(1 + tox_treat_1, 1 + tox_treat_0, 1 + tox_con_1, 1 + tox_con_0)
          
          allo_prob = safer_allo(eff_treat_vec,eff_con_vec, allo_prob_tox) 
          # will try without tuning, not expecting SAFER to be hugely relevant for now
          
        }
        else if (RAR == 4){
          # BRAR just based on efficacy
          allo_prob = 1 - rcpp_exact_beta(1 + eff_treat_suc,1 + eff_treat_fail,1 + eff_con_suc,1 + eff_con_fail)
          
          allo_prob = allo_prob^(IAn[i-1]/(2*N))/(allo_prob^(IAn[i-1]/(2*N)) + (1-allo_prob)^(IAn[i-1]/(2*N))) # tuning with n/2N
        }
        
        to_sim = (IAn[i] - IAn[i-1])
        
        to_sim_con = round(to_sim*(1 - allo_prob))
        to_sim_treat = to_sim -  to_sim_con
        
        
      }
      
      
      
      
      # simulating
      #print(eff_treat)
      if(is.na(to_sim_treat)){
        browser()
      } # error catching
  
      if (to_sim_treat < 0){
        browser()
      }    
    
      eff_treat_vec = append(eff_treat_vec, rbinom(to_sim_treat,1,eff_treat))  
      tox_treat_vec = append(tox_treat_vec, rbinom(to_sim_treat,1,tox_treat))
      eff_con_vec = append(eff_con_vec, rbinom(to_sim_con,1,eff_con))
      tox_con_vec = append(tox_con_vec, rbinom(to_sim_con,1,tox_con))
      
      # Doing these to try and avoid appending to see if it's faster
      #new_eff_treat_vec = rbinom(to_sim_treat,1,eff_treat)
      #eff_treat_vec = c(eff_treat_vec, new_eff_treat_vec) 
      
      #new_tox_treat_vec = rbinom(to_sim_treat,1,tox_treat)
      #tox_treat_vec = c(tox_treat_vec, new_tox_treat_vec) 
      
      #new_eff_con_vec = rbinom(to_sim_con,1,eff_con)
      #eff_con_vec = c(eff_con_vec, new_eff_con_vec) 
      
      #new_tox_con_vec = rbinom(to_sim_con,1,tox_con)
      #tox_con_vec = c(tox_con_vec, new_tox_con_vec) 
      # not quicker, just going back
      
      eff_treat_suc = sum(eff_treat_vec) # successes in the treatment
      eff_treat_fail = length(eff_treat_vec) - eff_treat_suc # failures in the treatment
      
      eff_con_suc = sum(eff_con_vec) # successes in the treatment
      eff_con_fail = length(eff_con_vec) - eff_con_suc # failures in the treatment
      
      #browser()
      
      # analysis
      if (test == 1){
        if (NI_margin == 0){
          # do rcpp beta
          
          post_prob = rcpp_exact_beta(1 + eff_treat_suc, 1+eff_treat_fail, 1 + eff_con_suc, 1 + eff_con_fail) # 1 + because Beta(1,1) prior
          
        }
        else{
          # do the sim beta
          if (norm_approx){
            post_prob = norm_diff_cdf(1 + eff_treat_suc, 1+eff_treat_fail, 1 + eff_con_suc, 1 + eff_con_fail , c = NI_margin )
          }
          else{
            #browser()
            post_prob = beta_sim(1 + eff_treat_suc, 1+eff_treat_fail, 1 + eff_con_suc, 1 + eff_con_fail , c = NI_margin ) # 1 + because Beta(1,1) prior
          }
        }
        
        
        if (post_prob > conf[i]){ # if we pass the decision boundary
          return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
          #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
        } 
        else if (fut_stop){
          eff_treat_suc_curtail = eff_treat_suc + round(eff_con*(N - IAn[i])/2) # if it does very well
          eff_con_suc_curtail = eff_con_suc + round(eff_con*(N - IAn[i])/2) # if it does as expected
          
          eff_treat_fail_curtail = eff_treat_fail + round((1 - eff_con)*(N - IAn[i])/2) # if it does very well
          eff_con_fail_curtail = eff_con_fail + round((1 - eff_con)*(N - IAn[i])/2) # if it does as expected
          
          if (norm_approx){
            post_prob_curtail = norm_diff_cdf(1 + eff_treat_suc_curtail, 1+eff_treat_fail_curtail, 1 + eff_con_suc_curtail, 1 + eff_con_fail_curtail , c = NI_margin ) # 1 + because Beta(1,1) prior
          }
          else{
            
            post_prob_curtail = beta_sim(1 + eff_treat_suc_curtail, 1+eff_treat_fail_curtail, 1 + eff_con_suc_curtail, 1 + eff_con_fail_curtail , c = NI_margin ) # 1 + because Beta(1,1) prior
          
          }
          # If the probability of it rejecting the null at the end is unlikely, then will just end the trial early
          if (post_prob_curtail < conf[length(IAn)] - 0.1){
            #browser()
            return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            
          }
          
        }
        if (safe_stop){
          
          tox_treat_suc = sum(tox_treat_vec) # successes in the treatment
          tox_treat_fail = length(tox_treat_vec) - tox_treat_suc # failures in the treatment
          
          tox_con_suc = sum(tox_con_vec) # successes in the treatment
          tox_con_fail = length(tox_con_vec) - tox_con_suc # failures in the treatment
          
          safe_post_prob = rcpp_exact_beta(1 + tox_treat_suc, 1+tox_treat_fail, 1 + tox_con_suc, 1 + tox_con_fail) # 1 + because Beta(1,1) prior
          
          if(safe_post_prob < safe_conf[i]){
            return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            
          }
          
        }
      }
      else if (test == 3){
        
        post_prob = bayes_factor_comp(1 + eff_treat_suc, 1+eff_treat_fail, 1 + eff_con_suc, 1 + eff_con_fail, norm_approx = norm_approx) # 1 + because Beta(1,1) prior
        
        #browser()
        
        if (post_prob < conf[i]){ # if we pass the decision boundary
          return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
          #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
        } 
        else if (fut_stop){
          eff_treat_suc_curtail = eff_treat_suc + round(eff_con*(N - IAn[i])/2) # if it does very well
          eff_con_suc_curtail = eff_con_suc + round(eff_con*(N - IAn[i])/2) # if it does as expected
          
          eff_treat_fail_curtail = eff_treat_fail + round((1 - eff_con)*(N - IAn[i])/2) # if it does very well
          eff_con_fail_curtail = eff_con_fail + round((1 - eff_con)*(N - IAn[i])/2) # if it does as expected
          
          post_prob_curtail = bayes_factor_comp(1 + eff_treat_suc_curtail, 1+eff_treat_fail_curtail, 1 + eff_con_suc_curtail, 1 + eff_con_fail_curtail, norm_approx = norm_approx)
  
          # If the probability of it rejecting the null at the end is unlikely, then will just end the trial early
          if (post_prob_curtail > conf[length(IAn)] + 0.1){
            return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            
          }
        }
        if (safe_stop){
          
          tox_treat_suc = sum(tox_treat_vec) # successes in the treatment
          tox_treat_fail = length(tox_treat_vec) - tox_treat_suc # failures in the treatment
          
          tox_con_suc = sum(tox_con_vec) # successes in the treatment
          tox_con_fail = length(tox_con_vec) - tox_con_suc # failures in the treatment
          
          safe_post_prob = rcpp_exact_beta(1 + tox_treat_suc, 1+tox_treat_fail, 1 + tox_con_suc, 1 + tox_con_fail) # 1 + because Beta(1,1) prior
          
          if(safe_post_prob < safe_conf[i]){
            return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            
          }
        }
      }
      
      else if (test == 4){
        
        # assuming normal then finding 95% approximate CI 
        num_IA = length(IAn)
        Pocock_bounds = numeric(num_IA)
        for (j in 1:num_IA){
          Pocock_bounds[j] = sqrt(num_IA/j)
        }
        
        
        n_treat = eff_treat_suc + eff_treat_fail
        mle_treat = eff_treat_suc/n_treat
        var_treat = (eff_treat_suc*eff_treat_fail) / n_treat^3
        
        n_con = eff_con_suc + eff_con_fail
        mle_con = eff_con_suc/n_con
        var_con = (eff_con_suc*eff_con_fail) / n_con^3
        
        var_diff = var_con + var_treat
        
        # Now, mle_treat - mle_con ~ N(delta, var_diff) under normal assumptions, where delta is the NI margin
        # So, mle_treat - mle_con +- 1.96*sqrt(var_diff) should be a 95% CI
        # If mle_treat - mle_con - 1.96*var_diff > delta, then we have NI 
        
        CI_conf = qnorm(1 - alpha)*Pocock_bounds[i]
        
        if (mle_con - mle_treat + CI_conf*sqrt(var_diff) < NI_margin){ # if we pass the decision boundary # NOTE not currently tuned for multiple IA 
          return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
        }
        else if (fut_stop){
          eff_treat_suc_curtail = eff_treat_suc + round(eff_con*(N - IAn[i])/2) # if it does very well
          eff_con_suc_curtail = eff_con_suc + round(eff_con*(N - IAn[i])/2) # if it does as expected
          
          eff_treat_fail_curtail = eff_treat_fail + round((1 - eff_con)*(N - IAn[i])/2) # if it does very well
          eff_con_fail_curtail = eff_con_fail + round((1 - eff_con)*(N - IAn[i])/2) # if it does as expected
          
          
          n_treat_c = eff_treat_suc_curtail + eff_treat_fail_curtail
          mle_treat_c  = eff_treat_suc_curtail/n_treat_c
          var_treat_c = (eff_treat_suc_curtail*eff_treat_fail_curtail) / n_treat_c^3
          
          n_con_c = eff_con_suc_curtail + eff_con_fail_curtail
          mle_con_C = eff_con_suc_curtail/n_con_c
          var_con_c = (eff_con_suc_curtail*eff_con_fail_curtail) / n_con_c^3
          
          var_diff_c = var_con_c + var_treat_c
          
          # Now, mle_treat - mle_con ~ N(delta, var_diff) under normal assumptions, where delta is the NI margin
          # So, mle_treat - mle_con +- 1.96*sqrt(var_diff) should be a 95% CI
          # If mle_treat - mle_con - 1.96*var_diff > delta, then we have NI 
          
          CI_conf = qnorm(1 - alpha)*Pocock_bounds[length(IAn)]
          
          #browser()
  
          # If the probability of it rejecting the null at the end is unlikely, then will just end the trial early
          if (mle_con_C - mle_treat_c + CI_conf*sqrt(var_diff_c) > NI_margin + 0.01){ # if we pass the decision boundary # NOTE not currently tuned for multiple IA 
            return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
          }
        }
        if (safe_stop){
          
          tox_treat_suc = sum(tox_treat_vec) # successes in the treatment
          tox_treat_fail = length(tox_treat_vec) - tox_treat_suc # failures in the treatment
          
          tox_con_suc = sum(tox_con_vec) # successes in the treatment
          tox_con_fail = length(tox_con_vec) - tox_con_suc # failures in the treatment
          
          safe_post_prob = rcpp_exact_beta(1 + tox_treat_suc, 1+tox_treat_fail, 1 + tox_con_suc, 1 + tox_con_fail) # 1 + because Beta(1,1) prior
          
          if(safe_post_prob < safe_conf[i]){
            return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            
          }
        }
      }
    }
    
    return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
    
  }
  
  
  OC_gen = function(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 1000, IAn = N, NI_margin = 0.06, conf = 0.99, n_sim = 10000, power = FALSE, safety_supe = FALSE, conf_tox = 0.99, conf_eff = 0.99, test = 1, alpha = 0.05, RAR = 0, fut_stop = FALSE, safe_stop = FALSE, safe_conf = rep(0.2,length(IAn)), norm_approx = TRUE, paired_seed = FALSE, set_seed = 2026){
    # finds the operating characterestics of the NI trial
    if (power == TRUE){
      power = 0
      for (i in 1:n_sim){
        
        if (paired_seed){
          set.seed(set_seed + i)
        }
        
        if (safety_supe){
          power = supe_safety_trial(eff_treat, tox_treat, eff_con, tox_con, N, IAn, conf_tox = conf_tox, conf_eff = conf_eff)[1] + power
        }
        else{
          power = NI_trial(eff_treat, tox_treat, eff_con, tox_con, N, IAn, NI_margin, conf, test = test, alpha = alpha)[1] + power
        }
      }
      power = power/n_sim
      return(power)
    }
    else{
      if (safety_supe){
        OC_vec = numeric(7)
      }
      else{
        OC_vec = numeric(8)
      }
      for (i in 1:n_sim){
        
        if (paired_seed){
          set.seed(set_seed + i)
        }
        
        if (safety_supe){
          OC_vec = supe_safety_trial(eff_treat, tox_treat, eff_con, tox_con, N, IAn, conf_tox = conf_tox, conf_eff = conf_eff) + OC_vec
        }
        else{
          OC_vec = NI_trial_RAR(eff_treat, tox_treat, eff_con, tox_con, N, IAn, NI_margin, conf, test = test, alpha = alpha, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx) + OC_vec
        }
      }
      OC_vec = OC_vec/n_sim
      return(OC_vec)
    }
  }
  
  # Then, want to make a plot which shows the sample size required to get 81% power with varying margins
  
  threshold_finder = function(tox_treat, eff_con, tox_con, N, IAn = N, NI_margin = 0.06, n_sim = 10000, alpha = 0.05, pre_conf = numeric(0), RAR = 0, fut_stop = FALSE, safe_stop = FALSE, safe_conf = rep(0.2,length(IAn)), norm_approx = TRUE){
    # finds the threshold that gives a certain type one error
    
    # If IAn is a vector longer than 1  (ex there are IA), then pre_conf is the confidence thresholds before this, so this function just optimises for one threhsold
    
    if (length(IAn) == 1){
      
      conf = 0.9
      while (OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim, conf = conf, power = FALSE, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1] > alpha){ # will continue looping until type one error acceptable
        conf = (4*conf + 1)/5 # can change the weighted ratio to speed up (or slow down for more accuracy)
        #browser()
        #print(conf)
      }
      
      upper = (5*conf - 1)/4
      lower = conf
      conf = (upper+lower)/2
      for (i in 1:5){
        # final bit of refinement to try and get as close to 0.05 type one as possible
        
        if (OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim, conf = conf, power = FALSE, test = 1, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1] > alpha){
          # Type one error too high
          # Gotta increase boundary
          #print(conf)
          
          
          lower = conf
          conf = (upper + conf)/2
          
        }
        else{
          # type one error low
          # gotta increase boundary
          
          upper = conf
          conf = (lower + conf)/2 
          
        }
      
        #print(conf)
      }
      return(conf)
    
    }
    else{
      conf = 0.9
      while (OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim, conf = c(pre_conf,conf), power = FALSE, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1] > alpha){ # will continue looping until type one error acceptable
        conf = (4*conf + 1)/5 # can change the weighted ratio to speed up (or slow down for more accuracy)
        #browser()
      }
      
      upper = (5*conf - 1)/4
      lower = conf
      conf = (upper+lower)/2
      for (i in 1:5){
        # final bit of refinement to try and get as close to 0.05 type one as possible
        
        if (OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim, conf = c(pre_conf,conf), power = FALSE, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1] > alpha){
          # Type one error too high
          # Gotta increase boundary
          
          
          
          lower = conf
          conf = (upper + conf)/2
          
        }
        else{
          # type one error low
          # gotta increase boundary
          
          upper = conf
          conf = (lower + conf)/2 
          
        }
      }
      
      
      return(c(pre_conf,conf))
    }
  }
  # conf = 0.959 controls type one error at 0.05
  
  threshold_finder_bisection = function(tox_treat, eff_con, tox_con, N, IAn = N, NI_margin = 0.06, n_sim = 10000, alpha = 0.05, pre_conf = numeric(0), RAR = 0, fut_stop = FALSE, safe_stop = FALSE, safe_conf = rep(0.2,length(IAn)), norm_approx = TRUE, iter = 11, recheck_factor = 2){
    # aims to improve on the above optimisation by using by considering bisection
    
    lower_conf = 0.9
    conf = 0.95
    upper_conf = 1
    
    error = 1.96*(sqrt(0.25/n_sim)) # roughly 95% MC error bound
    
    for (i in 1:iter){
      
      temp_alpha = OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim, c(pre_conf,conf), power = FALSE, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1] 
      
      if (temp_alpha > alpha + error){ # error too high, gotta make boundary more strict
        
        lower_conf = conf
        conf = (conf + upper_conf)/2
      }
      else if (temp_alpha < alpha - error){ # error too low, gotta make boundary less strict
        
        upper_conf = conf
        conf = (conf + lower_conf)/2
        
      }
      else{
        temp_alpha = OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim*recheck_factor, c(pre_conf,conf), power = FALSE, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1] 
        
        if (temp_alpha > alpha){ # error too high, gotta make boundary more strict
          
          lower_conf = conf
          conf = (conf + upper_conf)/2
        }
        else{ # error too low, gotta make boundary less strict
          
          upper_conf = conf
          conf = (conf + lower_conf)/2
          
        }
  
      }
    }
    
    
    return(c(pre_conf,conf))
    
  }
  
  threshold_finder_bayes_factor = function(tox_treat, eff_con, tox_con, N, IAn = N, NI_margin = 0.06, n_sim = 10000, alpha = 0.05, pre_conf = numeric(0), RAR = 0, fut_stop = FALSE, safe_stop = FALSE, safe_conf = rep(0.2,length(IAn)), norm_approx = TRUE){
    # finds the threshold that gives a certain type one error
    if (length(IAn) == 1){
    
      conf = 1
      while (OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim, conf = conf, power = FALSE, test = 3, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1] > alpha){ # will continue looping until type one error acceptable
        conf = (3*conf)/4 # can change the weighted ratio to speed up (or slow down for more accuracy)
        #browser()
        #print(conf)
      }
      upper = 4*conf/3
      lower = conf
      conf = (upper+lower)/2
      for (i in 1:5){
        # final bit of refinement to try and get as close to 0.05 type one as possible
        
        if (OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim, conf = conf, power = FALSE, test = 3, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1] > alpha){
          # Type one error too high
          # Gotta decrease boundary
          
          upper = conf
          conf = (lower + conf)/2 
          
        }
        else{
          # type one error low
          # gotta increase boundary
          
          lower = conf
          conf = (upper + conf)/2
          
        }
        
        
      }
      
      return(conf)
    
    }
    else{
      conf = 1
      while (OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim, conf = c(pre_conf,conf), power = FALSE, test = 3, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1] > alpha){ # will continue looping until type one error acceptable
        conf = (3*conf)/4 # can change the weighted ratio to speed up (or slow down for more accuracy)
        #browser()
        #print(conf)
      }
      
      upper = 4*conf/3
      lower = conf
      conf = (upper+lower)/2
      for (i in 1:5){
        # final bit of refinement to try and get as close to 0.05 type one as possible
        
        if (OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim, conf = c(pre_conf,conf), power = FALSE, test = 3, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1] > alpha){
          # Type one error too high
          # Gotta decrease boundary
          
          upper = conf
          conf = (lower + conf)/2 
          
        }
        else{
          # type one error low
          # gotta increase boundary
          
          lower = conf
          conf = (upper + conf)/2
          
        }
        
        
      }
      
      
      return(c(pre_conf,conf))
    }
    
    
    
  }
  # 0.316 is controls type one error and gets 90% power for N = 500
  # 0.133 is the boundary for N = 300, gets type one to 0.038
  threshold_finder_bayes_bisection = function(tox_treat, eff_con, tox_con, N, IAn = N, NI_margin = 0.06, n_sim = 10000, alpha = 0.05, pre_conf = numeric(0), RAR = 0, fut_stop = FALSE, safe_stop = FALSE, safe_conf = rep(0.2,length(IAn)), norm_approx = TRUE, iter = 11, recheck_factor = 2){
    # aims to improve on the above optimisation by using by considering bisection
    # should only really do this if around N = 320 as that what this covers
    
    # set up a while loop first to get the right region
    
    
    
    lower_conf = 0.1
    upper_conf = 1
    
    error = 1.96*(sqrt(0.25/n_sim)) # roughly 95% MC error bound
    
    temp_alpha = OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim, c(pre_conf,lower_conf), power = FALSE, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx, test = 3)[1] 
    
    while(temp_alpha > alpha + error ){ # gets it into the right magnitude 
      
      upper_conf = lower_conf
      lower_conf = lower_conf/10
      temp_alpha = OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim, c(pre_conf,upper_conf), power = FALSE, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx, test = 3)[1] 
      
      
    }
    
    conf = (lower_conf + upper_conf)/2
    
    
    for (i in 1:iter){
      
      temp_alpha = OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim, c(pre_conf,conf), power = FALSE, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx, test = 3)[1] 
      
      if (temp_alpha > alpha + error){ # error too high, gotta make boundary more strict
        # for Bayes factor, this means it going down
  
        upper_conf = conf
        conf = (conf + lower_conf)/2
      }
      else if (temp_alpha < alpha - error){ # error too low, gotta make boundary less strict
        
        lower_conf = conf
        conf = (conf + upper_conf)/2
        
        
        
      }
      else{
        temp_alpha = OC_gen(eff_con - NI_margin, tox_treat, eff_con, tox_con, N, IAn, NI_margin, n_sim = n_sim*recheck_factor, c(pre_conf,conf), power = FALSE, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx, test = 3)[1] 
        
        if (temp_alpha > alpha){ # error too high, gotta make boundary more strict
          
          
          
          upper_conf = conf
          conf = (conf + lower_conf)/2
        }
        else{ # error too low, gotta make boundary less strict
          
          lower_conf = conf
          conf = (conf + upper_conf)/2
          
        }
        
      }
    }
    
    
    return(c(pre_conf,conf))
    
  }
  
  
  
  ss_finder = function(tox_treat, eff_con, tox_con, NI_margin = 0.06, alpha = 0.05, power = 0.81, n_sim = 1000, current_ss = 100, test = 1, fut_stop = FALSE, safe_stop = FALSE, RAR = 0){
    # finds the sample size that gives target type one error and power (with no IA)
    
    #current_ss = 1000
    while(TRUE){
      
      # tox_treat = tox_con here as we don't want an overly generous H0
      temp_conf = threshold_finder_v2_no_safety(tox_treat = tox_con, eff_con, tox_con, n = current_ss, IAn = current_ss, NI_margin = NI_margin, n_sim = n_sim, alpha = alpha, test = test, fut_stop = fut_stop, safe_stop = safe_stop, RAR = RAR )
      #temp_conf = threshold_finder(tox_treat, eff_con, tox_con, N = current_ss, IAn = current_ss, NI_margin = NI_margin, n_sim = n_sim, alpha = alpha)
      temp_power = OC_gen(eff_con, tox_treat, eff_con, tox_con, N = current_ss, IAn = current_ss, NI_margin, n_sim = n_sim, conf = temp_conf, power = TRUE, test = test, fut_stop = fut_stop, safe_stop = safe_stop, RAR = RAR, alpha = temp_conf)
      if (temp_power > power){
        return(current_ss)
      }
      else{
        if (temp_power < 0.6*power){
          # big jump
          current_ss = current_ss + 200
        }
        else if (temp_power < 0.8*power){
          # middle jump jump
          current_ss = current_ss + 100
        }
        else{
          # small jump
          current_ss = current_ss + 50
        }
      }
      
    }
  }
  
  ss_finder_v2 = function(){
    
    poss_ss_vec = seq(320,1920,80)
  
    
    num_iter =  ceiling( log(length(poss_ss_vec), base = 2) )
    
    high_index = length(poss_ss_vec)
    low_index = 1
    
    
    
    for (i in 1:num_iter){
      
      temp_index = (low_index + high_index)/2
      
      temp_conf = threshold_finder_v2_no_safety(tox_treat = tox_con, eff_con, tox_con, n = poss_ss_vec[temp_index], IAn = poss_ss_vec[temp_index], NI_margin = NI_margin, n_sim = n_sim, alpha = alpha, test = test, fut_stop = fut_stop, safe_stop = safe_stop, RAR = RAR )
      #temp_conf = threshold_finder(tox_treat, eff_con, tox_con, N = current_ss, IAn = current_ss, NI_margin = NI_margin, n_sim = n_sim, alpha = alpha)
      temp_power = OC_gen(eff_con, tox_treat, eff_con, tox_con, N = poss_ss_vec[temp_index], IAn = poss_ss_vec[temp_index], NI_margin, n_sim = n_sim, conf = temp_conf, power = TRUE, test = test, fut_stop = fut_stop, safe_stop = safe_stop, RAR = RAR, alpha = temp_conf)
      
      
      if (temp_power > 1.05*power){ # power too high, can bring the sample size down
        
        high_index = temp_index
        
      }
      else if (temp_power < 0.95*power){ # power too low, can bring sample size up
        
        low_index = temp_index
      
      }
      
      
    }
    
    temp_index = (low_index + high_index)/2
    
    return(poss_ss_vec[temp_index])
    
    
  }
  
  plot_ss = function(tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, alpha = 0.05, power = 0.80, n_sim = 1000){
    NI_margin_vec = 3:8/100
    ss_vec = numeric(length(NI_margin_vec))
    
    for (i in 1:length(NI_margin_vec)){
      ss_vec[i] = ss_finder(tox_treat, eff_con, tox_con, alpha = alpha, power = power, n_sim = 1000, NI_margin = NI_margin_vec[i])
      print(i)
    }
    
    #browser()
    plot(x = NI_margin_vec, y = ss_vec,type = "l", xlab = "NI Margin", ylab = "Sample Size", main = "Sample size required for 80% power with varying NI margins.")
    
    #browser()
    
  }
  
  supe_safety_trial = function(eff_treat, tox_treat, eff_con, tox_con, N, IAn = N, conf_tox = 0.99, conf_eff = 0.99){
    ## The idea of this trial is a standard supe safety trial that has early stopping for efficacy
    
    eff_treat_vec = numeric(0)
    tox_treat_vec = numeric(0)
    eff_con_vec = numeric(0)
    tox_con_vec = numeric(0)
    
    for (i in 1:length(IAn)){  
      if (i == 1){ # how many to simulate
        
        to_sim = IAn[1]/2
        
      }
      else{
        to_sim = (IAn[i] - IAn[i-1])/2
      }
      
      # simulating
      
      eff_treat_vec = append(eff_treat_vec, rbinom(to_sim,1,eff_treat))
      tox_treat_vec = append(tox_treat_vec, rbinom(to_sim,1,tox_treat))
      eff_con_vec = append(eff_con_vec, rbinom(to_sim,1,eff_con))
      tox_con_vec = append(tox_con_vec, rbinom(to_sim,1,tox_con))
      
      
      eff_treat_suc = sum(eff_treat_vec) # successes in the treatment
      eff_treat_fail = IAn[i]/2 - eff_treat_suc # failures in the treatment
      
      eff_con_suc = sum(eff_con_vec) # successes in the treatment
      eff_con_fail = IAn[i]/2 - eff_con_suc # failures in the treatment
      
      tox_treat_suc = sum(tox_treat_vec) # toxics
      tox_treat_fail = IAn[i]/2 - tox_treat_suc # non-toxics
      
      tox_con_suc = sum(tox_con_vec) # toxics
      tox_con_fail = IAn[i]/2 - tox_con_suc # non toxics
      
      # plus 1 as beta(1,1) priors
      p_eff = rcpp_exact_beta(eff_treat_suc + 1,eff_treat_fail + 1,eff_con_suc + 1, eff_con_fail + 1)
      p_tox = rcpp_exact_beta(tox_treat_suc + 1,tox_treat_fail + 1,tox_con_suc + 1, tox_con_fail + 1)
      
      if(is.na(p_eff)){
        browser()
      }
     
      if (p_eff > conf_eff[i]){
        return_vec = c(FALSE, eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail, i, 101)
        return(return_vec)  # ENDS EARLY DUE TO FUTILITY
      }
      else if (p_tox > conf_tox[i]){
  
        return_vec = c(TRUE, eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail, i, 102)
        return(return_vec)  # ENDS EARLY FOR SUPE OF SAFETY (AND ACCEPTABLE EFF)
      }
       
    }
    
    
    # o/f fails because of not significance
    
    return_vec = c(FALSE, eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail, i, 103)
    return(return_vec)  # ENDS EARLY due to non significance
  }
  
  
  #microbenchmark(supe_safety_trial(0.8,0.2,0.9,0.5,1000))
  
  #microbenchmark(NI_trial(0.8,0.2,0.9,0.5,1000))
  
  threshold_safety_finder = function(alpha = 0.05, n = 1000, IAn = c(n), n_sim = 10000, alpha_spending = seq(0.1,1,length.out = (length(IAn))), pregiven_tox = numeric(0), pregiven_eff = numeric(0)){
    # Now, this is more complicated than the previous threshold finder as it's a grid search in 2D
    # But, as it the power is monotonic in both, we can be a bit more clever as we know the gradient sign
    
    # pregiven tox and pregiven eff are for the recursive calling
    
    if (length(IAn) == 1){
      # can do easy stuff by only considering one IA
      conf_eff_vec = seq(0.98,0.99,0.001)
      
      max_power = 0
      current_threshold = c(0,0)
      
      for (i in 1:length(conf_eff_vec)){
        
        print(i)
        
        temp_conf_tox = 0
        temp_calc = OC_gen(safety_supe = TRUE, eff_treat = 0.9, tox_treat = 0.5, conf_eff = conf_eff_vec[i], conf_tox = temp_conf_tox, N = n, IAn = IAn, n_sim = n_sim)
        if (temp_calc[1] < alpha){ # If with the most forgiving conf it fails to be over the type one error, then we're done
          temp_power = OC_gen(safety_supe = TRUE, eff_treat = 0.96, conf_eff = conf_eff_vec[i], conf_tox = temp_conf_tox, N = n, IAn = IAn, n_sim = n_sim)[1]
          if (temp_power > max_power){
            max_power = temp_power
            current_threshold = c(conf_eff_vec[i],temp_conf_tox)
          }
        }
        else{
          # So, we're above type one error, now to keep and as temp_conf_tox = 1 gives 0 type one error, we should be able to get close to 0.05 (remembering we're dealing with binary data)
          temp_conf_tox = 0.5
          upper = 1
          lower = 0 #lower and upper are what to halve to for the halving method
          
          temp_calc = OC_gen(safety_supe = TRUE, eff_treat = 0.9, tox_treat = 0.5, conf_eff = conf_eff_vec[i], conf_tox = temp_conf_tox, N = n, IAn = IAn, n_sim = n_sim)
          same_counter = 0
          while (same_counter < 3){ # when the result hasn't changed in 3 jumps, aka convergence.
            if (temp_calc[1] > alpha){
              # if our type one error is too high and we need conf_tox to be higher
              lower = temp_conf_tox # we could simply not record these?
              temp_conf_tox = (temp_conf_tox + upper)/2
              
            }
            else{
              # if our type one errr is too low and we need conf_tox to be lower
              upper = temp_conf_tox
              temp_conf_tox = (temp_conf_tox + lower)/2
              
              
            }
            
            new_calc = OC_gen(safety_supe = TRUE, eff_treat = 0.9, tox_treat = 0.5, conf_eff = conf_eff_vec[i], conf_tox = temp_conf_tox, N = n, IAn = IAn, n_sim = n_sim)
            if (new_calc[1] == temp_calc[1]){ ## WILL ONLY TERMINATE THIS WAY BECAUSE OF BINARY DATA
              same_counter = same_counter + 1
            }
            else{
              same_counter = 0
            }
            if (abs(new_calc[1] - alpha) < 0.0001){ #  or if it's close enough to 0.5
              same_counter = 3 # just end the loops there 
            }
            temp_calc = new_calc
            
          }
          
          temp_power = OC_gen(safety_supe = TRUE, eff_treat = 0.96, conf_eff = conf_eff_vec[i], conf_tox = temp_conf_tox, N = n, IAn = IAn, n_sim = n_sim)[1]
          if (temp_power > max_power){
            max_power = temp_power
            current_threshold = c(conf_eff_vec[i],temp_conf_tox)
          }
          
        }
        
      }
    
    
      return(c("Power" = max_power, "Efficacy conf" = current_threshold[1], "Toxicity conf" = current_threshold[2]))
    }
    else{
      # Could do recursive? Just call it again and again until we reach it
      # also we use alpha spending here, to denote how much alpha we want to spend at each stage
      
      if (length(pregiven_tox) == length(IAn) - 1){ # if we're ready to find the new threshold (as pre)
        conf_eff_vec = seq(0.98,0.99,0.001)
        
        max_power = 0
        current_threshold = c(0,0)
        
        for (i in 1:length(conf_eff_vec)){
          
          print(i)
          
          #browser()
          
          temp_conf_tox = c(pregiven_tox, 0)
          temp_conf_eff = c(pregiven_eff, conf_eff_vec[i])
          
          #temp_conf_tox = 0
          temp_calc = OC_gen(safety_supe = TRUE, eff_treat = 0.9, tox_treat = 0.5, conf_eff = temp_conf_eff, conf_tox = temp_conf_tox, N = n, IAn = IAn, n_sim = n_sim)
          if (temp_calc[1] < alpha){ # If with the most forgiving conf it fails to be over the type one error, then we're done
            temp_power = OC_gen(safety_supe = TRUE, eff_treat = 0.96, conf_eff = temp_conf_eff, conf_tox = temp_conf_tox, N = n, IAn = IAn, n_sim = n_sim)[1]
            if (temp_power > max_power){
              max_power = temp_power
              current_threshold = c(temp_conf_eff,temp_conf_tox)
            }
          }
          else{
            # So, we're above type one error, now to keep and as temp_conf_tox = 1 gives 0 type one error, we should be able to get close to 0.05 (remembering we're dealing with binary data)
            temp_conf_tox = c(pregiven_tox, 0.5)
            upper = c(pregiven_tox, 1)
            lower = c(pregiven_tox, 0) #lower and upper are what to halve to for the halving method
            
            temp_calc = OC_gen(safety_supe = TRUE, eff_treat = 0.9, tox_treat = 0.5, conf_eff = temp_conf_eff, conf_tox = temp_conf_tox, N = n, IAn = IAn, n_sim = n_sim)
            same_counter = 0
            while (same_counter < 3){ # when the result hasn't changed in 3 jumps, aka convergence.
              if (temp_calc[1] > alpha){
                # if our type one error is too high and we need conf_tox to be higher
                lower = temp_conf_tox # we could simply not record these?
                temp_conf_tox = (temp_conf_tox + upper)/2
                
              }
              else{
                # if our type one errr is too low and we need conf_tox to be lower
                upper = temp_conf_tox
                temp_conf_tox = (temp_conf_tox + lower)/2
                
                
              }
              
              new_calc = OC_gen(safety_supe = TRUE, eff_treat = 0.9, tox_treat = 0.5, conf_eff = temp_conf_eff, conf_tox = temp_conf_tox, N = n, IAn = IAn, n_sim = n_sim)
              if (new_calc[1] == temp_calc[1]){ ## WILL ONLY TERMINATE THIS WAY BECAUSE OF BINARY DATA
                same_counter = same_counter + 1
              }
              else{
                same_counter = 0
              }
              if (abs(new_calc[1] - alpha) < 0.0001){ #  or if it's close enough to 0.5
                same_counter = 3 # just end the loops there 
              }
              temp_calc = new_calc
              
            }
            
            temp_power = OC_gen(safety_supe = TRUE, eff_treat = 0.96, conf_eff = temp_conf_eff, conf_tox = temp_conf_tox, N = n, IAn = IAn, n_sim = n_sim)[1]
            if (temp_power > max_power){
              max_power = temp_power
              current_threshold =  c(temp_conf_eff,temp_conf_tox)
            }
            
          }
          
        }
        
        
        return(c("Power" = max_power, "Efficacy conf" = current_threshold[1:length(IAn)], "Toxicity conf" = current_threshold[(length(IAn) + 1):(2*length(IAn))]))
      }
      
      else{
        # recursive part, start at first IA and go up until we get it.
        
        for (i in 1:length(IAn)){
          
          temp_data = threshold_safety_finder(alpha = alpha*alpha_spending[i],n = IAn[i], IAn = IAn[1:i], pregiven_eff = pregiven_eff, pregiven_tox = pregiven_tox, n_sim = n_sim) # goes from bottom 
          
          pregiven_eff = temp_data[2:(i + 1)]
          pregiven_tox = temp_data[(i+2):(2*i + 1)]
          
          
          
        }
        
        return(c("Efficacy conf" = pregiven_eff, "Toxicity conf" = pregiven_tox))
      }
    }
  }
  
  
  threshold_finder_v2_no_safety = function(tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, alpha = 0.05, n = 300, IAn = c(n), alpha_spending = seq(0.1,1,length.out = (length(IAn))), pregiven_conf = numeric(0), test = 1, n_sim = 10000, NI_margin = 0.06, RAR = 0, fut_stop = FALSE, safe_stop = FALSE, safe_conf = rep(0.2,length(IAn)), norm_approx = TRUE, bisection = TRUE, iter = 11){
    ## Updates the threshold finder to call threshold finder and do it recursively that way
    
    if (length(IAn) == 1){
      alpha_spending = 1
    }
    
    conf_vec = numeric(0)
    
    for (i in 1:length(IAn)){
      
      ## Find the threshold of the first IA and then build up
      
      if (test == 1){
        
        if (bisection){
          
          # Assuming worst case scenario, toxicity not better in null
          conf_vec = threshold_finder_bisection(tox_con,eff_con,tox_con, N = IAn[i], IAn = IAn[1:i], alpha = alpha*alpha_spending[i], NI_margin = NI_margin, n_sim = n_sim, pre_conf = conf_vec, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx, iter = iter)
          
        }
        else{
          conf_vec = threshold_finder(tox_con,eff_con,tox_con, N = IAn[i], IAn = IAn[1:i], alpha = alpha*alpha_spending[i], NI_margin = NI_margin, n_sim = n_sim, pre_conf = conf_vec, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)
        }
        
      }
      else if (test == 3){
        
        
        if (bisection){
          conf_vec = threshold_finder_bayes_bisection(tox_con,eff_con,tox_con, N = IAn[i], IAn = IAn[1:i], alpha = alpha*alpha_spending[i], NI_margin = NI_margin, n_sim = n_sim, pre_conf = conf_vec, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx, iter = iter)
          
        }
        else{
          conf_vec = threshold_finder_bayes_factor(tox_con,eff_con,tox_con, N = IAn[i], IAn = IAn[1:i], alpha = alpha*alpha_spending[i], NI_margin = NI_margin, n_sim = n_sim, pre_conf = conf_vec, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)
        }
      }
      else if (test == 4){
        temp_alpha = 0.05
        
        high_alpha = 0.1
        low_alpha = 0
    
        # maximum of 20 simulations per optimisation
        for(i in 1:10){
          # bisection method, arbitarily chosen 10
          
          temp_res_alpha = OC_gen(eff_con - NI_margin, tox_treat = tox_con ,eff_con = eff_con, tox_con = tox_con, test = 4, N = n, IAn = IAn, alpha = temp_alpha, n_sim = n_sim, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1]
          
          if (temp_res_alpha > alpha + 0.001){
            # alpha too high
            
            high_alpha = temp_alpha
            temp_alpha = (low_alpha + temp_alpha) / 2
            
          }
          else if (temp_res_alpha < alpha - 0.001){
            # alpha too low
            
            low_alpha = temp_alpha
            temp_alpha = (high_alpha + temp_alpha) / 2
            
          }
          else{
            # alpha unclear, more res needed
            
            # same again but more simulations this time
            temp_res_alpha = OC_gen(eff_con - NI_margin, tox_treat = tox_con ,eff_con = eff_con, tox_con = tox_con, test = 4, N = n, IAn = IAn, alpha = temp_alpha, n_sim = n_sim*10, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1]
            
            if (temp_res_alpha > alpha){
              # alpha too high
              
              high_alpha = temp_alpha
              temp_alpha = (low_alpha + temp_alpha) / 2
              
            }
            else{
              # alpha too low
              
              low_alpha = temp_alpha
              temp_alpha = (high_alpha + temp_alpha) / 2
              
            
            }
            
          }
        }
        
        return(temp_alpha)
        
        
        # defunct increment code
        if (FALSE){
          while(TRUE){
            
            #browser()
            
            if (OC_gen(0.96 - NI_margin,0.5, test = 4, N = n, IAn = IAn, alpha = temp_alpha, n_sim = n_sim, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1] > alpha){
              
              temp_alpha = temp_alpha - 0.001
            }
            else{
              
              if(temp_alpha == 0.05){ #started too low, need to go higher
                while (OC_gen(0.96 - NI_margin, test = 4, N = n, IAn = IAn, alpha = temp_alpha + 0.001, n_sim = n_sim, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, norm_approx = norm_approx)[1] < alpha){
                  temp_alpha = temp_alpha + 0.01
                    
                }
              }
              
              return(temp_alpha)
            }
            
          }
        }
        
      }
      
      
      
    }
    
    return(conf_vec)
    
  }
  # c(0.9945024, 0.9590400) is the boundary for N = 300, IAn = 150, 300. For standard posterior
  # c(0.00563771, 0.13348389) is the boundary for "" for Bayes factor
  # c(0.9956020, 0.9790285, 0.9672320) for N = 300, IAn = 100,200,300 for standard post
  # c(0.00563771, 0.03167635, 0.10011292) for "" for Bayes factor
  bayes_factor_comp = function(a_1,b_1,a_2,b_2,c = 0.06, prior_a1 = 1, prior_b1 = 1, prior_a2 = 1, prior_b2 = 1, norm_approx = TRUE){
    # Computers the Bayes factor for the NI test
    # accurate until about the 1/million, can be improved to 1/10 million but then is noticeably slower
    
    H1_D = 1 - rcpp_exact_beta(a_1,b_1,a_2,b_2)
    
    H1 = 1 - rcpp_exact_beta(prior_a1,prior_b1,prior_a2,prior_b2)
    
    if (norm_approx){
      
      H0_D = 1 - norm_diff_cdf(a_1,b_1,a_2,b_2, c = c)
      
    }
    else{
      H0_D = 1 - beta_sim(a_1,b_1,a_2,b_2, c = c)
      
      counter = 1
      
      # If the posterior is 0, tries to get a better resolution before declaring it 0
      while(H0_D == 0){
        counter = counter + 1
        H0_D = 1 - beta_sim(a_1,b_1,a_2,b_2, c = c, n_sim = 10000*(10^counter))
        if (counter > 1){ # increase to increase the resolution by a magnitude 
          break
        }
      }
    }
    
    
    if (norm_approx){
      H0 = norm_diff_cdf(prior_a1,prior_b1,prior_a2,prior_b2, c = c)
    }
    else{
      H0 = 1 - beta_sim(prior_a1,prior_b1,prior_a2,prior_b2, c = c)
    }
    
    bayes_factor = (H0_D/H1_D)*(H1/H0) # small is evidence for the alternative
    
    #browser()
    
    if (is.na(bayes_factor)){
      browser()
    }
    
    return(bayes_factor)
  }
  
  
  if(FALSE){
    print("No IA Post")
    print(threshold_finder(0.4,0.96,0.5,N = 300, n_sim = 10000))
    
    print("No IA Bayes Factor")
    print(threshold_finder_bayes_factor(0.4,0.96,0.5,N = 300, n_sim = 10000))
    
    print("1 IA Post")
    print(threshold_finder_v2_no_safety(IAn = c(150,300)))
    
    print("1 IA Bayes")
    print(threshold_finder_v2_no_safety(IAn = c(150,300), test = 3))
    
    print("2 IA Post")
    print(threshold_finder_v2_no_safety(IAn = c(100,200,300)))
    
    print("2 IA Bayes")
    print(threshold_finder_v2_no_safety(IAn = c(100,200,300), test = 3))
    
    
    
  }
  
  
  # Look at strength of IA, how much alpha is allowed to be spent, and will
  # decide on how much alpha to be spent
  
  if (FALSE){
    
    #print("0.01")
    #print(threshold_finder_v2_no_safety(IAn = c(150,300), alpha_spending = c(0.01,1)))
    # c(0.999299, 0.948960)
    
    print("0.1")
    print(threshold_finder_v2_no_safety(IAn = c(150,300), alpha_spending = c(0.1,1)))
    # c(0.9947773, 0.9489600)
    
    print("0.2")
    print(threshold_finder_v2_no_safety(IAn = c(150,300), alpha_spending = c(0.2,1)))
    # c(0.9897995, 0.9489600)
    
    print("0.3")
    print(threshold_finder_v2_no_safety(IAn = c(150,300), alpha_spending = c(0.3,1), n_sim = 100000))
    # c(0.9840616, 0.9671040) 0.9840616 0.9671040
    
    print("0.4")
    print(threshold_finder_v2_no_safety(IAn = c(150,300), alpha_spending = c(0.4,1)))
    # c( 0.9800771, 0.9489600)
    
    print("0.5")
    print(threshold_finder_v2_no_safety(IAn = c(150,300), alpha_spending = c(0.5,1)))
    # c( 0.9750963, 0.9591680 )
    
    
  }
  
  # Want to compare the power and ESS of all of these, esp compared to the frequentist version
  
  # OC_gen(0.9, test = 4, N = 300, IAn = c(150,300)) # alpha = 0.0355 controls type one at 0.5
  
  if (FALSE){
    # Frequentist
    print(OC_gen(0.9, test = 4, N = 300, IAn = c(150,300), alpha = 0.0355, n = 100000))
    print(OC_gen(0.96, test = 4, N = 300, IAn = c(150,300), alpha = 0.0355, n = 100000))
    
    # Bayes Posterior
    # 0.01
    
    print(OC_gen(0.9, test = 1, N = 300, IAn = c(150,300), conf = c(0.999299, 0.948960)))
    print(OC_gen(0.96, test = 1, N = 300, IAn = c(150,300), conf = c(0.999299, 0.948960)))
    
    # 0.1
    
    print(OC_gen(0.9, test = 1, N = 300, IAn = c(150,300), conf = c(0.9947773, 0.9489600)))
    print(OC_gen(0.96, test = 1, N = 300, IAn = c(150,300), conf = c(0.9947773, 0.9489600)))
    
    # 0.2
    
    print(OC_gen(0.9, test = 1, N = 300, IAn = c(150,300), conf = c(0.9897995, 0.9489600)))
    print(OC_gen(0.96, test = 1, N = 300, IAn = c(150,300), conf = c(0.9897995, 0.9489600)))
    
    # 0.3
    
    print(OC_gen(0.9, test = 1, N = 300, IAn = c(150,300), conf = c(0.9840616, 0.9671040)))
    print(OC_gen(0.96, test = 1, N = 300, IAn = c(150,300), conf = c(0.9840616, 0.9671040)))
    
    # 0.4
    
    print(OC_gen(0.9, test = 1, N = 300, IAn = c(150,300), conf = c( 0.9800771, 0.9489600)))
    print(OC_gen(0.96, test = 1, N = 300, IAn = c(150,300), conf = c( 0.9800771, 0.9489600)))
    
    # 0.5
    
    print(OC_gen(0.9, test = 1, N = 300, IAn = c(150,300), conf = c( 0.9750963, 0.9591680 )))
    print(OC_gen(0.96, test = 1, N = 300, IAn = c(150,300), conf = c( 0.9750963, 0.9591680 )))
  }
  
  
  # Do some test which compares all of the three methods against different NI_margins
  # Probably going from 0.03 to 0.1 to test what maximum sample size is required for each power
  
  ss_finder_v2 = function(power = 0.8, type_1 = 0.05, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, NI_margin = 0.06, conf = c(0.9840616, 0.9671040), n_sim = 10000, safety_supe = FALSE, conf_tox = 0.99, conf_eff = 0.99, test = 1, alpha = 0.05){
    ## Finds the sample size required to get type_1 type one error with power power.
    # Will increment in 100s, and go down if required
    
    if (test == 1 | test == 3){
      # Bayesian Post
      n_start = 100
      
      while(TRUE){ # Will continue looping until break
        
        # Type one error check
        
        N = n_start
        print(N)
        
        conf = threshold_finder_v2_no_safety(n = N, IAn =  c(N/2,N), alpha_spending = c(0.3,1), NI_margin = NI_margin, test = test)
      
        
        # Power checl
        if(OC_gen(eff_con, tox_treat, eff_con, tox_con, N, c(N/2,N), NI_margin, n_sim = n_sim, conf = conf, power = TRUE, test = test) > power){
          # Go down in 20s and check which is best
          
          for (i in 1:4){
          
            
            N = n_start - 20
            print(N)
            
            if (OC_gen(eff_con, tox_treat, eff_con, tox_con, N - 20, c(N/2 - 10, N - 20), NI_margin, n_sim = n_sim, conf = conf, power = TRUE, test = test) > power){
              n_start = n_start - 20
            }
            else{
              return(c(n_start,conf))
            }
          }
          
          return(c(n_start,conf))
          
        }
        
        n_start = n_start + 100 # onto next 100
        
        
        
      }
    }
    else if (test == 4){
      # Bayesian Post
      n_start = 100
      
      while(TRUE){ # Will continue looping until break
        
        # Type one error check
        
        N = n_start
        print(N)
        
        alpha = threshold_finder_v2_no_safety(n = N, IAn =  c(N/2,N), alpha_spending = c(0.3,1), NI_margin = NI_margin, test = test)
        
        
        # Power checl
        if(OC_gen(eff_con, tox_treat, eff_con, tox_con, N, c(N/2,N), NI_margin, n_sim = n_sim, alpha = alpha, power = TRUE, test = test) > power){
          # Go down in 20s and check which is best
          
          for (i in 1:4){
            
            
            N = n_start - 20
            print(N)
            
            if (OC_gen(eff_con, tox_treat, eff_con, tox_con, N - 20, c(N/2 - 10, N - 20), NI_margin, n_sim = n_sim, alpha = alpha, power = TRUE, test = test) > power){
              n_start = n_start - 20
            }
            else{
              return(c(n_start,alpha))
            }
          }
          
          return(c(n_start,alpha))
          
        }
        
        n_start = n_start + 100 # onto next 100
        
        
        
      }
    }
    
  }
  
  ss_finder_v3 = function(power = 0.8, type_1 = 0.05, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, NI_margin = 0.06, n_sim = 10000, test = 1, alpha = 0.05, fut_stop = FALSE, safe_stop = FALSE, num_IA = 1, iter = 11){
    
    ## Will use bisection based methods to find the best sample size
    
    poss_ss = seq(100,1000, 20)
    
    low_index = 1
    high_index = length(poss_ss)
    
    cur_index = round((low_index + high_index)/2)
    
    for(i in 1:6){ # log_2(45) = 5.5, so 7 splits should be enough (90 being how many possible sample sizes there are)
      
      print(poss_ss[cur_index])
      
      temp_IAn = seq(poss_ss[cur_index]/num_IA, poss_ss[cur_index], poss_ss[cur_index]/num_IA) # definitely works for 1 or 2, will likely run into rounding errors for >=3
    
      temp_conf = threshold_finder_v2_no_safety(alpha = type_1, n = poss_ss[cur_index], IAn = temp_IAn, alpha_spending = seq(1/num_IA, 1, 1/num_IA), NI_margin = NI_margin, test = test)  
      
      #browser()
      
      temp_power = OC_gen(eff_con, tox_treat, eff_con, tox_con, poss_ss[cur_index], temp_IAn, NI_margin, n_sim = n_sim, conf = temp_conf, power = FALSE, test = test, fut_stop = fut_stop, safe_stop = safe_stop)[1]
      
      if (temp_power > power + 0.01){ # power too high, can lower ss
        
        high_index = cur_index
        cur_index = round((low_index + cur_index) / 2)
        
        
      }
      else if (temp_power < power - 0.01){ # power too low, needs to raise ss 
        low_index = cur_index
        cur_index = round((high_index + cur_index) / 2)
      }
      else{ # bang on and can end early
        return(poss_ss[cur_index])
      }
      
    }
    
    return(poss_ss[cur_index])
    
  }
  
  NI_margin_ss_plot = function(){
    # Will be a parallel function that finds the sample size required at different NI margins for the 
    # three different functions, plot them on the same plot, then return
    
    cores = detectCores()
    
    cl = makeCluster(cores/2)  ## To not get overloaded
    
    registerDoParallel(cl) 
    
    # So, want to cycle through 0.03-0.1 NI margins, across 3 different types of test,
    # so 8*3 types of tests, 24 tests, so split across cores hopefully shouldn't take too long
    
    par_df = data.frame("Sample Size" = numeric(24), "NI_margin" = rep((3:10)/100, 3),"Test" = rep(c(1,3,4), each = 8))
    
    par_ss = numeric(24)
    
    par_ss = foreach( i = 1:24, .combine = c, .export = c("ss_finder","bayes_factor_comp","OC_gen","threshold_finder_v2_no_safety","threshold_finder_bayes_factor", "threshold_finder", "rcpp_exact_beta_4_arm","rcpp_exact_beta_3_arm","rcpp_exact_beta","NI_trial_RAR","beta_sim","log_exact_beta","ss_finder_v3","mean_beta","var_beta","norm_diff_cdf","threshold_finder_bisection","threshold_finder_bayes_bisection")) %dopar%{
      
      # temp_test = par_df[i,3]
      # temp_NI = par_df[i,2]
      
      return(ss_finder_v3(NI_margin = par_df[i,2], test = par_df[i,3]))
      
    }
    
    
    stopCluster(cl)
    
    browser()
    
    par_df$Sample.Size = par_ss
    
    return(par_df)
    
    
  }
  
  #to_plot_df = NI_margin_ss_plot()
  #print(to_plot_df)
  
  #ggplot_NI = ggplot(to_plot_df, aes(x=NI_margin, y=Sample.Size)) +
  #  geom_line(aes(colour = factor(Test))) +
  #  labs( title = "Sample size required for 80% power against NI margin") +
  #  scale_color_hue(labels = c("Bayes Post", "Bayes Factor","Normal CI")) +
  #  scale_y_continuous(breaks = seq(100, 1000, by = 100))
    
    
  
  #microbenchmark(NI_trial_RAR(0.9,0.4,0.96,0.5,N = 320, IAn = c(80,160,240,320), RAR = 0))
  #microbenchmark(NI_trial_RAR(0.9,0.4,0.96,0.5,N = 320, IAn = c(80,160,240,320), RAR = 0), fut_stop = FALSE)
  #microbenchmark(NI_trial_RAR(0.9,0.4,0.96,0.5,N = 320, IAn = c(80,160,240,320), RAR = 1))
  #microbenchmark(NI_trial_RAR(0.9,0.4,0.96,0.5,N = 320, IAn = c(80,160,240,320), RAR = 1), fut_stop = FALSE)
  #microbenchmark(NI_trial_RAR(0.9,0.4,0.96,0.5,N = 320, IAn = c(80,160,240,320), RAR = 4))
  #microbenchmark(NI_trial_RAR(0.9,0.4,0.96,0.5,N = 320, IAn = c(80,160,240,320), RAR = 4), fut_stop = FALSE)
  
  if (FALSE){ # not normal approx
    print("Default")
    print(threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,360), norm_approx = FALSE))
    # c(0.9918396, 0.9790940, 0.9790940, 0.9831572)
    
    print("Futility stop")
    print(threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,360), fut_stop = TRUE, norm_approx = FALSE))
    # c(0.9926986, 0.9846384, 0.9673344, 0.9804047)
    
    print("Safety stop")
    print(threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,360), safe_stop = TRUE, norm_approx = FALSE))
    # c(0.9918396, 0.9815258, 0.9676358, 0.9573344)
    
    print("Both stop")
    print(threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,360), fut_stop = TRUE, safe_stop = TRUE, norm_approx = FALSE))
    # c(0.9918396, 0.9808481, 0.9671443, 0.9589466) # tinkered with a bit to get closer to 0.05 alpha
  }
  
  
  if (FALSE){ #  normal approx
    print("Default")
    print(threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,360)))
    # c(0.9914369, 0.9879623, 0.9831572, 0.9739466)
    
    print("Futility stop")
    print(threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,360), fut_stop = TRUE))
    # c(0.9955848, 0.9790940, 0.9790940, 0.9781572)
    
    print("Safety stop")
    print(threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,360), safe_stop = TRUE))
    # c(0.9914369, 0.9913765, 0.9789466, 0.9491680)
    
    print("Both stop")
    print(threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,360), fut_stop = TRUE, safe_stop = TRUE))
    # c(0.9955848, 0.9790940, 0.9711912, 0.9581572) # tinkered with a bit to get closer to 0.05 alpha
  }
  
  if (FALSE){
    print("100 Mono")
    temp_conf = threshold_finder(0.5,0.96,0.5, N = 100, n_sim = 10000)
    print(OC_gen(0.9, N = 100, conf = temp_conf))
    print(OC_gen(0.96, N = 100, conf = temp_conf))
    
    print("100 Bisection")
    temp_conf = threshold_finder_bisection(0.5,0.96,0.5, N = 100, n_sim = 10000)
    print(OC_gen(0.9, N = 100, conf = temp_conf))
    print(OC_gen(0.96, N = 100, conf = temp_conf))
    
    print("200 Mono")
    temp_conf = threshold_finder(0.5,0.96,0.5, N = 200, n_sim = 10000)
    print(OC_gen(0.9, N = 200, conf = temp_conf))
    print(OC_gen(0.96, N = 200, conf = temp_conf))
    
    print("200 Bisection")
    temp_conf = threshold_finder_bisection(0.5,0.96,0.5, N = 200, n_sim = 10000)
    print(OC_gen(0.9, N = 200, conf = temp_conf))
    print(OC_gen(0.96, N = 200, conf = temp_conf))
    
    print("300 Mono")
    temp_conf = threshold_finder(0.5,0.96,0.5, N = 300, n_sim = 10000)
    print(OC_gen(0.9, N = 300, conf = temp_conf))
    print(OC_gen(0.96, N = 300, conf = temp_conf))
    
    print("300 Bisection")
    temp_conf = threshold_finder_bisection(0.5,0.96,0.5, N = 300, n_sim = 10000)
    print(OC_gen(0.9, N = 300, conf = temp_conf))
    print(OC_gen(0.96, N = 300, conf = temp_conf))
    
    print("400 Mono")
    temp_conf = threshold_finder(0.5,0.96,0.5, N = 400, n_sim = 10000)
    print(OC_gen(0.9, N = 400, conf = temp_conf))
    print(OC_gen(0.96, N = 400, conf = temp_conf))
    
    print("400 Bisection")
    temp_conf = threshold_finder_bisection(0.5,0.96,0.5, N = 400, n_sim = 10000)
    print(OC_gen(0.9, N = 400, conf = temp_conf))
    print(OC_gen(0.96, N = 400, conf = temp_conf))
    
    
    
    
  }
  
  if (FALSE){
    ss_finder_v3(NI_margin = 0.03)
    ss_finder_v2(NI_margin = 0.03)
  }
  
  prop_plot = function(RAR = 1){
    # Will plot the proportion for each of the different statistics under the 4 ways of ending early
    
    par_df = data.frame("Proportion" = numeric(24), "Statistic" = rep(c(1,3,4), each = 4),"Ending_Early" = rep(c(1,2,3,4), 3))
    
    par_prop = numeric(12)
    
    cores = detectCores()
    
    cl = makeCluster(cores/2)  ## To not get overloaded
    
    registerDoParallel(cl) 
    
    #for (i in 1:12){
    par_prop = foreach( i = 1:12, .combine = c, .export = c("ss_finder","bayes_factor_comp","OC_gen","threshold_finder_v2_no_safety","threshold_finder_bayes_factor", "threshold_finder", "rcpp_exact_beta_4_arm","rcpp_exact_beta_3_arm","rcpp_exact_beta","NI_trial_RAR","beta_sim","log_exact_beta","ss_finder_v3","mean_beta","var_beta","norm_diff_cdf","threshold_finder_bisection","safer_allo","threshold_finder_bayes_bisection")) %dopar%{
      
      # temp_test = par_df[i,3]
      # temp_NI = par_df[i,2]
      #print(i)
      
      fut_stop = FALSE
      safe_stop = FALSE
      
      temp_end = par_df$Ending_Early[i]
      
      if (temp_end == 2 | temp_end == 4){
        fut_stop =TRUE
      }
      if (temp_end == 3 | temp_end == 4){
        safe_stop =TRUE
      }
  
      temp_conf = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = par_df$Statistic[i], fut_stop = fut_stop, safe_stop = safe_stop, RAR = RAR )
      
      
      temp_results = OC_gen(N = 320, IAn = c(80,160,240,320), test = par_df$Statistic[i], fut_stop = fut_stop, safe_stop = safe_stop, conf = temp_conf, RAR = RAR, alpha = temp_conf)
      
      temp_prop = (temp_results[2] + temp_results[3]) / temp_results[6]
      
      return(temp_prop)
      
      #par_prop[i] = temp_prop
    }
    
    
    stopCluster(cl)
  
    par_df$Proportion = par_prop
    
    #browser()
    
    #ggplot_prop = ggplot(par_df, aes(x=factor(Ending_Early), y=Proportion, fill = factor(Ending_Early))) +
    #  geom_bar(position = "dodge")
    
    return(par_df)
    
  }
  
  if (FALSE){
    # Optimisations
    
    conf_post_RAR0 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 0, test = 1)
    
    print(1)
    
    conf_post_RAR1 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 1, test = 1)
    
    print(2)
  
    conf_post_RAR2 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 2, test = 1)
    
    print(3)
    
    conf_post_RAR3 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 3, test = 1)
    
    print(4)
    
    conf_post_RAR4 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 4, test = 1)
    
    #conf_factor_RAR0 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 0, test = 3)
    
    print("1")
    
    #conf_factor_RAR1 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 1, test = 3)
    
    print("2")
    
    #conf_factor_RAR2 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 2, test = 3)
    
    print("3")
    
    #conf_factor_RAR3 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 3, test = 3)
    
    print("4")
    
    #conf_factor_RAR4 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 4, test = 3)
    
    conf_normal_RAR0 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 0, test = 4)
    
    print("5")
    
    conf_normal_RAR1 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 1, test = 4)
    
    print("6")
    
    conf_normal_RAR2 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 2, test = 4)
    
    print("7")
    
    conf_normal_RAR3 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 3, test = 4)
    
    print("8")
    
    conf_normal_RAR4 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 4, test = 4)
    
  } else{
    # Will fill automatically with what it was last time
    conf_post_RAR0 = c(0.9920654, 0.9802002, 0.9691650, 0.9546143)
    
    conf_post_RAR1 = c(0.9921631, 0.9799072, 0.9696045, 0.9582275)
    
    conf_post_RAR2 = c(0.9945557, 0.9783447, 0.9696533, 0.9574463)
    
    conf_post_RAR3 = c(0.9921143, 0.9801025, 0.9684326, 0.9555908)
    
    conf_post_RAR4 = c(0.9941650, 0.9799072, 0.9696045, 0.9595947)
    
    conf_factor_RAR0 = c( 0.008116943, 0.009997803, 0.101098633, 0.159106445)
    
    conf_factor_RAR1 = c( 0.005484619, 0.009997803, 0.117797852, 0.139770508)
    
    conf_factor_RAR2 = c( 0.007870850, 0.009997803, 0.122631836, 0.113842773)
    
    conf_factor_RAR3 = c( 0.005348389, 0.009997803, 0.105053711, 0.182836914)
    
    conf_factor_RAR4 = c( 0.005779053, 0.009997803, 0.106372070, 0.143725586)
    
    conf_normal_RAR0_FS = c(0.0426)
    
    conf_normal_RAR1_FS = c(0.0426)
    
    conf_normal_RAR2_FS = c(0.0426)
    
    conf_normal_RAR3_FS = c(0.0426)
    
    conf_normal_RAR4_FS = c(0.0426)
    
    #conf_list = list(conf_post_RAR0,conf_post_RAR1,conf_post_RAR2,conf_post_RAR3,conf_post_RAR4,conf_factor_RAR0,conf_factor_RAR1,conf_factor_RAR2,conf_factor_RAR3,conf_factor_RAR4,conf_normal_RAR0,conf_normal_RAR1,conf_normal_RAR2,conf_normal_RAR3,conf_normal_RAR4 )
    
  }
  
  red_green_vec = function(vec, high = TRUE){
    # Takes a vector and makes an empty string vector, except the highest is green and smallest is red
    # This is reversed if high = FALSE
    
    if (high){
      green_i = which.max(vec)
      red_i = which.min(vec)
    }
    else{
      green_i = which.min(vec)
      red_i = which.max(vec)
    }
    
    return_vec = numeric(length(vec))
    
    for (i in 1:length(vec)){
      if (i == green_i){
        return_vec[i] = "green"
      }
      else if (i == red_i){
        return_vec[i] = "red"
      }
      else{
        return_vec[i] = ""
      }
    }  
    
    return(return_vec)
    
  }
  
  table_gen = function(){
    # Makes a big table of all the data you could want
    # Will record
      # Type one error
      # Power
      # ESS | H0
      # ESS | H1
      # Prop of patients in supe treatment | H0
      # Prop of patients in supe treatment | H1
      # Number of toxics | H1
      # NUmber of deaths | H0 
    # For
      # 3 different types of endpoint
      # All types of early stopping
      # 5 different types of RAR (including ER)
    # so 15 settings in both the null and alternate hypothesis
    
    if (TRUE){
      # For ease of reference
      conf_post_RAR0 = c(0.9920654, 0.9802002, 0.9691650, 0.9546143)
      
      conf_post_RAR1 = c(0.9921631, 0.9799072, 0.9696045, 0.9582275)
      
      conf_post_RAR2 = c(0.9945557, 0.9783447, 0.9696533, 0.9574463)
      
      conf_post_RAR3 = c(0.9921143, 0.9801025, 0.9684326, 0.9555908)
      
      conf_post_RAR4 = c(0.9941650, 0.9799072, 0.9696045, 0.9595947)
      
      conf_factor_RAR0 = c( 0.008116943, 0.009997803, 0.101098633, 0.159106445)
      
      conf_factor_RAR1 = c( 0.005484619, 0.009997803, 0.117797852, 0.139770508)
      
      conf_factor_RAR2 = c( 0.007870850, 0.009997803, 0.122631836, 0.113842773)
      
      conf_factor_RAR3 = c( 0.005348389, 0.009997803, 0.105053711, 0.182836914)
      
      conf_factor_RAR4 = c( 0.005779053, 0.009997803, 0.106372070, 0.143725586)
      
      conf_normal_RAR0 = c(0.0426)
      
      conf_normal_RAR1 = c(0.0426)
      
      conf_normal_RAR2 = c(0.0426)
      
      # conf_normal_RAR3 = c(0.0426)
      
      conf_normal_RAR4 = c(0.0426)
      
      conf_list = list(conf_post_RAR0,conf_post_RAR1,conf_post_RAR2,conf_post_RAR3,conf_post_RAR4,conf_factor_RAR0,conf_factor_RAR1,conf_factor_RAR2,conf_factor_RAR3,conf_factor_RAR4,conf_normal_RAR0,conf_normal_RAR1,conf_normal_RAR2,conf_normal_RAR3,conf_normal_RAR4 )
      
    }
    
    
    cores = detectCores()
    
    cl = makeCluster(cores/2)  ## To not get overloaded
    
    registerDoParallel(cl) 
    
    
    #for (i in 1:15){
    par_return  = foreach( i = 1:15, .export = c("ss_finder","bayes_factor_comp","OC_gen","threshold_finder_v2_no_safety","threshold_finder_bayes_factor", "threshold_finder", "rcpp_exact_beta_4_arm","rcpp_exact_beta_3_arm","rcpp_exact_beta","NI_trial_RAR","beta_sim","log_exact_beta","ss_finder_v3","mean_beta","var_beta","norm_diff_cdf","threshold_finder_bisection","safer_allo","threshold_finder_bayes_bisection")) %dopar%{
      
      
      temp_RAR = (i-1)%%5 
      
      temp_test_vec = c(1,3,4)
      temp_test = (i-1)%/%5 + 1
      
      temp_test = temp_test_vec[temp_test]
      
      temp_conf = conf_list[[i]]
      temp_alpha = temp_conf
      
      #print(i)
      
      OC_H0 = OC_gen(0.9,0.5,0.96,0.5, N = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = temp_RAR, test = temp_test, conf = temp_conf, alpha = temp_alpha, n_sim = 100000)
      
      #print(-i)
      
      OC_H1 = OC_gen(0.96,0.4,0.96,0.5, N = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = temp_RAR, test = temp_test, conf = temp_conf, alpha = temp_alpha, n_sim = 100000)
      
      Type_1 = OC_H0[1]
      Power = OC_H1[1]
      ESS_H0 = OC_H0[6]
      ESS_H1 = OC_H1[6]
      
      prop_supe_H0 = (OC_H0[4] + OC_H0[5])/ESS_H0
      prop_supe_H1 = (OC_H1[2] + OC_H1[3])/ESS_H1
      
      toxic_H1 = OC_H1[8]
      deaths_H0 = OC_H0[6] - OC_H0[7]
      
      par_return = c(Type_1, Power, ESS_H0, ESS_H1, prop_supe_H0, prop_supe_H1, toxic_H1, deaths_H0)
      
      return(par_return)
      
      
    }
    
    stopCluster(cl)
    
    print("Done")
    
    table_df = data.frame("Type_I_Error" = numeric(15),"Power" = numeric(15),"ESS_H0" = numeric(15),"ESS_H1" = numeric(15),"prop_supe_H0" = numeric(15),"prop_supe_H1" = numeric(15),"toxic_H1" = numeric(15),"deaths_H0" = numeric(15), "RAR" = numeric(15), "test" = numeric(15))
    
    for (i in 1:15){
      temp_vec = par_return[[i]]
      
      table_df$Type_I_Error[i] = temp_vec[1]
      table_df$Power[i] = temp_vec[2]
      table_df$ESS_H0[i] = temp_vec[3]
      table_df$ESS_H1[i] = temp_vec[4]
      
      table_df$prop_supe_H0[i] = temp_vec[5]
      table_df$prop_supe_H1[i] = temp_vec[6]
      
      table_df$toxic_H1[i] = temp_vec[7]
      table_df$deaths_H0[i] = temp_vec[8]
      
      table_df$RAR[i] =  (i-1)%%5
      
      temp_test_vec = c(1,3,4)
      temp_test = (i-1)%/%5 + 1
      
      temp_test = temp_test_vec[temp_test]
      table_df$test[i] = temp_test
      
      
      
    }
    
    
    #kable(table_df, "pipe")
    
    #browser()
    
    #kable(table_df, col.names = gsub("[.]", " ", names(table_df)))
    
    to_return_table = table_df %>%
      kbl() %>%
      kable_styling() %>%
      kable_material_dark(c("striped", "hover"))%>%
      column_spec(1, background = red_green_vec(table_df$Type_I_Error, FALSE)) %>%
      column_spec(2, background = red_green_vec(table_df$Power, TRUE)) %>%
      column_spec(3, background = red_green_vec(table_df$ESS_H0, FALSE)) %>%
      column_spec(4, background = red_green_vec(table_df$ESS_H1, FALSE)) %>%
      column_spec(5, background = red_green_vec(table_df$prop_supe_H0, TRUE)) %>%
      column_spec(6, background = red_green_vec(table_df$prop_supe_H1, TRUE)) %>%
      column_spec(7, background = red_green_vec(table_df$toxic_H1, FALSE)) %>%
      column_spec(8, background = red_green_vec(table_df$deaths_H0, FALSE)) 
  
    print(to_return_table)
    return(to_return_table)
    
    
    
  }
  
  if (FALSE){
    
    
    #### Only Efficacy Stopping
    conf_normal_RAR0 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = FALSE, safe_stop = FALSE, RAR = 0, test = 4)
    
    print("1")
    
    conf_normal_RAR1 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = FALSE, safe_stop = FALSE, RAR = 1, test = 4)
    
    print("2")
    
    conf_normal_RAR2 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = FALSE, safe_stop = FALSE, RAR = 2, test = 4)
    
    print("3")
    
    conf_normal_RAR3 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = FALSE, safe_stop = FALSE, RAR = 3, test = 4)
    
    print("4")
    
    conf_normal_RAR4 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = FALSE, safe_stop = FALSE, RAR = 4, test = 4)
    
    
    #### Efficacy + Futility
    conf_normal_RAR0_F = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = FALSE, RAR = 0, test = 4)
    
    print("6")
    
    conf_normal_RAR1_F = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = FALSE, RAR = 1, test = 4)
    
    print("7")
    
    conf_normal_RAR2_F = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = FALSE, RAR = 2, test = 4)
    
    print("8")
    
    conf_normal_RAR3_F = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = FALSE, RAR = 3, test = 4)
    
    print("9")
    
    conf_normal_RAR4_F = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = FALSE, RAR = 4, test = 4)
    
    
    
    #### Efficacy + Safety 
    conf_normal_RAR0_S = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = FALSE, safe_stop = TRUE, RAR = 0, test = 4)
    
    print("11")
    
    conf_normal_RAR1_S = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = FALSE, safe_stop = TRUE, RAR = 1, test = 4)
    
    print("12")
    
    conf_normal_RAR2_S = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = FALSE, safe_stop = TRUE, RAR = 2, test = 4)
    
    print("13")
    
    conf_normal_RAR3_S = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = FALSE, safe_stop = TRUE, RAR = 3, test = 4)
    
    print("14")
    
    conf_normal_RAR4_S = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = FALSE, safe_stop = TRUE, RAR = 4, test = 4)
    
    
    
    #### Efficacy , safety, and futility
    conf_normal_RAR0_FS = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 0, test = 4)
    
    print("16")
    
    conf_normal_RAR1_FS = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 1, test = 4)
    
    print("17")
    
    conf_normal_RAR2_FS = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 2, test = 4)
    
    print("18")
    
    conf_normal_RAR3_FS = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 3, test = 4)
    
    print("19")
    
    conf_normal_RAR4_FS = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), fut_stop = TRUE, safe_stop = TRUE, RAR = 4, test = 4)
    
    
  } else{
    
    #### Only Efficacy Stopping
    conf_normal_RAR0 =  0.03061523
    
  
    
    conf_normal_RAR1 = 0.02885742
    
  
    
    conf_normal_RAR2 = 0.02905273
    
  
    
    conf_normal_RAR3 =  0.02924805
    
    conf_normal_RAR4 = 0.02963867
    
    
    #### Efficacy + Futility
    conf_normal_RAR0_F = 0.03168945
    
  
    
    conf_normal_RAR1_F = 0.03188477
    
    
    conf_normal_RAR2_F = 0.03061523
    
  
    
    conf_normal_RAR3_F = 0.03149414
  
    
    conf_normal_RAR4_F = 0.03149414
    
    
    
    #### Efficacy + Safety 
    conf_normal_RAR0_S = 0.0425293
    
    conf_normal_RAR1_S = 0.04243164
  
    
    conf_normal_RAR2_S = 0.04243164
    
    conf_normal_RAR3_S = 0.04262695
    
    conf_normal_RAR4_S = 0.04233398
    
    #### Efficacy , safety, and futility
    conf_normal_RAR0_FS = 0.04262695
    
    conf_normal_RAR1_FS = 0.04262695
  
    conf_normal_RAR2_FS =  0.04379883
  
    conf_normal_RAR3_FS = 0.04282227
    
    conf_normal_RAR4_FS = 0.04262695
  }
  
  normal_table_gen = function(){
    # Makes a big table of all the data you could want, this time focusing on the normal endpoint
    # Will record
    # Type one error
    # Power
    # ESS | H0
    # ESS | H1
    # Prop of patients in supe treatment | H0
    # Prop of patients in supe treatment | H1
    # Number of toxics | H1
    # NUmber of deaths | H0 
    # For
    # 3 different types of endpoint
    # All types of early stopping
    # 5 different types of RAR (including ER)
    # so 15 settings in both the null and alternate hypothesis
    
    if (TRUE){
      # For ease of reference
      #### Only Efficacy Stopping
      conf_normal_RAR0 =  0.03061523
      
      
      
      conf_normal_RAR1 = 0.02885742
      
      
      
      conf_normal_RAR2 = 0.02905273
      
      
      
      conf_normal_RAR3 =  0.02924805
      
      conf_normal_RAR4 = 0.02963867
      
      
      #### Efficacy + Futility
      conf_normal_RAR0_F = 0.03168945
      
      
      
      conf_normal_RAR1_F = 0.03188477
      
      
      conf_normal_RAR2_F = 0.03061523
      
      
      
      conf_normal_RAR3_F = 0.03149414
      
      
      conf_normal_RAR4_F = 0.03149414
      
      
      
      #### Efficacy + Safety 
      conf_normal_RAR0_S = 0.0425293
      
      conf_normal_RAR1_S = 0.04243164
      
      
      conf_normal_RAR2_S = 0.04243164
      
      conf_normal_RAR3_S = 0.04262695
      
      conf_normal_RAR4_S = 0.04233398
      
      #### Efficacy , safety, and futility
      conf_normal_RAR0_FS = 0.04262695
      
      conf_normal_RAR1_FS = 0.04262695
      
      conf_normal_RAR2_FS =  0.04379883
      
      conf_normal_RAR3_FS = 0.04282227
      
      conf_normal_RAR4_FS = 0.04262695
      
      conf_list = list(conf_normal_RAR0,conf_normal_RAR1,conf_normal_RAR2,conf_normal_RAR3,conf_normal_RAR4,conf_normal_RAR0_F,conf_normal_RAR1_F,conf_normal_RAR2_F,conf_normal_RAR3_F,conf_normal_RAR4_F,conf_normal_RAR0_S,conf_normal_RAR1_S,conf_normal_RAR2_S,conf_normal_RAR3_S,conf_normal_RAR4_S,conf_normal_RAR0_FS,conf_normal_RAR1_FS,conf_normal_RAR2_FS,conf_normal_RAR3_FS,conf_normal_RAR4_FS )
      
    }
    
    fut_vec = rep(c(0,1,0,1), each = 5)
    safe_vec = rep(c(0,0,1,1), each = 5)
    RAR_vec = rep(c(0,1,2,3,4), 4)
    
    cores = detectCores()
    
    cl = makeCluster(cores/2)  ## To not get overloaded
    
    registerDoParallel(cl) 
    
   
    
    
    #for (i in 1:20){
    par_return  = foreach( i = 1:20, .export = c("ss_finder","bayes_factor_comp","OC_gen","threshold_finder_v2_no_safety","threshold_finder_bayes_factor", "threshold_finder", "rcpp_exact_beta_4_arm","rcpp_exact_beta_3_arm","rcpp_exact_beta","NI_trial_RAR","beta_sim","log_exact_beta","ss_finder_v3","mean_beta","var_beta","norm_diff_cdf","threshold_finder_bisection","safer_allo","threshold_finder_bayes_bisection")) %dopar%{
      
      
      temp_RAR = RAR_vec[i]
      
      temp_fut = fut_vec[i]
      temp_safe = safe_vec[i]
      
      temp_alpha = conf_list[[i]]
      
      #print(i)
      
      OC_H0 = OC_gen(0.9,0.5,0.96,0.5, N = 320, IAn = c(80,160,240,320), fut_stop = temp_fut, safe_stop = temp_safe, RAR = temp_RAR, test = 4, alpha = temp_alpha, n_sim = 100000)
      
      #print(-i)
      
      OC_H1 = OC_gen(0.96,0.4,0.96,0.5, N = 320, IAn = c(80,160,240,320), fut_stop = temp_fut, safe_stop = temp_safe, RAR = temp_RAR, test = 4, conf = temp_conf, alpha = temp_alpha, n_sim = 100000)
      
      Type_1 = OC_H0[1]
      Power = OC_H1[1]
      ESS_H0 = OC_H0[6]
      ESS_H1 = OC_H1[6]
      
      prop_supe_H0 = (OC_H0[4] + OC_H0[5])/ESS_H0
      prop_supe_H1 = (OC_H1[2] + OC_H1[3])/ESS_H1
      
      toxic_H1 = OC_H1[8]
      deaths_H0 = OC_H0[6] - OC_H0[7]
      
      par_return = c(Type_1, Power, ESS_H0, ESS_H1, prop_supe_H0, prop_supe_H1, toxic_H1, deaths_H0)
      
      return(par_return)
      
      
    }
    
    stopCluster(cl)
    
    print("Done")
    
    table_df = data.frame("Type_I_Error" = numeric(20),"Power" = numeric(20),"ESS_H0" = numeric(20),"ESS_H1" = numeric(20),"prop_supe_H0" = numeric(20),"prop_supe_H1" = numeric(20),"toxic_H1" = numeric(20),"deaths_H0" = numeric(20), "RAR" = numeric(20), "Fut_Stop" = numeric(20), "Safe_Stop" = numeric(20))
    
    for (i in 1:20){
      temp_vec = par_return[[i]]
      
      table_df$Type_I_Error[i] = temp_vec[1]
      table_df$Power[i] = temp_vec[2]
      table_df$ESS_H0[i] = temp_vec[3]
      table_df$ESS_H1[i] = temp_vec[4]
      
      table_df$prop_supe_H0[i] = temp_vec[5]
      table_df$prop_supe_H1[i] = temp_vec[6]
      
      table_df$toxic_H1[i] = temp_vec[7]
      table_df$deaths_H0[i] = temp_vec[8]
      
      
      
      
    }
    
    table_df$RAR = RAR_vec
    table_df$Fut_Stop = fut_vec
    table_df$Safe_Stop = safe_vec
    
    
    
    #kable(table_df, "pipe")
    
    browser()
    
    #kable(table_df, col.names = gsub("[.]", " ", names(table_df)))
    
    to_return_table = table_df %>%
      kbl() %>%
      kable_styling() %>%
      kable_material_dark(c("striped", "hover"))%>%
      column_spec(1, background = red_green_vec(table_df$Type_I_Error, FALSE)) %>%
      column_spec(2, background = red_green_vec(table_df$Power, TRUE)) %>%
      column_spec(3, background = red_green_vec(table_df$ESS_H0, FALSE)) %>%
      column_spec(4, background = red_green_vec(table_df$ESS_H1, FALSE)) %>%
      column_spec(5, background = red_green_vec(table_df$prop_supe_H0, TRUE)) %>%
      column_spec(6, background = red_green_vec(table_df$prop_supe_H1, TRUE)) %>%
      column_spec(7, background = red_green_vec(table_df$toxic_H1, FALSE)) %>%
      column_spec(8, background = red_green_vec(table_df$deaths_H0, FALSE)) 
    
    print(to_return_table)
    return(to_return_table)
    
    
    
  }
  
  
  if (FALSE){
  
    prop_df_RAR0 = prop_plot(0)
    prop_df_RAR1 = prop_plot(1)
    prop_df_RAR2 = prop_plot(2)
    prop_df_RAR3 = prop_plot(3)
    prop_df_RAR4 = prop_plot(4)
    
    ggplot_prop0 = ggplot(prop_df_RAR0, aes(x=factor(Ending_Early), y=Proportion, fill = factor(Statistic))) +
      geom_bar(stat="identity", position=position_dodge()) +
      labs( title = "ER Proportion against different types of ending early") +
      coord_cartesian(ylim = c(0.45,0.55)) +
      scale_y_continuous(breaks = seq(0.45, 0.55, by = 0.025))
      
    print(ggplot_prop0)
    
    ggplot_prop1 = ggplot(prop_df_RAR1, aes(x=factor(Ending_Early), y=Proportion, fill = factor(Statistic))) +
      geom_bar(stat="identity", position=position_dodge()) +
      labs( title = "RITS Proportion against different types of ending early") +
      coord_cartesian(ylim = c(0.45,0.55)) +
      scale_y_continuous(breaks = seq(0.45, 0.55, by = 0.025))
    
    print(ggplot_prop1)
    
    ggplot_prop2 = ggplot(prop_df_RAR2, aes(x=factor(Ending_Early), y=Proportion, fill = factor(Statistic))) +
      geom_bar(stat="identity", position=position_dodge()) +
      labs( title = "Gaussian RITS Proportion against different types of ending early") +
      coord_cartesian(ylim = c(0.45,0.55)) +
      scale_y_continuous(breaks = seq(0.45, 0.55, by = 0.025))
    
    print(ggplot_prop2)
    
    ggplot_prop3 = ggplot(prop_df_RAR3, aes(x=factor(Ending_Early), y=Proportion, fill = factor(Statistic))) +
      geom_bar(stat="identity", position=position_dodge()) +
      labs( title = "SAFER Proportion against different types of ending early") +
      coord_cartesian(ylim = c(0.45,0.55)) +
      scale_y_continuous(breaks = seq(0.45, 0.55, by = 0.025))
    
    print(ggplot_prop3)
    
    ggplot_prop4 = ggplot(prop_df_RAR4, aes(x=factor(Ending_Early), y=Proportion, fill = factor(Statistic))) +
      geom_bar(stat="identity", position=position_dodge()) +
      labs( title = "Efficacy BRAR Proportion against different types of ending early") +
      coord_cartesian(ylim = c(0.45,0.55)) +
      scale_y_continuous(breaks = seq(0.45, 0.55, by = 0.025))
    
    print(ggplot_prop4)
  
  }
  
  T1 = Sys.time()
  #threshold_finder_v2_no_safety(n = 320, bisection = FALSE, norm_approx = FALSE)
  T1 = Sys.time() - T1
  
  T2 = Sys.time()
  
  #threshold_finder_v2_no_safety(n = 320)
  T2 = Sys.time() - T2
  
  print(T1)
  print(T2)
  
  
  
  ##### CTS Section ######
  
  
  NI_trial_RAR_cts = function(eff_treat = 0.96, tox_treat = 0.05, eff_con = 0.96, tox_con = 0.15, censor_time = 3, eff_treat_rate = 0, tox_treat_rate = 0, eff_con_rate = 0, tox_con_rate = 0,  N, IAn = N, NI_margin = 0.06, conf = c(0.99,0.98,0.97,0.96), test = 1, alpha = 0.05, RAR = 0, w = 0.5, sigma = 0.1, fut_stop = TRUE, safe_stop = TRUE, safe_conf = 0.75, safe_conf_2 =0.45, safe_conf_3 = 0.7, safe_conf_4 = 0.6, safe_test = 1, norm_approx = TRUE, PH_true = TRUE){
    ## Will compare the treatment and control in an NI trial context
    # eff_treat and tox_treat are the treatments efficacy and toxicity
    # eff_con and tox_con are the control's treatment and efficacy
    # N is maximum sample size
    # IAn is a vector of interims placements (choose divisible by 2 for ease)
    # NI margin is the level of clinically acceptable reduction in efficacy
    # conf is a vector (matching the length of IAn) of thresholds to compare the posterior probabilities against
    # beta 1,1 prior
    # w is the weighting for RITS. w close to 1 faovurs the efficacy, 0 favours the safety
    # sigma tunes gaussian RITS
    # fut_stop = TRUE means there is futility stopping (stochastic Curtailment approach)
    # PH_true determines whether PH is true, and whether we need to be concerned about changing safety
    
    # test = 1 is standard NI posterior test, adjusting for the NI margin
    # test = 2 is posterior that doesn't adjust for the NI
    # test = 3 is the Bayes factor
    # test = 4 is normal approximation, alpha is only used here to calculate the width of the CI
    
    # RAR = 0 is just equal randomisation
    # RAR = 1 is RITS
    # RAR = 2 is Gaussian RITS
    # RAR = 3 is safety
    # RAR = 4 is just standard efficacy based BRAR
    
    ## NOTE: Technically if we want our probabilities of x for safety we need to 
    # account for the fact it's censored by death, yielding to:
    
    #5% safety under 4% death: exp(0.0175) 
    
    #15% safety under 4% death: exp(0.0553) 
    
    #5% safety under 10% death: exp(0.0180) 
    
    #15% safety under 10% death: exp(0.0572) 
    
    
    ## First we want to derive the rates if they're not given:
    
    if (eff_treat_rate == 0){
      eff_treat_rate = -(1/censor_time)*log(eff_treat) # means death occurs with prob 1 eff_treat
    }
    if (eff_con_rate == 0){
      eff_con_rate = -(1/censor_time)*log(eff_con) # means death occurs with prob 1 - eff_con
    }
    if (tox_treat_rate == 0){
      tox_treat_rate = -(1/censor_time)*log(1 - tox_treat) # means toxic occurs with prob tox_treat
      if (PH_true == FALSE){
        tox_treat_rate_vec = seq(0.7*tox_treat_rate, 1.3*tox_treat_rate,length.out = 4) # safety improving over time
        tox_treat_rate_vec = sort(tox_treat_rate_vec, decreasing = TRUE) # high rate high toxicity, so start with highest
      }
    }
    if (tox_con_rate == 0){
      tox_con_rate = -(1/censor_time)*log(1 - tox_con) # means toxic occurs with prob tox_con
    }
    
    
    
    
    allo_prob = 0.5 # starts at fixed, just defining here to avoid errors later
    
    eff_treat_vec = numeric(0)
    tox_treat_vec = numeric(0)
    eff_con_vec = numeric(0)
    tox_con_vec = numeric(0)
    
    safe_censor_treat = numeric(0)
    safe_compete_treat = numeric(0)
    
    safe_censor_con = numeric(0)
    safe_compete_con = numeric(0)
    
    
    
    for (i in 1:length(IAn)){  
      if (i == 1){ # how many to simulate
        
        to_sim = IAn[1]/2 # burn in, no RAR
        to_sim_con = to_sim
        to_sim_treat = to_sim
        
      }
      else{
        # Gotta update the RAR so that's it's functional with cts safety
        
        
        if (RAR == 0){
          # no RAR
          allo_prob = 0.5
        }
        else if (RAR == 1){
          # RITS
          tox_treat_1 = sum(tox_treat_vec)
          tox_treat_0 = length(tox_treat_vec) - tox_treat_1 
          
          tox_con_1 = sum(tox_con_vec)
          tox_con_0 = length(tox_con_vec) - tox_con_1
          
          
          allo_prob_eff = 1 - rcpp_exact_beta(1 + eff_treat_suc,1 + eff_treat_fail,1 + eff_con_suc,1 + eff_con_fail)
          # high favours treatment
          
          #allo_prob_tox = rcpp_exact_beta(1 + tox_treat_1, 1 + tox_treat_0, 1 + tox_con_1, 1 + tox_con_0)
          
          # In the cts case, we have the probability that an IG is bigger than another IG, which we already use in safe_test = 1
           
          
          treat_event_time = sum(tox_treat_vec) # time at risk for the treatment
          treat_num_event = sum(1 - safe_censor_treat) # number of events occurring
          
          con_event_time = sum(tox_con_vec) # time at risk for the treatment
          con_num_event = sum(1 - safe_censor_con) # number of events occurring
          
          allo_prob_tox = 1 - Bayes_cts_test(treat_event_time, treat_num_event, con_event_time, con_num_event)
          # low means favour the control
          
          # CTS Above
          
          allo_prob = w*allo_prob_eff + (1 - w)*allo_prob_tox
          
          allo_prob = allo_prob^(IAn[i-1]/(2*N))/(allo_prob^(IAn[i-1]/(2*N)) + (1-allo_prob)^(IAn[i-1]/(2*N))) # tuning with n/2N
          
          #browser()
          
        }
        else if (RAR == 2){
          # Gaussian RITS
          
          tox_treat_1 = sum(tox_treat_vec)
          tox_treat_0 = length(tox_treat_vec) - tox_treat_1 
          
          tox_con_1 = sum(tox_con_vec)
          tox_con_0 = length(tox_con_vec) - tox_con_1
          
          
          allo_prob_eff = 1 - rcpp_exact_beta(1 + eff_treat_suc,1 + eff_treat_fail,1 + eff_con_suc,1 + eff_con_fail)
          #allo_prob_tox = rcpp_exact_beta(1 + tox_treat_1, 1 + tox_treat_0, 1 + tox_con_1, 1 + tox_con_0)
          
          # In the cts case, we have the probability that an IG is bigger than another IG, which we already use in safe_test = 1
          
          
          treat_event_time = sum(tox_treat_vec) # time at risk for the treatment
          treat_num_event = sum(1 - safe_censor_treat) # number of events occurring
          
          con_event_time = sum(tox_con_vec) # time at risk for the treatment
          con_num_event = sum(1 - safe_censor_con) # number of events occurring
          
          allo_prob_tox = 1 - Bayes_cts_test(treat_event_time, treat_num_event, con_event_time, con_num_event)
          # low means favour the control
          
          
          g_x = exp( - ((allo_prob_eff - 0.5)^2) / sigma )
          
          allo_prob = (1 - g_x)*allo_prob_eff + g_x*allo_prob_tox
          
          allo_prob = allo_prob^(IAn[i-1]/(2*N))/(allo_prob^(IAn[i-1]/(2*N)) + (1-allo_prob)^(IAn[i-1]/(2*N))) # tuning with n/2N
          
          
        }
        else if (RAR == 3){
          # SAFER
          tox_treat_1 = sum(tox_treat_vec)
          tox_treat_0 = length(tox_treat_vec) - tox_treat_1 
          
          tox_con_1 = sum(tox_con_vec)
          tox_con_0 = length(tox_con_vec) - tox_con_1
          
          #allo_prob_tox = rcpp_exact_beta(1 + tox_treat_1, 1 + tox_treat_0, 1 + tox_con_1, 1 + tox_con_0)
          
          # In the cts case, we have the probability that an IG is bigger than another IG, which we already use in safe_test = 1
          
          
          treat_event_time = sum(tox_treat_vec) # time at risk for the treatment
          treat_num_event = sum(1 - safe_censor_treat) # number of events occurring
          
          con_event_time = sum(tox_con_vec) # time at risk for the treatment
          con_num_event = sum(1 - safe_censor_con) # number of events occurring
          
          allo_prob_tox = 1 - Bayes_cts_test(treat_event_time, treat_num_event, con_event_time, con_num_event)
          # low means favour the control
          
          allo_prob = safer_allo(eff_treat_vec,eff_con_vec, allo_prob_tox) 
          # will try without tuning, not expecting SAFER to be hugely relevant for now
          
        }
        else if (RAR == 4){
          # BRAR just based on efficacy
          allo_prob = 1 - rcpp_exact_beta(1 + eff_treat_suc,1 + eff_treat_fail,1 + eff_con_suc,1 + eff_con_fail)
          
          allo_prob = allo_prob^(IAn[i-1]/(2*N))/(allo_prob^(IAn[i-1]/(2*N)) + (1-allo_prob)^(IAn[i-1]/(2*N))) # tuning with n/2N
        }
        
        to_sim = (IAn[i] - IAn[i-1])
        
        to_sim_con = round(to_sim*(1 - allo_prob))
        to_sim_treat = to_sim -  to_sim_con
        
        
      }
      
      
      
      
      # simulating
      #print(eff_treat)
      if(is.na(to_sim_treat)){
        browser()
      } # error catching
      
      if (to_sim_treat < 0){
        browser()
      }    
      # Censoring time, censor_time. If the death/safety time exceeds this, then it is censored
      # and an event doesn't happen
      
      # otherwise, H0: is given by control ~ exp() efficacy, exp() safety; treat ~ exp() efficacy, exp() safety
      # H1: control ~ exp() efficacy, exp() safety; treat ~ exp() efficacy, exp() safety
      
      ### 
      # Treatment
      
      death_times_treat = rexp(to_sim_treat, eff_treat_rate)
      
      if (PH_true == FALSE){ # if we're violating PH with a changing safety
        toxic_times_treat = rexp(to_sim_treat, tox_treat_rate_vec[i])
      }
      else{
        toxic_times_treat = rexp(to_sim_treat, tox_treat_rate)
      }
      
      
      
      toxic_events_treat = (toxic_times_treat < death_times_treat) & (toxic_times_treat < 3) # safety event
      toxic_time_censor_treat = 1 - toxic_events_treat # base censoring, 1 means not a safety event
      
      toxic_compete_censor_treat = (death_times_treat < toxic_times_treat) & (death_times_treat < 3) # death happens first, so toxics are censored
      toxic_times_treat[which(toxic_compete_censor_treat)] = death_times_treat[which(toxic_compete_censor_treat)]
      toxic_compete_censor_treat = toxic_compete_censor_treat*2 #(So that censoring due to death becomes 2)
      
      
      
      death_treat_vec = 1 - (death_times_treat < censor_time)*1 #(*1 so it's numeric, for easier calcs ), 1 denotes event here
      
      #toxic_time_censor_treat = (toxic_times_treat > 3) # *1 later for numerics, but we want boolean now for which,  1 denotes censored here
      toxic_times_treat[which(toxic_times_treat > 3)] = 3 # equal to censored time
      
      # recording
      
      safe_censor_treat = append(safe_censor_treat, toxic_time_censor_treat)
      safe_compete_treat = append(safe_compete_treat, toxic_compete_censor_treat) # has a 2 whenever it's censored due to death, will be used for AJE 
      eff_treat_vec = append(eff_treat_vec, death_treat_vec)  
      tox_treat_vec = append(tox_treat_vec, toxic_times_treat)
      
      ###
      # Control
      
      
      death_times_con = rexp(to_sim_con, eff_con_rate)
      toxic_times_con = rexp(to_sim_con, tox_con_rate)
      
      toxic_events_con = (toxic_times_con < death_times_con) & (toxic_times_con < 3) # safety event
      toxic_time_censor_con = 1 - toxic_events_con # base censoring, 1 means not a safety event
      
      toxic_compete_censor_con = (death_times_con < toxic_times_con) & (death_times_con < 3) # death happens first, so toxics are censored
      toxic_times_con[which(toxic_compete_censor_con)] = death_times_con[which(toxic_compete_censor_con)]
      toxic_compete_censor_con = toxic_compete_censor_con*2 #(So that censoring due to death becomes 2)
      
      death_con_vec = 1 - (death_times_con < censor_time)*1 #(*1 so it's numeric, for easier calcs ), 1 denotes event here
      
      
      #toxic_time_censor_con = (toxic_times_con > 3) # once again *1 for numerics, 1 denotes censored here
      toxic_times_con[which(toxic_times_con > 3)] = 3 # equal to censored time
      
      
      # recording
      
      safe_censor_con = append(safe_censor_con, toxic_time_censor_con)
      safe_compete_con = append(safe_compete_con, toxic_compete_censor_con) # has a 2 whenever it's censored due to death, will be used for AJE 
      eff_con_vec = append(eff_con_vec, death_con_vec)
      tox_con_vec = append(tox_con_vec, toxic_times_con )
      
      # browser() # testing to see if the generation works as expected
      
      
      # Doing these to try and avoid appending to see if it's faster
      #new_eff_treat_vec = rbinom(to_sim_treat,1,eff_treat)
      #eff_treat_vec = c(eff_treat_vec, new_eff_treat_vec) 
      
      #new_tox_treat_vec = rbinom(to_sim_treat,1,tox_treat)
      #tox_treat_vec = c(tox_treat_vec, new_tox_treat_vec) 
      
      #new_eff_con_vec = rbinom(to_sim_con,1,eff_con)
      #eff_con_vec = c(eff_con_vec, new_eff_con_vec) 
      
      #new_tox_con_vec = rbinom(to_sim_con,1,tox_con)
      #tox_con_vec = c(tox_con_vec, new_tox_con_vec) 
      # not quicker, just going back
      
      eff_treat_suc = sum(eff_treat_vec) # successes in the treatment
      eff_treat_fail = length(eff_treat_vec) - eff_treat_suc # failures in the treatment
      
      eff_con_suc = sum(eff_con_vec) # successes in the treatment
      eff_con_fail = length(eff_con_vec) - eff_con_suc # failures in the treatment
      
      #browser()
      
      # Data ready for returning, so don't have to repeat it so often
      
      # Number of toxic events in general 
      num_toxics = sum(1 - safe_censor_con) + sum(1 - safe_censor_treat)
      num_toxics_con = sum(1 - safe_censor_con) 
      num_toxics_treat = sum(1 - safe_censor_treat)
      
      # Average time until a toxic event in the control if it occurs
      mean_toxic_time_con = mean(tox_con_vec[which(safe_censor_con == 0)])
      
      # Average time until a toxic event in the treatment if it occurs
      mean_toxic_time_treat = mean(tox_treat_vec[which(safe_censor_treat == 0)])
      
      
      # analysis
      if (test == 1){
        if (NI_margin == 0){
          # do rcpp beta
          
          post_prob = rcpp_exact_beta(1 + eff_treat_suc, 1+eff_treat_fail, 1 + eff_con_suc, 1 + eff_con_fail) # 1 + because Beta(1,1) prior
          
        }
        else{
          # do the sim beta
          if (norm_approx){
            post_prob = norm_diff_cdf(1 + eff_treat_suc, 1+eff_treat_fail, 1 + eff_con_suc, 1 + eff_con_fail , c = NI_margin )
          }
          else{
            #browser()
            post_prob = beta_sim(1 + eff_treat_suc, 1+eff_treat_fail, 1 + eff_con_suc, 1 + eff_con_fail , c = NI_margin ) # 1 + because Beta(1,1) prior
          }
        }
        
        
        if (post_prob > conf[i]){ # if we pass the decision boundary
           
  
          
          #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
          return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
          #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
          
          
        } 
        else if (fut_stop){
          eff_treat_suc_curtail = eff_treat_suc + round(eff_con*(N - IAn[i])/2) # if it does very well
          eff_con_suc_curtail = eff_con_suc + round(eff_con*(N - IAn[i])/2) # if it does as expected
          
          eff_treat_fail_curtail = eff_treat_fail + round((1 - eff_con)*(N - IAn[i])/2) # if it does very well
          eff_con_fail_curtail = eff_con_fail + round((1 - eff_con)*(N - IAn[i])/2) # if it does as expected
          
          if (norm_approx){
            post_prob_curtail = norm_diff_cdf(1 + eff_treat_suc_curtail, 1+eff_treat_fail_curtail, 1 + eff_con_suc_curtail, 1 + eff_con_fail_curtail , c = NI_margin ) # 1 + because Beta(1,1) prior
          }
          else{
            
            post_prob_curtail = beta_sim(1 + eff_treat_suc_curtail, 1+eff_treat_fail_curtail, 1 + eff_con_suc_curtail, 1 + eff_con_fail_curtail , c = NI_margin ) # 1 + because Beta(1,1) prior
            
          }
          # If the probability of it rejecting the null at the end is unlikely, then will just end the trial early
          if (post_prob_curtail < conf[length(IAn)] - 0.1){
            #browser()
            
            
            
            
            #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
            #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
            
            #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            
          }
          
        }
        if (safe_stop){
          
          if (safe_test == 1){ # BOP2 esque Bayesian Posterior testing 
            
            
            treat_event_time = sum(tox_treat_vec) # time at risk for the treatment
            treat_num_event = sum(1 - safe_censor_treat) # number of events occurring
            
            con_event_time = sum(tox_con_vec) # time at risk for the treatment
            con_num_event = sum(1 - safe_censor_con) # number of events occurring
            
            safe_prob = Bayes_cts_test(treat_event_time, treat_num_event, con_event_time, con_num_event)
            
            #browser()
            
            if (safe_prob > safe_conf){ # wth high prob is the treat more toxic
              
  
              
              
              #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
              #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
              
              #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            }
            
          }
          
          else if (safe_test == 2){ # Cox PH test (Likelihood Ratio)
            
            
            
            p_test = Cox_ph_test(tox_treat_vec, safe_censor_treat, tox_con_vec, safe_censor_con, zero_censor = FALSE)
            
            #browser()
            
            if (p_test < safe_conf_2){ # wth high prob is the treat more toxic
              
              treat_rate = sum(tox_treat_vec)/sum(1 - safe_censor_treat)
              con_rate = sum(tox_con_vec)/sum(1 - safe_censor_con)
              
              if ( treat_rate < con_rate){ # if the treatment has more safety events more frequently
                
  
                
                
                #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
                return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
                #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
                
                #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              }
              
            }
            
             
          }
          else if (safe_test == 3){ # AJE cumulative incidence testing
            
            # Cumualtive incidence test
            
            AJE_censor_treat = numeric(length(tox_treat_vec))
            AJE_censor_con = numeric(length(tox_con_vec))
            
            for (i in 1:length(AJE_censor_treat)){
              AJE_censor_treat[i] = max(safe_censor_treat[i],safe_compete_treat[i] )
            }
            for (i in 1:length(AJE_censor_con)){
              AJE_censor_con[i] = max(safe_censor_con[i],safe_compete_con[i] )
            }
            
            
            
            z_AJE = AJE_test(tox_treat_vec, safe_censor_treat, tox_con_vec, safe_censor_con, zero_censor = FALSE)
            
            #browser()
            
            if (z_AJE > safe_conf_3){ # high negative z test means that the rate of the treatment toxics is faster than the control
              
              
  
              
              
              #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
              #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
              
              #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              
            }
            
          }
          else if (safe_test == 4){ # Cox PH test (Likelihood Ratio)
            
            tau_1 = max(tox_treat_vec)
            tau_2 = max(tox_con_vec)
            
            tau = min(c(tau_1,tau_2))
            
            p_test = RMST_test(tox_treat_vec, safe_censor_treat, tox_con_vec, safe_censor_con, zero_censor = FALSE, tau = tau)
            
            #browser()
            
            if (p_test < -safe_conf_4){ # high negative z test means that the rate of the treatment toxics is faster than the control
              
  
              
              
              #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
              #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
              
              #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              
            }
            
            
          }
          
          
          
        }
      }
      else if (test == 3){
        
        post_prob = bayes_factor_comp(1 + eff_treat_suc, 1+eff_treat_fail, 1 + eff_con_suc, 1 + eff_con_fail, norm_approx = norm_approx) # 1 + because Beta(1,1) prior
        
        #browser()
        
        if (post_prob < conf[i]){ # if we pass the decision boundary
          
  
          
          
          #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
          return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
          #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
          
          #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
          #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
        } 
        else if (fut_stop){
          eff_treat_suc_curtail = eff_treat_suc + round(eff_con*(N - IAn[i])/2) # if it does very well
          eff_con_suc_curtail = eff_con_suc + round(eff_con*(N - IAn[i])/2) # if it does as expected
          
          eff_treat_fail_curtail = eff_treat_fail + round((1 - eff_con)*(N - IAn[i])/2) # if it does very well
          eff_con_fail_curtail = eff_con_fail + round((1 - eff_con)*(N - IAn[i])/2) # if it does as expected
          
          post_prob_curtail = bayes_factor_comp(1 + eff_treat_suc_curtail, 1+eff_treat_fail_curtail, 1 + eff_con_suc_curtail, 1 + eff_con_fail_curtail, norm_approx = norm_approx)
          
          # If the probability of it rejecting the null at the end is unlikely, then will just end the trial early
          if (post_prob_curtail > conf[length(IAn)] + 0.1){
            
  
            
            #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
            #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
            
            
            #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            
          }
        }
        if (safe_stop){
          
          if (safe_test == 1){ # BOP2 esque Bayesian Posterior testing 
            
            
            treat_event_time = sum(tox_treat_vec) # time at risk for the treatment
            treat_num_event = sum(1 - safe_censor_treat) # number of events occurring
            
            con_event_time = sum(tox_con_vec) # time at risk for the treatment
            con_num_event = sum(1 - safe_censor_con) # number of events occurring
            
            safe_prob = Bayes_cts_test(treat_event_time, treat_num_event, con_event_time, con_num_event)
            
            # browser()
            
            if (safe_prob > safe_conf){ # wth high prob is the treat more toxic
              
  
              
              
              #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
              #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
              
              
              #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            }
            
          }
          else if (safe_test == 2){ # Cox PH test (Likelihood Ratio)
            
            
            
            p_test = Cox_ph_test(tox_treat_vec, safe_censor_treat, tox_con_vec, safe_censor_con, zero_censor = FALSE)
            
            if (p_test < safe_conf_2){ # wth high prob is the treat more toxic
              
              treat_rate = sum(tox_treat_vec)/sum(1 - safe_censor_treat)
              con_rate = sum(tox_con_vec)/sum(1 - safe_censor_con)
              
              if ( treat_rate < con_rate){ # if the treatment has more safety events more frequently
                
                
  
                
                
                #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
                return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
                #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
                
                
                #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              }
              
            }
          }
          else if (safe_test == 3){ # AJE cumulative incidence testing
            
            # Cumualtive incidence test
            
            AJE_censor_treat = numeric(length(tox_treat_vec))
            AJE_censor_con = numeric(length(tox_con_vec))
            
            for (i in 1:length(AJE_censor_treat)){
              AJE_censor_treat[i] = max(safe_censor_treat[i],safe_compete_treat[i] )
            }
            for (i in 1:length(AJE_censor_con)){
              AJE_censor_con[i] = max(safe_censor_con[i],safe_compete_con[i] )
            }
            
            z_AJE = AJE_test(tox_treat_vec, safe_censor_treat, tox_con_vec, safe_censor_con, zero_censor = FALSE)
            
            #browser()
            
            if (z_AJE > safe_conf_3){ # high negative z test means that the rate of the treatment toxics is faster than the control
              
              
  
              
              
              #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
              #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
              
              
              
              #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              
            }
            
          }
          else if (safe_test == 4){ # Cox PH test (Likelihood Ratio)
            
            tau_1 = max(tox_treat_vec)
            tau_2 = max(tox_con_vec)
            
            tau = min(c(tau_1,tau_2))
            
            p_test = RMST_test(tox_treat_vec, safe_censor_treat, tox_con_vec, safe_censor_con, zero_censor = FALSE, tau = tau)
            
            #browser()
            
            if (p_test < -safe_conf_4){ # high negative z test means that the rate of the treatment toxics is faster than the control
              
              
  
              
              
              #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
              #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
              
              
              #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              
            }
            
            
          }
        }
      }
      
      else if (test == 4){
        
        # assuming normal then finding 95% approximate CI 
        num_IA = length(IAn)
        Pocock_bounds = numeric(num_IA)
        for (j in 1:num_IA){
          Pocock_bounds[j] = sqrt(num_IA/j)
        }
        
        
        n_treat = eff_treat_suc + eff_treat_fail
        mle_treat = eff_treat_suc/n_treat
        var_treat = (eff_treat_suc*eff_treat_fail) / n_treat^3
        
        n_con = eff_con_suc + eff_con_fail
        mle_con = eff_con_suc/n_con
        var_con = (eff_con_suc*eff_con_fail) / n_con^3
        
        var_diff = var_con + var_treat
        
        # Now, mle_treat - mle_con ~ N(delta, var_diff) under normal assumptions, where delta is the NI margin
        # So, mle_treat - mle_con +- 1.96*sqrt(var_diff) should be a 95% CI
        # If mle_treat - mle_con - 1.96*var_diff > delta, then we have NI 
        
        CI_conf = qnorm(1 - alpha)*Pocock_bounds[i]
        
        if (mle_con - mle_treat + CI_conf*sqrt(var_diff) < NI_margin){ # if we pass the decision boundary # NOTE not currently tuned for multiple IA 
          
  
          
          
          #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
          return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
          #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
          
          
          
          #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
        }
        else if (fut_stop){
          eff_treat_suc_curtail = eff_treat_suc + round(eff_con*(N - IAn[i])/2) # if it does very well
          eff_con_suc_curtail = eff_con_suc + round(eff_con*(N - IAn[i])/2) # if it does as expected
          
          eff_treat_fail_curtail = eff_treat_fail + round((1 - eff_con)*(N - IAn[i])/2) # if it does very well
          eff_con_fail_curtail = eff_con_fail + round((1 - eff_con)*(N - IAn[i])/2) # if it does as expected
          
          
          n_treat_c = eff_treat_suc_curtail + eff_treat_fail_curtail
          mle_treat_c  = eff_treat_suc_curtail/n_treat_c
          var_treat_c = (eff_treat_suc_curtail*eff_treat_fail_curtail) / n_treat_c^3
          
          n_con_c = eff_con_suc_curtail + eff_con_fail_curtail
          mle_con_C = eff_con_suc_curtail/n_con_c
          var_con_c = (eff_con_suc_curtail*eff_con_fail_curtail) / n_con_c^3
          
          var_diff_c = var_con_c + var_treat_c
          
          # Now, mle_treat - mle_con ~ N(delta, var_diff) under normal assumptions, where delta is the NI margin
          # So, mle_treat - mle_con +- 1.96*sqrt(var_diff) should be a 95% CI
          # If mle_treat - mle_con - 1.96*var_diff > delta, then we have NI 
          
          CI_conf = qnorm(1 - alpha)*Pocock_bounds[length(IAn)]
          
          #browser()
          
          # If the probability of it rejecting the null at the end is unlikely, then will just end the trial early
          if (mle_con_C - mle_treat_c + CI_conf*sqrt(var_diff_c) > NI_margin + 0.01){ # if we pass the decision boundary # NOTE not currently tuned for multiple IA 
            
  
            
            
            #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
            #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
            
            
            
            
            #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
          }
        }
        if (safe_stop){
          
          if (safe_test == 1){ # BOP2 esque Bayesian Posterior testing 
            
            
            treat_event_time = sum(tox_treat_vec) # time at risk for the treatment
            treat_num_event = sum(1 - safe_censor_treat) # number of events occurring
            
            con_event_time = sum(tox_con_vec) # time at risk for the treatment
            con_num_event = sum(1 - safe_censor_con) # number of events occurring
            
            safe_prob = Bayes_cts_test(treat_event_time, treat_num_event, con_event_time, con_num_event)
            
            # browser()
            
            if (safe_prob > safe_conf){ # wth high prob is the treat more toxic
              
              
  
              
              #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
              #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
              
              
              
              
              #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
            }
            
          }
          else if (safe_test == 2){ # Cox PH test (Likelihood Ratio)
            
            
            
            p_test = Cox_ph_test(tox_treat_vec, safe_censor_treat, tox_con_vec, safe_censor_con, zero_censor = FALSE)
            
            if (p_test < safe_conf_2){ # wth high prob is the treat more toxic
              
              treat_rate = sum(tox_treat_vec)/sum(1 - safe_censor_treat)
              con_rate = sum(tox_con_vec)/sum(1 - safe_censor_con)
              
              if ( treat_rate < con_rate){ # if the treatment has more safety events more frequently
                
  
                
                
                #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
                return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
                #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
                
                
                
                
                #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
              }
              
            }
          }
          else if (safe_test == 3){ # AJE cumulative incidence testing
            
            # Cumualtive incidence test
            
            AJE_censor_treat = numeric(length(tox_treat_vec))
            AJE_censor_con = numeric(length(tox_con_vec))
            
            for (j in 1:length(AJE_censor_treat)){
              AJE_censor_treat[j] = max(safe_censor_treat[j],safe_compete_treat[j] )
            }
            for (j in 1:length(AJE_censor_con)){
              AJE_censor_con[j] = max(safe_censor_con[j],safe_compete_con[j] )
            }
            
            #print(tox_treat_vec)
            #print(safe_censor_treat)
            #print(tox_con_vec)
            #print(safe_censor_con)
            
            # Gives an error when both are fully
            
            if (identical(safe_censor_treat, safe_censor_con) == FALSE){
              # We only need to consider safety stopping if the data isn't identical
            
              z_AJE = AJE_test(tox_treat_vec, safe_censor_treat, tox_con_vec, safe_censor_con, zero_censor = FALSE)
              
              #browser()
              
              
              if (z_AJE > safe_conf_3){ # high negative z test means that the rate of the treatment toxics is faster than the control
                
    
                #browser()
                
                #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
                return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
                #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
                
                
                
                #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
                
              }
            }
            
          }
          else if (safe_test == 4){ # Cox PH test (Likelihood Ratio)
            
            tau_1 = max(tox_treat_vec)
            tau_2 = max(tox_con_vec)
            
            tau = min(c(tau_1,tau_2))
            
            p_test = RMST_test(tox_treat_vec, safe_censor_treat, tox_con_vec, safe_censor_con, zero_censor = FALSE, tau = tau)
            
            #browser()
            
            if(is.na(p_test)){
              # Happens when the sets are the same, if they're not we need to investigate , otherwise we pass the test
              if (identical(safe_censor_treat, safe_censor_con) == FALSE){
                browser()
              }
              
            }
            else{
              
            
              if (p_test < -safe_conf_4){ # high negative z test means that the rate of the treatment toxics is faster than the control
                
    
                
                
                #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
                return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
                #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
                
                
                
                
                #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
                
              }
            }
            
            
          }
        }
      }
    }
    
    # Number of toxic events in general 
    num_toxics = sum(1 - safe_censor_con) + sum(1 - safe_censor_treat)
    num_toxics_con = sum(1 - safe_censor_con) 
    num_toxics_treat = sum(1 - safe_censor_treat)
    
    # Average time until a toxic event in the control if it occurs
    mean_toxic_time_con = mean(tox_con_vec[which(safe_censor_con == 0)])
    
    # Average time until a toxic event in the treatment if it occurs
    mean_toxic_time_treat = mean(tox_treat_vec[which(safe_censor_treat == 0)])
    
    
    #return(c("Null Reject" = 1,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
    return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = num_toxics, "Average time until treatment toxic" = mean_toxic_time_treat, "Average time until control toxic" = mean_toxic_time_con, "Treatment toxics" = num_toxics_treat, "Control toxics" = num_toxics_con))
    #return(c(1,eff_treat_suc, eff_treat_fail, eff_con_suc, eff_con_fail,IAn[i]))
    
    
    
    
    #return(c("Null Reject" = 0,"Treatment success" = eff_treat_suc,"Treatment failures" = eff_treat_fail,"Control success" = eff_con_suc,"Control failures" = eff_con_fail,"Sample size" = IAn[i], "Total number of  successes" = eff_treat_suc + eff_con_suc, "Total number of toxics" = sum(tox_con_vec) + sum(tox_treat_vec)))
    
  }
  
  Bayes_cts_test = function(event_time_treat, num_events_treat, event_time_con, num_events_con, prior_a = 0.001, prior_b = 0.001){
    ## Does the BOP2 model of IG(0.001,0.001) prior, Exponential likelhiood to 
      # get an IG posterior, which we calculate the probability of being better of here
    
    # We have two variables, treat: IG(num_events_treat + prior_a, event_time_treat + prior_b)
      # control: IG(num_events_con + prior_a, event_time_con + prior_b)
    
    a_t = num_events_treat + prior_a
    b_t = event_time_treat + prior_b
    
    a_c = num_events_con + prior_a
    b_c = event_time_con + prior_b
    
    # Then, we want to utilise the fact that P(x > Y ) = P( 1/X < 1 / Y), where X,Y are gamma
      # 1/X,1/Y are IG. Then, the probability of a Gamma(a_1,b_1) > Gamma(a_2,b_2) = P(Beta(a_2,a_1) < b_2/(b_1 + b_2))
    
    post_prob = pbeta(b_c/(b_c + b_t), a_c, a_t)
    
    return(post_prob)
    
    
    
  }
  
  Cox_ph_test = function(event_time_treat, censor_treat, event_time_con, censor_con, zero_censor = FALSE){
    # Zero censor is a variable which determines whether 0 is used to mean an event or censored
    # Surv() takes 0 to be censored, so zero_censor == false means we have to flip
    
    
    
    
    time_vec_df = c(event_time_treat,event_time_con)
    
    censor_vec_df = c(censor_treat,censor_con)
    
    if (zero_censor == FALSE){
      # Surv() needs zeros to denote censoring, so this flips it 
      
      censor_vec_df = 1 - censor_vec_df
    }
    
    treat_vec_df = c(rep(1, length(censor_treat)),rep(0, length(censor_con)) )
    
    test_df = data.frame("Time" = time_vec_df, "Censor" = censor_vec_df, "Treatment" = factor(treat_vec_df))
    
    surv_model = coxph(Surv(Time, Censor) ~ Treatment, data = test_df) 
    
    sum_surv_model = summary(surv_model)
    
    lrt_test = sum_surv_model$logtest[3] # likelihood ratio test
    
    return(lrt_test)
    
  }
  
  AJE_test = function(event_time_treat, censor_treat, event_time_con, censor_con, zero_censor = FALSE){
    # Does the cumulative incidence test mentioned in Fine and Gray 1999
    # Robust to competing risks
    
    # censor here should note 0,1,2, the conditions of event
    
    time_vec_df = c(event_time_treat,event_time_con)
    
    censor_vec_df = c(censor_treat,censor_con)
    
    treat_vec_df = c(rep(1, length(censor_treat)),rep(0, length(censor_con)) )
    
    test_df = data.frame("Time" = time_vec_df, "Censor" = censor_vec_df, "Treatment" = factor(treat_vec_df))
    
    censor_code = 1
    fail_code = 0
    
    if (zero_censor){
      censor_code = 0
      fail_code = 1
    }
    
    crr_results = crr(ftime = test_df$Time, fstatus = test_df$Censor, cov1 = test_df$Treatment, cencode =  censor_code, failcode = fail_code)
    
    z_value = summary(crr_results)$coef[4]
    
    #browser()
    
    return(z_value)
  }
  
  RMST_test = function(event_time_treat, censor_treat, event_time_con, censor_con, tau = 2.5, zero_censor = TRUE){
    # Does Restricted Mean Survival Time testing
    # Robust to proportional hazards
    
    
    
    time_vec_df = c(event_time_treat,event_time_con)
    
    censor_vec_df = c(censor_treat,censor_con)
    
    if (zero_censor == FALSE){
      # Surv() needs zeros to denote censoring, so this flips it 
      
      censor_vec_df = 1 - censor_vec_df
    }
    
    treat_vec_df = c(rep(1, length(censor_treat)),rep(0, length(censor_con)) )
    
    test_df = data.frame("Time" = time_vec_df, "Censor" = censor_vec_df, "Treatment" = factor(treat_vec_df))
    
    
    rmst_result = test_df |>
      rmst(
        var_label_tte = "Time",
        var_label_event = "Censor",
        var_label_group = "Treatment",
        tau = tau,
        reference = "0"
      )
    
    return(rmst_result$z)
  }
  
  ## Parameter record
  # alpha = 0.045 gives 0.07 type one error without any safety stopping
  # safe_conf = 0.75 gives a type one of 0.0497 (Safe_test = 1)
  # safe_conf_2 = 0.45 gives a type one error of 0.0498
  # safe_conf_3 = 0.7 gives a type one of 0.0501 (safe_test =3)
  # safe_conf_4 = 0.6, gives a type one of 0.0508 (safe_test = 4) 
  
  OC_gen_cts = function(eff_treat = 0.96, tox_treat = 0.05, eff_con = 0.96, tox_con = 0.15, IAn = c(80,160,240,320), N = 320, n_sim = 10000, paired_seed = FALSE, set_seed = 2026, NI_margin = 0.06, conf = c(0.99,0.98,0.97,0.96), test = 4, alpha = 0.045, RAR = 2, w = 0.5, sigma = 0.1, fut_stop = TRUE, safe_stop = TRUE, safe_conf = 0.75, safe_conf_2 =0.45, safe_conf_3 = 0.7, safe_conf_4 = 0.6, safe_test = 1, norm_approx = TRUE, PH_true = TRUE){
    
    temp_vec = NI_trial_RAR_cts(eff_treat = eff_treat, tox_treat = tox_treat, eff_con = eff_con, tox_con = tox_con, N = N, IAn = IAn)
    result_vec = numeric(length(temp_vec))  # so that this still works if we change the length of the output
    
    for (i in 1:n_sim){
      
      if (paired_seed){
        set.seed(set_seed + i) # makes it so that the inter comparison variance is reduced, every set receives the same set of patients
      }
      
      result_vec = result_vec + NI_trial_RAR_cts(eff_treat = eff_treat, tox_treat = tox_treat, eff_con = eff_con, tox_con = tox_con, N = N, IAn = IAn, NI_margin = NI_margin, conf = conf, test = test, alpha = alpha, RAR = RAR, w = w, sigma = sigma, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, safe_conf_2 =safe_conf_2, safe_conf_3 = safe_conf_3, safe_conf_4 = safe_conf_4, safe_test = safe_test, norm_approx = norm_approx, PH_true = PH_true)
  
      
      
      
    }
    
    prop = (result_vec[2] + result_vec[3])/result_vec[6]
    return(c(result_vec/n_sim,"Proportion" = prop))
    
  }
  
  safety_opt = function(iter = 10, target_alpha = 0.05, safe_test = 1, eff_treat = 0.96, tox_treat = 0.15, eff_con = 0.96, tox_con = 0.15, IAn = c(80,160,240,320), N = 320, NI_margin = 0.06, conf = c(0.99,0.98,0.97,0.96), test = 4, alpha = 0.045, RAR = 2, w = 0.5, sigma = 0.1, fut_stop = TRUE, safe_stop = TRUE, norm_approx = TRUE, PH_true = TRUE, binary = FALSE){
    
    # defining them all here so we can just copy and paste the same function without errors
    safe_conf = 1
    safe_conf_2 = 1
    safe_conf_3 = 1
    safe_conf_4 = 1
    
    safe_conf_up = 1
    safe_conf_low = 0
    
    if (binary){
      
      for (i in 1:iter){
        
        safe_conf = (safe_conf_up + safe_conf_low) / 2
        
        temp_alpha = OC_gen(eff_treat = eff_con - NI_margin, tox_treat = tox_treat, eff_con = eff_con, tox_con = tox_con, N = N, IAn = IAn, NI_margin = NI_margin, conf = conf, test = test, alpha = alpha, RAR = RAR, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = rep(safe_conf, length(IAn)), norm_approx = norm_approx)[1]
        
        if (temp_alpha < 0.975*target_alpha){ # alpha too low, test needs to be more lenient
          # so the confidence goes towards 1
          
          
          safe_conf_up = safe_conf
          
        } else if (temp_alpha > 1.025*target_alpha){ # more strict as we don't want to be above the type one error
          
          safe_conf_low = safe_conf
          
          
        } # if neither condition is met, we do not change anything, as we are close to 0.05
        
      }
      
      return(safe_conf)
      
      
    }
    
    if (safe_test == 1){
      
      
      for (i in 1:iter){
        
        safe_conf = (safe_conf_up + safe_conf_low) / 2
        
        temp_alpha = OC_gen_cts(eff_treat = eff_treat - NI_margin, tox_treat = tox_treat, eff_con = eff_con, tox_con = tox_con, N = N, IAn = IAn, NI_margin = NI_margin, conf = conf, test = test, alpha = alpha, RAR = RAR, w = w, sigma = sigma, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, safe_conf_2 =safe_conf_2, safe_conf_3 = safe_conf_3, safe_conf_4 = safe_conf_4, safe_test = safe_test, norm_approx = norm_approx, PH_true = PH_true)[1]
        
        if (temp_alpha < 0.975*target_alpha){ # alpha too low, test needs to be more lenient
          # so the confidence goes towards 1
          
          safe_conf_low = safe_conf
          
        } else if (temp_alpha > 1.025*target_alpha){ # more strict as we don't want to be above the type one error
          
  
          safe_conf_up = safe_conf
          
        } # if neither condition is met, we do not change anything, as we are close to 0.05
        
      }
      
      return(safe_conf)
      
    }
    if (safe_test == 2){
      
      
      for (i in 1:iter){
        
        safe_conf_2 = (safe_conf_up + safe_conf_low) / 2
        
        temp_alpha = OC_gen_cts(eff_treat = eff_treat - NI_margin, tox_treat = tox_treat, eff_con = eff_con, tox_con = tox_con, N = N, IAn = IAn, NI_margin = NI_margin, conf = conf, test = test, alpha = alpha, RAR = RAR, w = w, sigma = sigma, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, safe_conf_2 =safe_conf_2, safe_conf_3 = safe_conf_3, safe_conf_4 = safe_conf_4, safe_test = safe_test, norm_approx = norm_approx, PH_true = PH_true)[1]
        
        if (temp_alpha < 0.975*target_alpha){ # alpha too low, test needs to be more lenient
          # so the confidence goes towards 1
          
          safe_conf_up = safe_conf_2
          
        } else if (temp_alpha > 1.025*target_alpha){ # more strict as we don't want to be above the type one error
          
          safe_conf_low = safe_conf_2
          
        } # if neither condition is met, we do not change anything, as we are close to 0.05
        
      }
      
      return(safe_conf_2)
      
    }
    if (safe_test == 3){
      
      
      for (i in 1:iter){
        
        safe_conf_3 = (safe_conf_up + safe_conf_low) / 2
        
        temp_alpha = OC_gen_cts(eff_treat = eff_treat - NI_margin, tox_treat = tox_treat, eff_con = eff_con, tox_con = tox_con, N = N, IAn = IAn, NI_margin = NI_margin, conf = conf, test = test, alpha = alpha, RAR = RAR, w = w, sigma = sigma, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, safe_conf_2 =safe_conf_2, safe_conf_3 = safe_conf_3, safe_conf_4 = safe_conf_4, safe_test = safe_test, norm_approx = norm_approx, PH_true = PH_true)[1]
        
        if (temp_alpha < 0.975*target_alpha){ # alpha too low, test needs to be more lenient
          # so the confidence goes towards 1
          safe_conf_low = safe_conf_3
          
        } else if (temp_alpha > 1.025*target_alpha){ # more strict as we don't want to be above the type one error
          
          safe_conf_up = safe_conf_3
          
        } # if neither condition is met, we do not change anything, as we are close to 0.05
        
      }
      
      return(safe_conf_3)
      
    }
    if (safe_test == 4){
      
      
      for (i in 1:iter){
        
        safe_conf_4 = (safe_conf_up + safe_conf_low) / 2
        
        temp_alpha = OC_gen_cts(eff_treat = eff_treat - NI_margin, tox_treat = tox_treat, eff_con = eff_con, tox_con = tox_con, N = N, IAn = IAn, NI_margin = NI_margin, conf = conf, test = test, alpha = alpha, RAR = RAR, w = w, sigma = sigma, fut_stop = fut_stop, safe_stop = safe_stop, safe_conf = safe_conf, safe_conf_2 =safe_conf_2, safe_conf_3 = safe_conf_3, safe_conf_4 = safe_conf_4, safe_test = safe_test, norm_approx = norm_approx, PH_true = PH_true)[1]
        
        if (temp_alpha < 0.975*target_alpha){ # alpha too low, test needs to be more lenient
          # so the confidence goes towards 1
          safe_conf_low = safe_conf_4
          
        } else if (temp_alpha > 1.025*target_alpha){ # more strict as we don't want to be above the type one error
          
          safe_conf_up = safe_conf_4
          
        } # if neither condition is met, we do not change anything, as we are close to 0.05
        
      }
      
      return(safe_conf_4)
      
    }
    
  }
  
  heatmap_plot = function(w = 0.5, sigma = 0.1){
    ## Makes a heatmap plot of 4 different RAR and puts them on the same plot
    
    n = 160 # total per arm at end of trial without RAR
    
    base_success = round(n*0.9) # 0.9 is the experimental success rate
    base_toxic = round(n*0.5) # 0.5 is the control toxic rate
    
    suc_vec = seq(base_success - 10, base_success + 10, 1)
    tox_vec = rep(seq(base_toxic - 10, base_toxic + 10, 1), length(suc_vec))
    suc_vec = rep(suc_vec, each = length(suc_vec))
    
    n_sim = length(suc_vec)
    
    RITS_df = data.frame(suc_vec,tox_vec,allo_prob = n_sim)
    G_RITS_df = data.frame(suc_vec,tox_vec,allo_prob = n_sim)
    SAFER_df = data.frame(suc_vec,tox_vec,allo_prob = n_sim)
    BRAR_df = data.frame(suc_vec,tox_vec,allo_prob = n_sim)
    
    for (i in 1:n_sim){
      # RITS
      
      temp_suc = RITS_df$suc_vec[i]
      temp_tox = RITS_df$tox_vec[i]
      
      allo_prob_eff = 1 - rcpp_exact_beta(1 + temp_suc,1 + n - temp_suc,1 + base_success,1 + n - base_success)
      allo_prob_tox = rcpp_exact_beta(1 + temp_tox, 1 + n -temp_tox, 1 + base_toxic, 1 + n - base_toxic)
      
      allo_prob = w*allo_prob_eff + (1 - w)*allo_prob_tox
      
      allo_prob = allo_prob^(1/(2))/(allo_prob^(1/(2)) + (1-allo_prob)^(1/(2))) # tuning with n/2N
      
      RITS_df$allo_prob[i] = allo_prob
      
      
    }
    for (i in 1:n_sim){
      # Gaussian RITS
      
      temp_suc = G_RITS_df$suc_vec[i]
      temp_tox = G_RITS_df$tox_vec[i]
      
      allo_prob_eff = 1 - rcpp_exact_beta(1 + temp_suc,1 + n - temp_suc,1 + base_success,1 + n - base_success)
      allo_prob_tox = rcpp_exact_beta(1 + temp_tox, 1 + n -temp_tox, 1 + base_toxic, 1 + n - base_toxic)
      
      
      
      g_x = exp( - ((allo_prob_eff - 0.5)^2) / sigma )
      
      allo_prob = (1 - g_x)*allo_prob_eff + g_x*allo_prob_tox
      
      allo_prob = allo_prob^(1/(2))/(allo_prob^(1/(2)) + (1-allo_prob)^(1/(2))) # tuning with n/2N
      
      G_RITS_df$allo_prob[i] = allo_prob
      
    }
    for (i in 1:n_sim){
      # SAFER
  
      temp_suc = SAFER_df$suc_vec[i]
      temp_tox = SAFER_df$tox_vec[i]
      
      allo_prob_eff = 1 - rcpp_exact_beta(1 + temp_suc,1 + n - temp_suc,1 + base_success,1 + n - base_success)
      allo_prob_tox = rcpp_exact_beta(1 + temp_tox, 1 + n -temp_tox, 1 + base_toxic, 1 + n - base_toxic)
      
      eff_treat_vec = sample(c(rep(1,temp_suc),rep(0, n - temp_suc)))
      eff_con_vec = sample(c(rep(1,base_success),rep(0, n - base_success)))
      
      allo_prob = safer_allo(eff_treat_vec,eff_con_vec, allo_prob_tox) 
      # will try without tuning, not expecting SAFER to be hugely relevant for now
      
      SAFER_df$allo_prob[i] = allo_prob
      
    }
    for (i in 1:n_sim){
      # BRAR just based on efficacy
      temp_suc = BRAR_df$suc_vec[i]
      temp_tox = BRAR_df$tox_vec[i]
      
      allo_prob = 1 - rcpp_exact_beta(1 + temp_suc,1 + n - temp_suc,1 + base_success,1 + n - base_success)
  
      
      allo_prob = allo_prob^(1/(2))/(allo_prob^(1/(2)) + (1-allo_prob)^(1/(2))) # tuning with n/2N
      
      BRAR_df$allo_prob[i] = allo_prob
    }
    
    #browser()
    par(mfrow = c(2, 2))
    
    RITS_plot = ggplot(RITS_df, aes(suc_vec, tox_vec, fill= allo_prob)) + 
      geom_tile() +
      scale_fill_gradientn(colours=c("darkred", "red", "pink", 
                                     "white",
                                     "lightblue", "blue", "darkblue"),
                           values=(c(0, 0.2, 0.4,
                                            0.5,
                                            0.6, 0.8, 1)),
                           guide="colorbar") +
      labs(
        x = "Number of successes in the experimental arm",
        y  = "Number of toxics in the experimental arm",
        title = "Allocation to the experimental- RITS",
        subtitle = "Control: 144 successes, 80 toxics",
        fill = "Allocation prob") 
      
    
    print(RITS_plot)
    
    G_RITS_plot = ggplot(G_RITS_df, aes(suc_vec, tox_vec, fill= allo_prob)) + 
      geom_tile() +
      scale_fill_gradientn(colours=c("darkred", "red", "pink", 
                                     "white",
                                     "lightblue", "blue", "darkblue"),
                           values=(c(0, 0.2, 0.4,
                                     0.5,
                                     0.6, 0.8, 1)),
                           guide="colorbar") +
      labs(
        x = "Number of successes in the experimental arm",
        y  = "Number of toxics in the experimental arm",
        title = "Allocation to the experimental- Gaussian RITS",
        subtitle = "Control: 144 successes, 80 toxics",
        fill = "Allocation prob") 
    
    
    print(G_RITS_plot)
    
    SAFER_plot = ggplot(SAFER_df, aes(suc_vec, tox_vec, fill= allo_prob)) + 
      geom_tile() +
      scale_fill_gradientn(colours=c("darkred", "red", "pink", 
                                     "white",
                                     "lightblue", "blue", "darkblue"),
                           values=(c(0, 0.2, 0.4,
                                     0.5,
                                     0.6, 0.8, 1)),
                           guide="colorbar") +
      labs(
        x = "Number of successes in the experimental arm",
        y  = "Number of toxics in the experimental arm",
        title = "Allocation to the experimental- SAFER",
        subtitle = "Control: 144 successes, 80 toxics",
        fill = "Allocation prob") 
    
    
    print(SAFER_plot)
    
    BRAR_plot = ggplot(BRAR_df, aes(suc_vec, tox_vec, fill= allo_prob)) + 
      geom_tile() +
      scale_fill_gradientn(colours=c("darkred", "red", "pink", 
                                     "white",
                                     "lightblue", "blue", "darkblue"),
                           values=(c(0, 0.2, 0.4,
                                     0.5,
                                     0.6, 0.8, 1)),
                           guide="colorbar") +
      labs(
        x = "Number of successes in the experimental arm",
        y  = "Number of toxics in the experimental arm",
        title = "Allocation to the experimental- BRAR",
        subtitle = "Control: 144 successes, 80 toxics",
        fill = "Allocation prob") 
    
    
    print(BRAR_plot)
    
    par(mfrow = c(1, 1))
    
    grid.arrange(RITS_plot, G_RITS_plot, SAFER_plot, BRAR_plot, nrow = 2, ncol = 2)
    
  }
  
  if (FALSE){
    
    safe_df_H0 = data.frame("Alpha" = numeric(4),"ESS" = numeric(4) , "Proportion" = numeric(4), "Total_Toxics" = numeric(4), "Safe_test" = c(1,2,3,4) )
    safe_df_H1 = data.frame("Power" = numeric(4),"ESS" = numeric(4) , "Proportion" = numeric(4), "Total_Toxics" = numeric(4), "Safe_test" = c(1,2,3,4) )
    
    for (i in 1:4){
      
      print(i)
      
      temp_H0 = OC_gen_cts(eff_con = 0.96, tox_con = 0.15, eff_treat = 0.9, tox_treat = 0.15, safe_test = i, paired_seed = TRUE, n_sim = 10000, NI_margin = 0.06, fut_stop = TRUE)
      temp_H1 = OC_gen_cts(eff_con = 0.96, tox_con = 0.15, eff_treat = 0.96, tox_treat = 0.05, safe_test = i, paired_seed = TRUE, n_sim = 10000, NI_margin = 0.06, fut_stop = TRUE)
      
      
      safe_df_H0$Alpha[i] = temp_H0[1]
      safe_df_H0$ESS[i] = temp_H0[6]
      safe_df_H0$Proportion[i] = temp_H0[13]
      safe_df_H0$Total_Toxics[i] = temp_H0[8]
      
      safe_df_H1$Power[i] = temp_H1[1]
      safe_df_H1$ESS[i] = temp_H1[6]
      safe_df_H1$Proportion[i] = temp_H1[13]
      safe_df_H1$Total_Toxics[i] = temp_H1[8]
      
      print("Normal")
      print(safe_df_H0)
      print(safe_df_H1)
      
      
    }
  }
  
  if (FALSE){
    
    safe_df_H0 = data.frame("Alpha" = numeric(4),"ESS" = numeric(4) , "Proportion" = numeric(4), "Total_Toxics" = numeric(4), "Safe_test" = c(1,2,3,4) )
    safe_df_H1 = data.frame("Power" = numeric(4),"ESS" = numeric(4) , "Proportion" = numeric(4), "Total_Toxics" = numeric(4), "Safe_test" = c(1,2,3,4) )
    
    for (i in 1:4){
      
      print(i)
      
      temp_H0 = OC_gen_cts(eff_con = 0.96, tox_con = 0.15, eff_treat = 0.9, tox_treat = 0.05, safe_test = i, paired_seed = TRUE, n_sim = 10000, NI_margin = 0.06, fut_stop = TRUE)
      temp_H1 = OC_gen_cts(eff_con = 0.96, tox_con = 0.15, eff_treat = 0.96, tox_treat = 0.15, safe_test = i, paired_seed = TRUE, n_sim = 10000, NI_margin = 0.06, fut_stop = TRUE)
      
      safe_df_H0$Alpha[i] = temp_H0[1]
      safe_df_H0$ESS[i] = temp_H0[6]
      safe_df_H0$Proportion[i] = temp_H0[13]
      safe_df_H0$Total_Toxics[i] = temp_H0[8]
      
      safe_df_H1$Power[i] = temp_H1[1]
      safe_df_H1$ESS[i] = temp_H1[6]
      safe_df_H1$Proportion[i] = temp_H1[13]
      safe_df_H1$Total_Toxics[i] = temp_H1[8]
      
      print("Worst case safety stopping")
      print(safe_df_H0)
      print(safe_df_H1)
      
      
    }
  }
  
  if (FALSE){
    
    safe_df_H0 = data.frame("Alpha" = numeric(4),"ESS" = numeric(4) , "Proportion" = numeric(4), "Total_Toxics" = numeric(4), "Safe_test" = c(1,2,3,4) )
    safe_df_H1 = data.frame("Power" = numeric(4),"ESS" = numeric(4) , "Proportion" = numeric(4), "Total_Toxics" = numeric(4), "Safe_test" = c(1,2,3,4) )
    
    for (i in 1:1){
      
      print(i)
      
      temp_H0 = OC_gen_cts(eff_con = 0.96, tox_con = 0.5, eff_treat = 0.9, tox_treat = 0.5, safe_test = i, paired_seed = TRUE, n_sim = 100000, NI_margin = 0.06, fut_stop = TRUE)
      temp_H1 = OC_gen_cts(eff_con = 0.96, tox_con = 0.5, eff_treat = 0.96, tox_treat = 0.4, safe_test = i, paired_seed = TRUE, n_sim = 100000, NI_margin = 0.06, fut_stop = TRUE)
      
      safe_df_H0$Alpha[i] = temp_H0[1]
      safe_df_H0$ESS[i] = temp_H0[6]
      safe_df_H0$Proportion[i] = temp_H0[13]
      safe_df_H0$Total_Toxics[i] = temp_H0[8]
      
      safe_df_H1$Power[i] = temp_H1[1]
      safe_df_H1$ESS[i] = temp_H1[6]
      safe_df_H1$Proportion[i] = temp_H1[13]
      safe_df_H1$Total_Toxics[i] = temp_H1[8]
      
      print("Medium toxic rate")
      print(safe_df_H0)
      print(safe_df_H1)
      
      
    }
  }
  
  if (FALSE){
    
    safe_df_H0 = data.frame("Alpha" = numeric(4),"ESS" = numeric(4) , "Proportion" = numeric(4), "Total_Toxics" = numeric(4), "Safe_test" = c(1,2,3,4) )
    safe_df_H1 = data.frame("Power" = numeric(4),"ESS" = numeric(4) , "Proportion" = numeric(4), "Total_Toxics" = numeric(4), "Safe_test" = c(1,2,3,4) )
    
    for (i in 1:4){
      
      print(i)
      
      temp_H0 = OC_gen_cts(eff_con = 0.96, tox_con = 0.95, eff_treat = 0.9, tox_treat = 0.95, safe_test = i, paired_seed = TRUE, n_sim = 10000, NI_margin = 0.06, fut_stop = TRUE)
      temp_H1 = OC_gen_cts(eff_con = 0.96, tox_con = 0.95, eff_treat = 0.96, tox_treat = 0.85, safe_test = i, paired_seed = TRUE, n_sim = 10000, NI_margin = 0.06, fut_stop = TRUE)
      
      safe_df_H0$Alpha[i] = temp_H0[1]
      safe_df_H0$ESS[i] = temp_H0[6]
      safe_df_H0$Proportion[i] = temp_H0[13]
      safe_df_H0$Total_Toxics[i] = temp_H0[8]
      
      safe_df_H1$Power[i] = temp_H1[1]
      safe_df_H1$ESS[i] = temp_H1[6]
      safe_df_H1$Proportion[i] = temp_H1[13]
      safe_df_H1$Total_Toxics[i] = temp_H1[8]
      
      print("High toxic rate")
      print(safe_df_H0)
      print(safe_df_H1)
      
      
    }
  }
  
  if (FALSE){
    
    safe_df_H0_NPH = data.frame("Alpha" = numeric(4),"ESS" = numeric(4) , "Proportion" = numeric(4), "Total_Toxics" = numeric(4), "Safe_test" = c(1,2,3,4) )
    safe_df_H1_NPH = data.frame("Power" = numeric(4),"ESS" = numeric(4) , "Proportion" = numeric(4), "Total_Toxics" = numeric(4), "Safe_test" = c(1,2,3,4) )
    
    for (i in 1:4){
      
      print(i)
      
      temp_H0 = OC_gen_cts(eff_treat = 0.9, tox_treat = 0.15, safe_test = i, paired_seed = TRUE, n_sim = 100000, PH_true = FALSE)
      temp_H1 = OC_gen_cts(safe_test = i, paired_seed = TRUE, n_sim = 100000, PH_true = FALSE)
      
      safe_df_H0_NPH$Alpha[i] = temp_H0[1]
      safe_df_H0_NPH$ESS[i] = temp_H0[6]
      safe_df_H0_NPH$Proportion[i] = temp_H0[13]
      safe_df_H0_NPH$Total_Toxics[i] = temp_H0[8]
      
      safe_df_H1_NPH$Power[i] = temp_H1[1]
      safe_df_H1_NPH$ESS[i] = temp_H1[6]
      safe_df_H1_NPH$Proportion[i] = temp_H1[13]
      safe_df_H1_NPH$Total_Toxics[i] = temp_H1[8]
      
      
    }
  }

}

#### Paper figures ####

## Figure 1
fig_1 = FALSE
if (fig_1){
  set.seed = 2026
  
  plot_ss(0.9, n_sim = 10000)
}

## Table 1
tab_1_opt = FALSE
tab_1 = FALSE
if (tab_1_opt){
  set.seed(2026)
  para_bop2 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1) # BOP2
  para_bf = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3) # Bayes Factor
  para_na = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4) # Normal Approx
}else{
  para_bop2 = c(0.9937744, 0.9823486, 0.9818115, 0.9702881)
  para_bf = c(0.004935303, 0.009997803, 0.009997803, 0.114282227)
  para_na = 0.0308
  
}
if (tab_1){
  set.seed(2026)
  row_H0_bop2 = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 1, conf = para_bop2, paired_seed = TRUE)
  row_H1_bop2 = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 1, conf = para_bop2, paired_seed = TRUE)
  
  row_H0_bf = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 3, conf = para_bf, paired_seed = TRUE)
  row_H1_bf = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 3, conf = para_bf, paired_seed = TRUE)
  
  row_H0_na = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 4, alpha = para_na, paired_seed = TRUE)
  row_H1_na = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 4, alpha = para_na, paired_seed = TRUE)
}

## tab+fig 2 

tab_2 = FALSE
fig_2 = FALSE
tab_fig_2_opt = FALSE
if (tab_fig_2_opt){
  set.seed(2026)
  para_bop2 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1) # BOP2
  para_bf = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3) # Bayes Factor
  para_na = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4) # Normal Approx
  
  print(1)
  
  para_bop2_safe = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, safe_stop = TRUE) # BOP2
  para_bf_safe = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, safe_stop = TRUE) # Bayes Factor
  para_na_safe = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, safe_stop = TRUE) # Normal Approx
  
  print(2)
  
  para_bop2_fut = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, fut_stop = TRUE) # BOP2
  para_bf_fut = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, fut_stop = TRUE) # Bayes Factor
  para_na_fut = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, fut_stop = TRUE) # Normal Approx
  
  print(3)
  
  para_bop2_safe_fut = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, safe_stop = TRUE, fut_stop = TRUE) # BOP2
  para_bf_safe_fut = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, safe_stop = TRUE, fut_stop = TRUE) # Bayes Factor
  para_na_safe_fut = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, safe_stop = TRUE, fut_stop = TRUE) # Normal Approx
  
  print(4)
  
}else{
  para_bop2 = c(0.9937744, 0.9823486, 0.9818115, 0.9702881)
  para_bf = c(0.004935303, 0.009997803, 0.009997803, 0.114282227)
  para_na = 0.0308
  
  para_bop2_safe = c(0.9921631, 0.9802490, 0.9700928, 0.9558838)
  para_bf_safe = c(0.008116943, 0.009997803, 0.100219727, 0.101098633)
  para_na_safe = 0.0425293
  
  para_bop2_fut = c(0.9928467, 0.9823486, 0.9790283, 0.9705811)
  para_bf_fut = c(0.007176514, 0.009997803, 0.100219727, 0.032478027)
  para_na_fut = 0.03208008
  
  para_bop2_safe_fut = c(0.9922119, 0.9799561, 0.9681885, 0.9558350)
  para_bf_safe_fut = c(0.005493408, 0.009997803, 0.109887695, 0.137573242)
  para_na_safe_fut = 0.0425293
  
  
  
}
if (fig_2){
  
  safe_stop_vec = c(0,1,0,1)
  fut_stop_vec = c(0,0,1,1)
  
  stopping_vec_df = rep(c("No futility or safety","Only safety","Only futility","Both futility and safety"),3)
  test_vec_df = rep(c("Post Prob","Bayes Factor","CI testing"), each = 4)
  
  plot_df_h0 = data.frame("Stopping" = factor(stopping_vec_df), "Test" = factor(test_vec_df), ESS = numeric(12))
  plot_df_h1 = data.frame("Stopping" = factor(stopping_vec_df), "Test" = factor(test_vec_df), ESS = numeric(12))
  
  bop2_ess_vec_h0 = numeric(4)
  bop2_ess_vec_h1 = numeric(4)
  
  bf_ess_vec_h0 = numeric(4)
  bf_ess_vec_h1 = numeric(4)
  
  na_ess_vec_h0 = numeric(4)
  na_ess_vec_h1 = numeric(4)
  
  bop2_conf_list = list(para_bop2,para_bop2_safe,para_bop2_fut,para_bop2_safe_fut)
  bf_conf_list = list(para_bf,para_bf_safe,para_bf_fut,para_bf_safe_fut)
  na_conf_list = list(para_na,para_na_safe,para_na_fut,para_na_safe_fut)
  
  for (i in 1:4){
    # paired seed sets the seed default in each one to 2026 
    
    plot_df_h0$ESS[i] = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 1, conf = bop2_conf_list[[i]], paired_seed = TRUE, safe_stop = safe_stop_vec[i], fut_stop = fut_stop_vec[i])[6]
    plot_df_h1$ESS[i] = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 1, conf = bop2_conf_list[[i]], paired_seed = TRUE, safe_stop = safe_stop_vec[i], fut_stop = fut_stop_vec[i])[6]
    
    plot_df_h0$ESS[i+4] = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 3, conf = bf_conf_list[[i]], paired_seed = TRUE, safe_stop = safe_stop_vec[i], fut_stop = fut_stop_vec[i])[6]
    plot_df_h1$ESS[i+4] = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 3, conf = bf_conf_list[[i]], paired_seed = TRUE, safe_stop = safe_stop_vec[i], fut_stop = fut_stop_vec[i])[6]
    
    plot_df_h0$ESS[i+8] = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 4, alpha = na_conf_list[[i]], paired_seed = TRUE, safe_stop = safe_stop_vec[i], fut_stop = fut_stop_vec[i])[6]
    plot_df_h1$ESS[i+8] = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 4, alpha = na_conf_list[[i]], paired_seed = TRUE, safe_stop = safe_stop_vec[i], fut_stop = fut_stop_vec[i])[6]
    
    
    
    
  }
  
  plot_h0 = ggplot(plot_df_h0, aes(fill=Test, y=ESS, x=reorder(Stopping,-ESS))) + 
    geom_bar(position="dodge", stat="identity") + 
    ggtitle("ESS changes across different methods of early stopping - H0") + 
    xlab("Stopping type")
  
  print(plot_h0)
  
  plot_h1 = ggplot(plot_df_h1, aes(fill=Test, y=ESS, x=reorder(Stopping,-ESS))) + 
    geom_bar(position="dodge", stat="identity") + 
    ggtitle("ESS changes across different methods of early stopping - H1") + 
    xlab("Stopping type")
  
  print(plot_h1)
  
  
}
if (tab_2){
  
  safe_stop_vec = c(0,1,0,1)
  fut_stop_vec = c(0,0,1,1)
  
  stopping_vec_df = rep(c("No futility or safety","Only safety","Only futility","Both futility and safety"),3)
  test_vec_df = rep(c("Post Prob","Bayes Factor","CI testing"), each = 4)
  null_reject_df = rep(1,12)
  num_toxics_df = rep(1,12)
  num_successes_df = rep(1,12)
  
  plot_df_h0 = data.frame("Stopping" = factor(stopping_vec_df), "Test" = factor(test_vec_df), ESS = rep(1,12), null_reject_df, num_toxics_df, num_successes_df)
  plot_df_h1 = data.frame("Stopping" = factor(stopping_vec_df), "Test" = factor(test_vec_df), ESS = rep(1,12), null_reject_df, num_toxics_df, num_successes_df)
  
  #bop2_ess_vec_h0 = numeric(4)
  #bop2_ess_vec_h1 = numeric(4)
  
  #bf_ess_vec_h0 = numeric(4)
  #bf_ess_vec_h1 = numeric(4)
  
  na_ess_vec_h0 = numeric(4)
  na_ess_vec_h1 = numeric(4)
  
  #bop2_conf_list = list(para_bop2,para_bop2_safe,para_bop2_fut,para_bop2_safe_fut)
  #bf_conf_list = list(para_bf,para_bf_safe,para_bf_fut,para_bf_safe_fut)
  na_conf_list = list(para_na,para_na_safe,para_na_fut,para_na_safe_fut)
  
  for (i in 1:4){
    # paired seed sets the seed default in each one to 2026 
    
    #plot_df_h0$ESS[i] = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 1, conf = bop2_conf_list[[i]], paired_seed = TRUE, safe_stop = safe_stop_vec[i], fut_stop = fut_stop_vec[i])[6]
    #plot_df_h1$ESS[i] = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 1, conf = bop2_conf_list[[i]], paired_seed = TRUE, safe_stop = safe_stop_vec[i], fut_stop = fut_stop_vec[i])[6]
    
    #plot_df_h0$ESS[i+4] = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 3, conf = bf_conf_list[[i]], paired_seed = TRUE, safe_stop = safe_stop_vec[i], fut_stop = fut_stop_vec[i])[6]
    #plot_df_h1$ESS[i+4] = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 3, conf = bf_conf_list[[i]], paired_seed = TRUE, safe_stop = safe_stop_vec[i], fut_stop = fut_stop_vec[i])[6]
    
    temp_h0 = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 4, alpha = na_conf_list[[i]], paired_seed = TRUE, safe_stop = safe_stop_vec[i], fut_stop = fut_stop_vec[i])
    temp_h1 = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = 4, alpha = na_conf_list[[i]], paired_seed = TRUE, safe_stop = safe_stop_vec[i], fut_stop = fut_stop_vec[i])
    
    plot_df_h0$null_reject_df[i+8] = temp_h0[1]
    plot_df_h1$null_reject_df[i+8] = temp_h1[1]
    
    plot_df_h0$num_successes_df[i+8] = temp_h0[7]
    plot_df_h1$num_successes_df[i+8] = temp_h1[7]
    
    plot_df_h0$num_toxics_df[i+8] = temp_h0[8]
    plot_df_h1$num_toxics_df[i+8] = temp_h1[8]
    
    plot_df_h0$ESS[i+8] = temp_h0[6]
    plot_df_h1$ESS[i+8] = temp_h1[6]
    
  }
  
  
  print(plot_df_h0)
  
  
  print(plot_df_h1)
  
  
}

## Fig 3
fig_3 = FALSE
if (fig_3){
  heatmap_plot()
}

## Fig 4
fig_4_opt = FALSE
fig_4 = FALSE

if (fig_4_opt){
  set.seed(2026)
  # RAR = 0 is just equal randomisation
  # RAR = 1 is RITS
  # RAR = 2 is Gaussian RITS
  # RAR = 3 is safety
  # RAR = 4 is just standard efficacy based BRAR
  
  # ER
  print(0)
  
  para_bop2_safe_fut = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, safe_stop = TRUE, fut_stop = TRUE) # BOP2
  para_bf_safe_fut = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, safe_stop = TRUE, fut_stop = TRUE) # Bayes Factor
  para_na_safe_fut = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, safe_stop = TRUE, fut_stop = TRUE) # Normal Approx
  
  print(1)
  
  # RITS
  
  para_bop2_safe_fut_1 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, safe_stop = TRUE, fut_stop = TRUE, RAR = 1) # BOP2
  para_bf_safe_fut_1 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, safe_stop = TRUE, fut_stop = TRUE, RAR = 1) # Bayes Factor
  para_na_safe_fut_1 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, safe_stop = TRUE, fut_stop = TRUE, RAR = 1) # Normal Approx
  
  print(2)
  
  # G RITS
  
  para_bop2_safe_fut_2 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, safe_stop = TRUE, fut_stop = TRUE, RAR = 2) # BOP2
  para_bf_safe_fut_2 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, safe_stop = TRUE, fut_stop = TRUE, RAR = 2) # Bayes Factor
  para_na_safe_fut_2 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, safe_stop = TRUE, fut_stop = TRUE, RAR = 2) # Normal Approx
  
  print(3)
  
  # SAFER
  
  para_bop2_safe_fut_3 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, safe_stop = TRUE, fut_stop = TRUE, RAR = 3) # BOP2
  para_bf_safe_fut_3 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, safe_stop = TRUE, fut_stop = TRUE, RAR = 3) # Bayes Factor
  para_na_safe_fut_3 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, safe_stop = TRUE, fut_stop = TRUE, RAR = 3) # Normal Approx
  
  print(4)
  
  # BRAR
  
  para_bop2_safe_fut_4 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, safe_stop = TRUE, fut_stop = TRUE, RAR = 4) # BOP2
  para_bf_safe_fut_4 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, safe_stop = TRUE, fut_stop = TRUE, RAR = 4) # Bayes Factor
  para_na_safe_fut_4 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, safe_stop = TRUE, fut_stop = TRUE, RAR = 4) # Normal Approx
  
  print(5)
  
}else{
  
  
  # ER
  
  para_bop2_safe_fut = c(0.9937744, 0.9799072, 0.9692627, 0.9625244) # Post prob
  para_bf_safe_fut = c(0.007945557, 0.009997803, 0.109008789, 0.138891602)  # Bayes Factor
  para_na_safe_fut = 0.04262695 # Normal Approx
  
  
  # RITS
  
  para_bop2_safe_fut_1 = c(0.9920166, 0.9812256, 0.9725830, 0.9587158) # Post prob
  para_bf_safe_fut_1 = c(0.006605225, 0.009997803, 0.112963867, 0.163940430)  # Bayes Factor
  para_na_safe_fut_1 = 0.04301758 # Normal Approx
  
  # G RITS
  
  para_bop2_safe_fut_2 = c(0.9920166, 0.9803955, 0.9718018, 0.9561279) # Post prob
  para_bf_safe_fut_2 = c(0.005141846, 0.009997803, 0.114721680, 0.163500977)  # Bayes Factor
  para_na_safe_fut_2 = 0.0425293 # Normal Approx
  
  # SAFER
  
  para_bop2_safe_fut_3 = c(0.9920166, 0.9825928, 0.9686768, 0.9574463) # Post prob
  para_bf_safe_fut_3 = c(0.007607178, 0.009997803, 0.108569336, 0.163940430)  # Bayes Factor
  para_na_safe_fut_3 = 0.04282227 # Normal Approx
  
  # BRAR
  
  para_bop2_safe_fut_4 = c(0.9921143, 0.9773193, 0.9702881, 0.9617432) # Post prob
  para_bf_safe_fut_4 = c(0.005489014, 0.009997803, 0.111645508, 0.158666992)  # Bayes Factor
  para_na_safe_fut_4 = 0.04233398 # Normal Approx
  
}
if (fig_4){
  ## We want ESS and number of patients allocated to the inferior arm (ESS*(1-prop))
  # Note paired_seed = TRUE is setting all the seeds 
  
  # Simulating over 3 tests, 5 RAR types
  
  test_vec = c(1,3,4) # As test = 2 is defunct legacy code
  
  bop2_list = list(para_bop2_safe_fut,para_bop2_safe_fut_1,para_bop2_safe_fut_2,para_bop2_safe_fut_3,para_bop2_safe_fut_4)
  bf_list = list(para_bf_safe_fut,para_bf_safe_fut_1,para_bf_safe_fut_2,para_bf_safe_fut_3,para_bf_safe_fut_4)
  na_list = list(para_na_safe_fut, para_na_safe_fut_1, para_na_safe_fut_2, para_na_safe_fut_3, para_na_safe_fut_4)
  
  conf_list = list(bop2_list, bf_list, na_list) #  for ease of indexing
  
  ess_h0_list = vector(mode = "list", length = 3)
  ess_h1_list = vector(mode = "list", length = 3)
  num_inf_h0_list = vector(mode = "list", length = 3)
  num_inf_h1_list = vector(mode = "list", length = 3)
  
  
  for (i in 1:3){
    # Iterating over tests
    
    temp_ess_h0 = numeric(5)
    temp_ess_h1 = numeric(5)
    
    temp_inf_num_h0 = numeric(5)
    temp_inf_num_h1 = numeric(5)
    
    for (j in 0:4){
      # Iterating over RAR
      print(j)
      temp_conf = (conf_list[[i]])[[j + 1]]
      print(temp_conf)
      
      if (i == 3){
        temp_H0 = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = test_vec[i], alpha = temp_conf, paired_seed = TRUE, safe_stop = TRUE, fut_stop = TRUE, RAR = j)
        
        temp_H1 = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = test_vec[i], alpha = temp_conf, paired_seed = TRUE, safe_stop = TRUE, fut_stop = TRUE, RAR = j)
        
      }
      else{
        temp_H0 = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = test_vec[i], conf = temp_conf, paired_seed = TRUE, safe_stop = TRUE, fut_stop = TRUE, RAR = j)
        
        temp_H1 = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = test_vec[i], conf = temp_conf, paired_seed = TRUE, safe_stop = TRUE, fut_stop = TRUE, RAR = j)
        
      }
      
      
      temp_ess_h0[j+1] = temp_H0[6]
      temp_ess_h1[j+1] = temp_H1[6]
      
      temp_inf_num_h0[j+1] = temp_H0[2] + temp_H0[3]
      temp_inf_num_h1[j+1] = temp_H1[4] + temp_H1[5]
      
      
    }
    
    
    ess_h0_list[[i]] = temp_ess_h0
    ess_h1_list[[i]] = temp_ess_h1
    num_inf_h0_list[[i]] = temp_inf_num_h0
    num_inf_h1_list[[i]] = temp_inf_num_h1
    
    
    
    
    
  }
  
  
  plot(y = ess_h0_list[[1]], x = 0:4, ylim = c(70,180), main = "ESS and patients in the inferior arm under H0", xlab = "RAR type", ylab = "Number of patients", xaxt = "n" )
  lines(y = ess_h0_list[[1]], x = 0:4, col = "red")
  lines(y = ess_h0_list[[2]], x = 0:4, type = "p", pch = 2 )
  lines(y = ess_h0_list[[2]], x = 0:4, col = "red")
  lines(y = ess_h0_list[[3]], x = 0:4, type = "p", pch = 3 )
  lines(y = ess_h0_list[[3]], x = 0:4, col = "red")
  
  lines(y = num_inf_h0_list[[1]], x = 0:4, type = "p", pch = 1 )
  lines(y = num_inf_h0_list[[1]], x = 0:4, col = "blue")
  lines(y = num_inf_h0_list[[2]], x = 0:4, type = "p", pch = 2 )
  lines(y = num_inf_h0_list[[2]], x = 0:4, col = "blue")
  lines(y = num_inf_h0_list[[3]], x = 0:4, type = "p", pch = 3 )
  lines(y = num_inf_h0_list[[3]], x = 0:4, col = "blue")
  
  legend(1, 130, legend=c("ESS", "Patients in the inferior arm"), 
         fill = c("red","blue"))
  legend(0, 130, legend=c("BOP2", "CI test", "Bayes Factor"), 
         pch = c(1,3,2))
  axis(1, at=0:4, labels=c("ER","RITS","Gaussian RITS","SAFER","BRAR"))
  
  plot(y = ess_h1_list[[1]], x = 0:4, ylim = c(100,230), main = "ESS and patients in the inferior arm under H1", xlab = "RAR type", ylab = "Number of patients", xaxt = "n" )
  lines(y = ess_h1_list[[1]], x = 0:4, col = "red")
  lines(y = ess_h1_list[[2]], x = 0:4, type = "p", pch = 2 )
  lines(y = ess_h1_list[[2]], x = 0:4, col = "red")
  lines(y = ess_h1_list[[3]], x = 0:4, type = "p", pch = 3 )
  lines(y = ess_h1_list[[3]], x = 0:4, col = "red")
  
  lines(y = num_inf_h1_list[[1]], x = 0:4, type = "p", pch = 1 )
  lines(y = num_inf_h1_list[[1]], x = 0:4, col = "blue")
  lines(y = num_inf_h1_list[[2]], x = 0:4, type = "p", pch = 2 )
  lines(y = num_inf_h1_list[[2]], x = 0:4, col = "blue")
  lines(y = num_inf_h1_list[[3]], x = 0:4, type = "p", pch = 3 )
  lines(y = num_inf_h1_list[[3]], x = 0:4, col = "blue")
  
  legend(1, 180, legend=c("ESS", "Patients in the inferior arm"), 
         fill = c("red","blue"))
  legend(0, 180, legend=c("BOP2", "CI test", "Bayes Factor"), 
         pch = c(1,3,2))
  axis(1, at=0:4, labels=c("ER","RITS","Gaussian RITS","SAFER","BRAR"))
  
  
  
  #temp = OC_gen()
  
  #ESS[i] = temp[x]
  #Prop[i] = temp[y]
  
  # num_inf_H0 = ESS*prop
  # num_inf_H1 = ESS*(1-prop)
  
}


## Table 3+4
tab_3_4_opt = FALSE
tab_3_4 = FALSE

if (tab_3_4_opt){
  set.seed(2026)
  # RAR = 0 is just equal randomisation
  # RAR = 1 is RITS
  # RAR = 2 is Gaussian RITS
  # RAR = 3 is safety
  # RAR = 4 is just standard efficacy based BRAR
  
  # ER
  print(0)
  
  para_bop2_safe_fut = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, safe_stop = TRUE, fut_stop = TRUE) # BOP2
  para_bf_safe_fut = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, safe_stop = TRUE, fut_stop = TRUE) # Bayes Factor
  para_na_safe_fut = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, safe_stop = TRUE, fut_stop = TRUE) # Normal Approx
  
  print(1)
  
  # RITS
  
  para_bop2_safe_fut_1 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, safe_stop = TRUE, fut_stop = TRUE, RAR = 1) # BOP2
  para_bf_safe_fut_1 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, safe_stop = TRUE, fut_stop = TRUE, RAR = 1) # Bayes Factor
  para_na_safe_fut_1 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, safe_stop = TRUE, fut_stop = TRUE, RAR = 1) # Normal Approx
  
  print(2)
  
  # G RITS
  
  para_bop2_safe_fut_2 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, safe_stop = TRUE, fut_stop = TRUE, RAR = 2) # BOP2
  para_bf_safe_fut_2 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, safe_stop = TRUE, fut_stop = TRUE, RAR = 2) # Bayes Factor
  para_na_safe_fut_2 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, safe_stop = TRUE, fut_stop = TRUE, RAR = 2) # Normal Approx
  
  print(3)
  
  # SAFER
  
  para_bop2_safe_fut_3 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, safe_stop = TRUE, fut_stop = TRUE, RAR = 3) # BOP2
  para_bf_safe_fut_3 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, safe_stop = TRUE, fut_stop = TRUE, RAR = 3) # Bayes Factor
  para_na_safe_fut_3 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, safe_stop = TRUE, fut_stop = TRUE, RAR = 3) # Normal Approx
  
  print(4)
  
  # BRAR
  
  para_bop2_safe_fut_4 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 1, safe_stop = TRUE, fut_stop = TRUE, RAR = 4) # BOP2
  para_bf_safe_fut_4 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 3, safe_stop = TRUE, fut_stop = TRUE, RAR = 4) # Bayes Factor
  para_na_safe_fut_4 = threshold_finder_v2_no_safety(n = 320, IAn = c(80,160,240,320), test = 4, safe_stop = TRUE, fut_stop = TRUE, RAR = 4) # Normal Approx
  
  print(5)
  
}else{
  
  
  # ER
  
  para_bop2_safe_fut = c(0.9937744, 0.9799072, 0.9692627, 0.9625244) # Post prob
  para_bf_safe_fut = c(0.007945557, 0.009997803, 0.109008789, 0.138891602)  # Bayes Factor
  para_na_safe_fut = 0.04262695 # Normal Approx
  
  
  # RITS
  
  para_bop2_safe_fut_1 = c(0.9920166, 0.9812256, 0.9725830, 0.9587158) # Post prob
  para_bf_safe_fut_1 = c(0.006605225, 0.009997803, 0.112963867, 0.163940430)  # Bayes Factor
  para_na_safe_fut_1 = 0.04301758 # Normal Approx
  
  # G RITS
  
  para_bop2_safe_fut_2 = c(0.9920166, 0.9803955, 0.9718018, 0.9561279) # Post prob
  para_bf_safe_fut_2 = c(0.005141846, 0.009997803, 0.114721680, 0.163500977)  # Bayes Factor
  para_na_safe_fut_2 = 0.0425293 # Normal Approx
  
  # SAFER
  
  para_bop2_safe_fut_3 = c(0.9920166, 0.9825928, 0.9686768, 0.9574463) # Post prob
  para_bf_safe_fut_3 = c(0.007607178, 0.009997803, 0.108569336, 0.163940430)  # Bayes Factor
  para_na_safe_fut_3 = 0.04282227 # Normal Approx
  
  # BRAR
  
  para_bop2_safe_fut_4 = c(0.9921143, 0.9773193, 0.9702881, 0.9617432) # Post prob
  para_bf_safe_fut_4 = c(0.005489014, 0.009997803, 0.111645508, 0.158666992)  # Bayes Factor
  para_na_safe_fut_4 = 0.04233398 # Normal Approx
  
}
if (tab_3_4){
  ## We want ESS and number of patients allocated to the inferior arm (ESS*(1-prop))
  # Note paired_seed = TRUE is setting all the seeds 
  
  # Simulating over 3 tests, 5 RAR types
  
  test_vec = c(1,3,4) # As test = 2 is defunct legacy code
  
  bop2_list = list(para_bop2_safe_fut,para_bop2_safe_fut_1,para_bop2_safe_fut_2,para_bop2_safe_fut_3,para_bop2_safe_fut_4)
  bf_list = list(para_bf_safe_fut,para_bf_safe_fut_1,para_bf_safe_fut_2,para_bf_safe_fut_3,para_bf_safe_fut_4)
  na_list = list(para_na_safe_fut, para_na_safe_fut_1, para_na_safe_fut_2, para_na_safe_fut_3, para_na_safe_fut_4)
  
  conf_list = list(bop2_list, bf_list, na_list) #  for ease of indexing
  
  h0_full_list = vector(mode = "list", length = 5)
  h1_full_list  = vector(mode = "list", length = 5)
  
  
  
  for (j in 0:4){
    # Iterating over tests
    
    temp_h0_list = vector(mode = "list", length = 3)
    temp_h1_list = vector(mode = "list", length = 3)
    
    for (i in 1:3){
      # Iterating over RAR
      print(i)
      temp_conf = (conf_list[[i]])[[j + 1]]
      
      if (i == 3){
        temp_H0 = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = test_vec[i], alpha = temp_conf, paired_seed = TRUE, safe_stop = TRUE, fut_stop = TRUE, RAR = j)
        
        temp_H1 = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = test_vec[i], alpha = temp_conf, paired_seed = TRUE, safe_stop = TRUE, fut_stop = TRUE, RAR = j)
        
      }
      else{
        temp_H0 = OC_gen(eff_treat = 0.9, tox_treat = 0.5, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = test_vec[i], conf = temp_conf, paired_seed = TRUE, safe_stop = TRUE, fut_stop = TRUE, RAR = j)
        
        temp_H1 = OC_gen(eff_treat = 0.96, tox_treat = 0.4, eff_con = 0.96, tox_con = 0.5, N = 320, IAn = c(80,160,240,320), test = test_vec[i], conf = temp_conf, paired_seed = TRUE, safe_stop = TRUE, fut_stop = TRUE, RAR = j)
        
      }
      
      
      temp_h0_list[[i]] = temp_H0
      temp_h1_list[[i]] = temp_H1
      
      
      
      
    }
    
    print(j*10)
    
    print("H0 - BOP2")
    to_print = c(temp_h0_list[[1]][1],temp_h0_list[[1]][6],(temp_h0_list[[1]][4] + temp_h0_list[[1]][5])/temp_h0_list[[1]][6],temp_h0_list[[1]][3] + temp_h0_list[[1]][5] )
    print(to_print)
    
    print("H0 - BF")
    to_print = c(temp_h0_list[[2]][1],temp_h0_list[[2]][6],(temp_h0_list[[2]][4] + temp_h0_list[[2]][5])/temp_h0_list[[2]][6],temp_h0_list[[2]][3] + temp_h0_list[[2]][5] )
    print(to_print)
    
    
    print("H0 - NA")
    to_print = c(temp_h0_list[[3]][1],temp_h0_list[[3]][6],(temp_h0_list[[3]][4] + temp_h0_list[[3]][5])/temp_h0_list[[3]][6],temp_h0_list[[3]][3] + temp_h0_list[[3]][5] )
    print(to_print)
    
    
    
    
    print("H1 - BOP2")
    to_print = c(temp_h1_list[[1]][1],temp_h1_list[[1]][6],(temp_h1_list[[1]][2] + temp_h1_list[[1]][3])/temp_h1_list[[1]][6],temp_h1_list[[1]][8] )
    print(to_print)
    
    print("H1 - BF")
    to_print = c(temp_h1_list[[2]][1],temp_h1_list[[2]][6],(temp_h1_list[[2]][2] + temp_h1_list[[2]][3])/temp_h1_list[[2]][6],temp_h1_list[[2]][8] )
    print(to_print)
    
    
    print("H1 - NA")
    to_print = c(temp_h1_list[[3]][1],temp_h1_list[[3]][6],(temp_h1_list[[3]][2] + temp_h1_list[[3]][3])/temp_h1_list[[3]][6],temp_h1_list[[3]][8] )
    print(to_print)
    
    
    h0_full_list[[j+1]] = temp_h0_list
    h1_full_list[[j+1]] = temp_h1_list
    
    
    
    
    
  }
  
  
  
  
  #temp = OC_gen()
  
  #ESS[i] = temp[x]
  #Prop[i] = temp[y]
  
  # num_inf_H0 = ESS*prop
  # num_inf_H1 = ESS*(1-prop)
  
  browser()
  
}


# Table 5
tab_5_opt = FALSE
tab_5 = FALSE

if (tab_5_opt){
  set.seed(2026)
  safe_conf = safety_opt(safe_test = 1)
  safe_conf_2 = safety_opt(safe_test = 2)
  safe_conf_3 = safety_opt(safe_test = 3)
  safe_conf_4 = safety_opt(safe_test = 4)
} else{
  safe_conf = 0.7734375
  safe_conf_2 = 0.4335938
  safe_conf_3 = 0.828125
  safe_conf_4 = 0.6523438
}
if (tab_5){
  
  print("H0 then H1 bop2")
  print(OC_gen_cts(eff_treat = 0.9, tox_treat = 0.15, safe_stop = TRUE, safe_test = 1, safe_conf = safe_conf, paired_seed = TRUE, n_sim = 100000))
  print(OC_gen_cts(eff_treat = 0.96, tox_treat = 0.05, safe_stop = TRUE, safe_test = 1, safe_conf = safe_conf, paired_seed = TRUE, n_sim = 100000))
  
  print("H0 then H1 cox ph")
  print(OC_gen_cts(eff_treat = 0.9, tox_treat = 0.15, safe_stop = TRUE, safe_test = 2, safe_conf = safe_conf_2, paired_seed = TRUE, n_sim = 100000))
  print(OC_gen_cts(eff_treat = 0.96, tox_treat = 0.05, safe_stop = TRUE, safe_test = 2, safe_conf = safe_conf_2, paired_seed = TRUE, n_sim = 100000))
  
  print("H0 then H1 AJE")
  print(OC_gen_cts(eff_treat = 0.9, tox_treat = 0.15, safe_stop = TRUE, safe_test = 3, safe_conf = safe_conf_3, paired_seed = TRUE, n_sim = 100000))
  print(OC_gen_cts(eff_treat = 0.96, tox_treat = 0.05, safe_stop = TRUE, safe_test = 3, safe_conf = safe_conf_3, paired_seed = TRUE, n_sim = 100000))
  
  print("H0 then H1 RMST")
  print(OC_gen_cts(eff_treat = 0.9, tox_treat = 0.15, safe_stop = TRUE, safe_test = 4, safe_conf = safe_conf_4, paired_seed = TRUE, n_sim = 100000))
  print(OC_gen_cts(eff_treat = 0.96, tox_treat = 0.05, safe_stop = TRUE, safe_test = 4, safe_conf = safe_conf_4, paired_seed = TRUE, n_sim = 100000))
  
}


# Tables 6+7
tab_6_7_opt = FALSE
tab_6_7 = FALSE

if (tab_6_7_opt){
  set.seed(2026)
  
  conf_compete = threshold_finder_v2_no_safety(eff_con = 0.53, safe_stop = FALSE, alpha = 0.07, fut_stop = FALSE, test = 4, RAR = 2, n = 320, IAn = c(80,160,240,320))
  
  safe_conf = safety_opt(safe_test = 1, fut_stop = FALSE, eff_con = 0.53, alpha = conf_compete)
  safe_conf_2 = safety_opt(safe_test = 2, fut_stop = FALSE, eff_con = 0.53, alpha = conf_compete)
  safe_conf_3 = safety_opt(safe_test = 3, fut_stop = FALSE, eff_con = 0.53, alpha = conf_compete)
  safe_conf_4 = safety_opt(safe_test = 4, fut_stop = FALSE, eff_con = 0.53, alpha = conf_compete)
  
} else{
  conf_compete = 0.05698242
  
  safe_conf = 0.7480469
  safe_conf_2 = 0.4921875
  safe_conf_3 = 0.609375
  safe_conf_4 = 0.6171875
}
if (tab_6_7){
  
  print("PH then competing risks bop2")
  print(OC_gen_cts(eff_treat = 0.9, tox_treat = 0.15, safe_stop = TRUE, safe_test = 1, safe_conf = safe_conf, paired_seed = TRUE, n_sim = 10000, PH_true = FALSE ))
  print(OC_gen_cts(eff_treat = 0.47, tox_treat = 0.15, eff_con = 0.53, safe_stop = TRUE, safe_test = 1, safe_conf = safe_conf, paired_seed = TRUE, n_sim = 10000, fut_stop = FALSE, alpha = conf_compete))
  
  print("PH then competing risks cox ph")
  print(OC_gen_cts(eff_treat = 0.9, tox_treat = 0.15, safe_stop = TRUE, safe_test = 2, safe_conf = safe_conf_2, paired_seed = TRUE, n_sim = 10000, PH_true = FALSE))
  print(OC_gen_cts(eff_treat = 0.47, tox_treat = 0.15, eff_con = 0.53, safe_stop = TRUE, safe_test = 2, safe_conf = safe_conf_2, paired_seed = TRUE, n_sim = 10000, fut_stop = FALSE, alpha = conf_compete))
  
  print("PH then competing risks AJE")
  print(OC_gen_cts(eff_treat = 0.9, tox_treat = 0.15, safe_stop = TRUE, safe_test = 3, safe_conf = safe_conf_3, paired_seed = TRUE, n_sim = 10000, PH_true = FALSE))
  print(OC_gen_cts(eff_treat = 0.47, tox_treat = 0.15, eff_con = 0.53, safe_stop = TRUE, safe_test = 3, safe_conf = safe_conf_3, paired_seed = TRUE, n_sim = 10000, fut_stop = FALSE, alpha = conf_compete))
  
  print("PH then competing risks RMST")
  print(OC_gen_cts(eff_treat = 0.9, tox_treat = 0.15, safe_stop = TRUE, safe_test = 4, safe_conf = safe_conf_4, paired_seed = TRUE, n_sim = 10000, PH_true = FALSE))
  print(OC_gen_cts(eff_treat = 0.47, tox_treat = 0.15, eff_con = 0.53, safe_stop = TRUE, safe_test = 4, safe_conf = safe_conf_4, paired_seed = TRUE, n_sim = 10000, fut_stop = FALSE, alpha = conf_compete))
  
}



# Table 8
tab_8_opt = FALSE
tab_8 = FALSE

if (tab_8_opt){
  set.seed(2026)
  
  conf_tab_8 = threshold_finder_v2_no_safety(eff_con = 0.96, safe_stop = FALSE, alpha = 0.07, fut_stop = TRUE, test = 4, RAR = 2, n = 320, IAn = c(80,160,240,320))
  
  safe_conf = safety_opt(safe_test = 1, fut_stop = TRUE, binary = TRUE, alpha = conf_tab_8)
  safe_conf_cts = safety_opt(safe_test = 1, fut_stop = TRUE, alpha = conf_tab_8)
} else{
  
  conf_tab_8 = 0.04282227
  
  safe_conf = 0.2207031
  safe_conf_cts = 0.7753906
}
if (tab_8){
  
  print("H0 then H1 binary")
  print(OC_gen(eff_treat = 0.9, tox_treat = 0.15, eff_con = 0.96, tox_con = 0.15, N = 320, IAn = c(80,160,240,320), safe_stop = TRUE, safe_conf = rep(safe_conf, 4), paired_seed = TRUE, n_sim = 100000, alpha = conf_tab_8, fut_stop = TRUE, RAR = 2, test = 4))
  print(OC_gen(eff_treat = 0.96, tox_treat = 0.05, eff_con = 0.96, tox_con = 0.15, N = 320, IAn = c(80,160,240,320), safe_stop = TRUE, safe_conf = rep(safe_conf, 4), paired_seed = TRUE, n_sim = 100000, alpha = conf_tab_8, fut_stop = TRUE, RAR = 2, test = 4))
  
  
  
  
  
  print("H0 then H1 cts")
  print(OC_gen_cts(eff_treat = 0.9, tox_treat = 0.15, eff_con = 0.96, tox_con = 0.15, n = 320, IAn = c(80,160,240,320), safe_stop = TRUE, safe_test = 1, safe_conf = safe_conf_cts, paired_seed = TRUE, n_sim = 100000, alpha = conf_tab_8, fut_stop = TRUE))
  print(OC_gen_cts(eff_treat = 0.96, tox_treat = 0.05, eff_con = 0.96, tox_con = 0.15, n = 320, IAn = c(80,160,240,320), safe_stop = TRUE, safe_test = 1, safe_conf = safe_conf_cts, paired_seed = TRUE, n_sim = 100000, alpha = conf_tab_8, fut_stop = TRUE))
  
  browser()
}




# Appendix Figures 5+6

fig_5_6 = FALSE

if (fig_5_6){
  
  set.seed(2026)
  
  n = 160 # per arm
  
  
  q_vec = seq(1, 99, 1)/100
  
  post_prob_real = numeric(length(q_vec))
  post_prob_normal = numeric(length(q_vec))
  post_prob_sim = numeric(length(q_vec))
  
  for(i in 1:length(q_vec)){
    
    post_prob_real[i] = qq_inv_function(q_vec[i], method = 1, n = n)
    
    
    post_prob_normal[i] = qq_inv_function(q_vec[i], method = 2, n = n)
    
    
    post_prob_sim[i] = qq_inv_function(q_vec[i], method = 3, n = n )
    
    
  }
  
  plot(x = post_prob_real, y = post_prob_normal, main = "QQ plot - Normal approximation vs True", xlab = "Theoretical quantities", ylab = "Normal approximated quantities")
  lines(y = c(0,1000), x = c(0,1000), col = "red")
  legend(x = post_prob_real[1], y = post_prob_normal[length(q_vec)], legend = c("Y = X"), col = "red", lty = 1  )
  
  plot(x = post_prob_real, y = post_prob_sim, main = "QQ plot - Simulation approximation vs True", xlab = "Theoretical quantities", ylab = "Simulated approximated quantities")
  lines(y = c(0,1000), x = c(0,1000), col = "red")
  legend(x = post_prob_real[1], y = post_prob_sim[length(q_vec)], legend = c("Y = X"), col = "red", lty = 1  )
  
  
  
  
  
}


# Appendix Figure 7

fig_7 = FALSE

if (fig_7){
  
  set.seed(2026)
  
  eff_con_vec = seq(0.12,0.96,0.06)
  
  ss_vec = numeric(length(eff_con_vec))
  
  for (i in 1:length(eff_con_vec)){
    
    print(i)
    
    ss_vec[i] = ss_finder(0.4,eff_con_vec[i],0.5, test = 4, fut_stop = FALSE, safe_stop = FALSE, RAR = 2, power = 0.8, n_sim = 10000)
    
  }
  
  plot(y = ss_vec, x = eff_con_vec, lty = 1, xlab = "Control efficacy", ylab = "Sample size", main = "Plot of sample size required for 80% power with varying efficacy")
  lines(y = ss_vec, x = eff_con_vec) 
  
}