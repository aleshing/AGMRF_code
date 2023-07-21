library(dplyr)
library(INLA)
library(tidyverse)

#### Set simulation parameters ####
num_reps <- 200
num_years <- 30
num_par_settings <- 6
num_prec_settings <- 3

# NA indicates the usage of the triangle trend
par_settings <- data.frame(mu_1 = c(rep(-2, num_par_settings - 2), NA, NA),
                           mu_2 = c(c(-2, -2, -1, -1), NA, NA),
                           tau_1 = rep(20, num_par_settings),
                           tau_2 = c(20, 10, 20, 10, 20, 10))
triangle <- c(rep(-2, 8), 
              seq(from = -2, to = -0.75, length.out = 5)[2:4], 
              -0.75, 
              seq(from = -2, to = -0.75, length.out = 5)[4:2], 
              rep(-2, 15))
prec_settings <- c(150, 300, 75)
trend_settings <- c("Trend: Constant", "Trend: Constant", "Trend: Level Change",
                    "Trend: Level Change", "Trend: Triangle", "Trend: Triangle")
re_settings <- c("Random Effect: Equal Precisions", 
                 "Random Effect: Unequal Precisions",
                 "Random Effect: Equal Precisions", 
                 "Random Effect: Unequal Precisions",
                 "Random Effect: Equal Precisions", 
                 "Random Effect: Unequal Precisions")
par_settings_descrip <- 
  c("Trend: Constant; Random Effect: Equal Precisions",
    "Trend: Constant; Random Effect: Unequal Precisions",
    "Trend: Level Change; Random Effect: Equal Precisions",
    "Trend: Level Change; Random Effect: Unequal Precisions",
    "Trend: Triangle; Random Effect: Equal Precisions",
    "Trend: Triangle; Random Effect: Unequal Precisions")

results <- data.frame(rep = NA, prec_setting = NA, par_setting = NA, 
                      trend_settings = NA, re_settings = NA, model = NA,
                      rmse_logit_risk = NA, rmse_logit_obs = NA, rmse_risk = NA, 
                      rmse_obs = NA, cpo = NA, dic = NA)

for(i in 1:num_prec_settings){
    for(j in 1:num_par_settings){
        results_proposed <- data.frame(rep = 1:num_reps, prec_setting = i, 
                                       par_setting = par_settings_descrip[j], 
                                       trend_settings = trend_settings[j], 
                                       re_settings = re_settings[j],
                                       model = "Proposed",
                                       rmse_logit_risk = 0, 
                                       rmse_logit_obs = 0,
                                       rmse_risk = 0, 
                                       rmse_obs = 0, cpo = 0, dic = 0)
        results_smoothed_direct <- 
          data.frame(rep = 1:num_reps, 
                     prec_setting = i, 
                     par_setting = par_settings_descrip[j], 
                     trend_settings = trend_settings[j], 
                     re_settings = re_settings[j],
                     model = "Smoothed Direct",
                     rmse_logit_risk = 0, 
                     rmse_logit_obs = 0,
                     rmse_risk = 0, 
                     rmse_obs = 0, cpo = 0, dic = 0)
        for(k in 1:num_reps){
            load(paste0("../Results/Fits_proposed/fit_prec_", i, "_par_", j, 
                        "_rep_", k))
            load(paste0("../Results/Fits_smoothed_direct/fit_prec_", i, "_par_", 
                        j, "_rep_", k))
            load(paste0("../Data/dat_prec_", i, "_par_", j, "_rep_", k))
            
            results_proposed$rmse_logit_risk[k] <- 
                sqrt(mean((fit$summary.fitted.values$mean[1:num_years] - 
                               dat$risk) ^ 2))
            results_smoothed_direct$rmse_logit_risk[k] <- 
                sqrt(mean((fit_bym2$summary.fitted.values$mean[1:num_years] - 
                               dat$risk) ^ 2))
            results_proposed$rmse_logit_obs[k] <- 
                sqrt(mean((fit$summary.fitted.values$mean[1:num_years] - 
                               dat$logit.est) ^ 2))
            results_smoothed_direct$rmse_logit_obs[k] <- 
                sqrt(mean((fit_bym2$summary.fitted.values$mean[1:num_years] - 
                               dat$logit.est) ^ 2))
            
            results_proposed$cpo[k] <- -mean(log(fit$cpo$cpo[1:num_years]))
            results_smoothed_direct$cpo[k] <- 
                -mean(log(fit_bym2$cpo$cpo[1:num_years]))
            
            results_proposed$dic[k] <- fit$dic$dic
            results_smoothed_direct$dic[k] <- fit_bym2$dic$dic
            
           
           
        }
        print(paste0("mu_1 = ", par_settings$mu_1[j], 
                     ", mu_2 = ", par_settings$mu_2[j], 
                     ", tau_1 = ",  par_settings$tau_1[j], 
                     ", tau_2 = ",  par_settings$tau_2[j]))
        print("Proposed:")
        print(colMeans(results_proposed[reps, 7:12]))
        print("Smoothed Direct:")
        print(colMeans(results_smoothed_direct[reps, 7:12]))
        
        results <- rbind(results, 
                         results_proposed[reps,], 
                         results_smoothed_direct[reps,])
    }
}

results <- results[-1, ]
save(results, file = paste0("../Results/model_comp"))
