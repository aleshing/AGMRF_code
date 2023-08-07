library(SUMMER)
library(readstata13)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(xtable)

#### Load results ####
load(file = "../Data/generated_data/direct_estimates.RData")
load(file = "../Data/generated_data/structure_matrices.RData")
load(file = "../Results/smoothed_direct_estimates.RData")
load(file = "../Results/smoothed_direct_linear.RData")
load(file = "../Results/smoothed_direct_adaptive_linear.RData")
source("../../helper_files/calc_theta_prior.R")

#### Compare prior and posterior for theta ####
pc.u.theta <- 0.75
pc.alpha.theta <- 0.75
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
inla.mmarginal(data.frame(x = theta_prior$theta, 
                          y = exp(theta_prior$log_density)))
theta_post <- 
  inla.tmarginal(fun = SUMMER::expit, 
                 marginal = fit_adaptive_linear$marginals.hyperpar$`Theta3 for region.struct`)
theta_dens <- rbind(data.frame(theta = theta_prior$theta, 
                               density = exp(theta_prior$log_density), 
                               Which = "Prior"),
                    data.frame(theta = theta_post[, 1], 
                               density = theta_post[, 2], 
                               Which = "Post"))
theta_dens %>% ggplot(aes(x = theta, y = density, colour = Which)) + 
  geom_line() + xlab("Theta") + ylab("Density") 
ggsave(filename = "../Results/theta_prior_post_comp_linear.pdf",
       height = 5, width = 8)

#### Parameter Summaries ####
survey_years <- c(1992, 2000, 2005, 2008, 2010, 2015)

smoothed_summaries <- 
  data.frame(Parameter = c("mu", "beta", "tau", "phi", "theta",
                           paste0("nu_", survey_years)),
             Mean = NA, SD = NA, `2.5% Quantile` = NA, Median = NA, 
             `97.5% Quantile` = NA, Mode = NA)
smoothed_summaries[1:2, 3:8] <- 
  fit_linear$summary.fixed[, 1:6] %>% unname()
smoothed_summaries[3:4, 3:8] <- 
  fit_linear$summary.hyperpar[1:2, 1:6] %>% unname() 
smoothed_summaries[6:11, 3:8] <- 
  fit_linear$summary.random$survey.id[, 2:7] %>% unname() 


tau_post <- 
  inla.tmarginal(fun = exp, 
                 marginal = fit_adaptive_linear$marginals.hyperpar$`Theta1 for region.struct`)
tau_sum <- inla.zmarginal(tau_post)
tau_mode <- inla.mmarginal(tau_post)
phi_post <- 
  inla.tmarginal(fun = SUMMER::expit, 
                 marginal = fit_adaptive_linear$marginals.hyperpar$`Theta2 for region.struct`)
phi_sum <- inla.zmarginal(phi_post)
phi_mode <- inla.mmarginal(phi_post)
theta_sum <- inla.zmarginal(theta_post)
theta_mode <- inla.mmarginal(theta_post)

proposed_summaries <- 
  data.frame(Parameter = c("mu", "beta", "tau", "phi", "theta",
                           paste0("nu_", survey_years)),
             Mean = NA, SD = NA, `2.5% Quantile` = NA, Median = NA, 
             `97.5% Quantile` = NA, Mode = NA)
proposed_summaries[1:2, 3:8] <- 
  fit_adaptive_linear$summary.fixed[, 1:6] %>% unname()
proposed_summaries[3, 3:8] <- 
  c(tau_sum[c(1:3, 5, 7)], tau_mode)
proposed_summaries[4, 3:8] <- 
  c(phi_sum[c(1:3, 5, 7)], phi_mode)
proposed_summaries[5, 3:8] <- 
  c(theta_sum[c(1:3, 5, 7)], theta_mode)
proposed_summaries[6:11, 3:8] <- 
  fit_adaptive_linear$summary.random$survey.id[, 2:7] %>% unname() 

summaries <- rbind(cbind(Model = "Smoothed Direct Linear", 
                         smoothed_summaries),
                   cbind(Model = "Proposed Linear", 
                         proposed_summaries))

print(xtable(summaries), include.rownames = FALSE)

#### Compare smoothed direct estimates ####
out_combined <- rbind(cbind(out_linear, 
                            Model = "Smoothed Direct Linear"),
                      cbind(out_adaptive_linear, 
                            Model = "Proposed Linear"))
out_combined <- out_combined %>% 
  mutate(years = region, years = as.numeric(years),
         Estimate = median)
out_combined_for_plot <- out_combined %>% 
  select(years, Estimate, upper, lower, Model) %>% 
  rbind(data.frame(years = igme$years, Estimate = igme$est, upper = igme$upper,
                   lower = igme$lower, Model = "IGME"),
        data.frame(years = as.numeric(meta_analysis_estimates$years), 
                   Estimate = meta_analysis_estimates$est, upper = NA,
                   lower = NA, Model = "Survey Meta")) %>%
  mutate(years = if_else(Model == "Smoothed Direct Linear", 
                         years - 0.25, years),
         years = if_else(Model == "Proposed Linear", 
                         years + 0.25, years))

out_combined_for_plot %>%
  ggplot(aes(x = years, y = Estimate, color = Model)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) + 
  geom_point() + 
  ylab("U5MR") + xlab("Year") + 
  scale_color_manual(values = c("Smoothed Direct Linear" = "grey40", 
                                "Proposed Linear" = "black", 
                                "IGME" = "blue", "Survey Meta" = "red")) 
ggsave(filename = "../Results/smoothed_direct_comp_1_rwanda_linear.pdf", 
       height = 5, width = 8)




