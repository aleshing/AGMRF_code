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
load(file = "../Results/smoothed_direct_conflict_int.RData")
load(file = "../Results/smoothed_direct_adaptive_conflict_int.RData")
source("../../helper_files/calc_theta_prior.R")

#### Compare prior and posterior for theta ####
pc.u.theta <- 0.75
pc.alpha.theta <- 0.75
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
inla.mmarginal(data.frame(x = theta_prior$theta, y = exp(theta_prior$log_density)))
theta_post <- 
  inla.tmarginal(fun = SUMMER::expit, 
                 marginal = fit_adaptive_conflict_int$marginals.hyperpar$`Theta3 for region.struct`)
theta_dens <- rbind(data.frame(theta = theta_prior$theta, 
                               density = exp(theta_prior$log_density), 
                               Which = "Prior"),
                    data.frame(theta = theta_post[, 1], 
                               density = theta_post[, 2], 
                               Which = "Post"))
theta_dens %>% ggplot(aes(x = theta, y = density, colour = Which)) + 
  geom_line() + xlab("Theta") + ylab("Density") 
ggsave(filename = "../Results/theta_prior_post_comp_conflict_int.pdf",
       height = 5, width = 8)

#### Parameter Summaries ####
smoothed_summaries <- 
  data.frame(Period = c("Conflict", "Non-Conflict", rep(NA, 3)),
             Parameter = c(rep("mu", 2), "tau", "phi", "theta"),
             Mean = NA, SD = NA, `2.5% Quantile` = NA, Median = NA, `97.5% Quantile` = NA, Mode = NA)
smoothed_summaries[1:2, 2:7] <- 
  fit_conflict_int$summary.fixed[, 1:6] %>% unname()
smoothed_summaries[3:4, 2:7] <- 
  fit_conflict_int$summary.hyperpar[1:2, 1:6] %>% unname() 


tau_post <- 
  inla.tmarginal(fun = exp, 
                 marginal = fit_adaptive_conflict_int$marginals.hyperpar$`Theta1 for region.struct`)
tau_sum <- inla.zmarginal(tau_post)
tau_mode <- inla.mmarginal(tau_post)
phi_post <- 
  inla.tmarginal(fun = SUMMER::expit, 
                 marginal = fit_adaptive_conflict_int$marginals.hyperpar$`Theta2 for region.struct`)
phi_sum <- inla.zmarginal(phi_post)
phi_mode <- inla.mmarginal(phi_post)
theta_sum <- inla.zmarginal(theta_post)
theta_mode <- inla.mmarginal(theta_post)

proposed_summaries <- 
  data.frame(Period = c("Conflict", "Non-Conflict", rep(NA, 3)),
             Parameter = c(rep("mu", 2), "tau", "phi", "theta"),
             Mean = NA, SD = NA, `2.5% Quantile` = NA, Median = NA, `97.5% Quantile` = NA, Mode = NA)
proposed_summaries[1:2, 2:7] <- 
  fit_adaptive_conflict_int$summary.fixed[, 1:6] %>% unname()
proposed_summaries[3, 2:7] <- 
  c(tau_sum[c(1:3, 5, 7)], tau_mode)
proposed_summaries[4, 2:7] <- 
  c(phi_sum[c(1:3, 5, 7)], phi_mode)
proposed_summaries[5, 2:7] <- 
  c(theta_sum[c(1:3, 5, 7)], theta_mode)

summaries <- rbind(cbind(Model = "Smoothed Direct Conflict-Intercept", smoothed_summaries),
                   cbind(Model = "Proposed Conflict-Intercept", proposed_summaries))

print(xtable(summaries), include.rownames = FALSE)

#### Compare smoothed direct estimates ####
out_combined <- rbind(cbind(out_conflict_int, Model = "Smoothed Direct Conflict-Intercept"),
                      cbind(out_adaptive_conflict_int, Model = "Proposed Conflict-Intercept"))
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
  mutate(years = if_else(Model == "Smoothed Direct Conflict-Intercept", years - 0.25, years),
         years = if_else(Model == "Proposed Conflict-Intercept", years + 0.25, years))

plot_direct <- 
  ggplot(direct_estimates, 
         aes(x = as.numeric(years), y = mean, color = surveyYears)) + 
  geom_point() + geom_line() + ylab("U5MR") + xlab("Year") + 
  labs(color = "Survey Year")
ggsave(filename = "../Results/direct_estimates.pdf", height = 5, width = 8)

igme_meta <- igme %>% select(years, est, lower, upper) %>%
  mutate(source = "IGME") %>%
  rbind(meta_analysis_estimates %>%
          mutate(lower = NA, upper = NA, source = "Survey Meta"))

plot_igme <- 
  ggplot(igme_meta, 
         aes(x = as.numeric(years), y = est, color = source)) + 
  ylab("U5MR") + xlab("Year") + 
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  geom_point() + 
  labs(color = "Estimate") 
ggsave(filename = "../Results/igme_meta.pdf", height = 5, width = 8)

ggarrange(plot_direct, plot_igme,
          ncol = 1, nrow = 2, align = "v")
ggsave(filename = "../Results/direct_igme.pdf", height = 10, width = 8)

ggplot(direct_estimates, 
       aes(x = as.numeric(years), y = mean, color = surveyYears)) + 
  geom_point(alpha = 0.5) + geom_line(alpha = 0.5) + ylab("U5MR") + 
  xlab("Year") + 
  labs(color = "Survey Year") +
  annotate("point",  x = igme$years, y = igme$est, col = "blue") + 
  annotate("linerange", x = igme$years, ymin = igme$lower, ymax = igme$upper, 
           col = "blue") +
  annotate("point", x = as.numeric(meta_analysis_estimates$years), 
           y = meta_analysis_estimates$est, col = "red")
ggsave(filename = "../Results/direct_estimates_plus.pdf", height = 5, width = 8)

out_combined_for_plot %>%
  ggplot(aes(x = years, y = Estimate, color = Model)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) + 
  geom_point() + 
  ylab("U5MR") + xlab("Year") + 
  scale_color_manual(values = c("Smoothed Direct Conflict-Intercept" = "grey40", 
                                "Proposed Conflict-Intercept" = "black", 
                                "IGME" = "blue", "Survey Meta" = "red")) 
ggsave(filename = "../Results/smoothed_direct_comp_1_rwanda_conflict_int.pdf", height = 5, width = 8)




