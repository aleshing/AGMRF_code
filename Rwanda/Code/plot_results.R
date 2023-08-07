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
source("../../helper_files/calc_theta_prior.R")

#### Plot logit precisions ####

direct_estimates %>% 
  ggplot(aes(x = as.numeric(years), y = logit.prec, col = surveyYears)) +
  geom_point() + geom_line() + xlab("Year") + 
  ylab(expression(hat(V)["is"]^-1)) + labs(col = "Survey Year")
ggsave(filename = "../Results/survey_precisions.pdf", height = 5, width = 8)

#### Compare different priors ####
pc.u.theta <- 0.75
pc.alpha.theta <- 0.25
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
theta_dens <- data.frame(Theta = theta_prior$theta, 
                         Density = exp(theta_prior$log_density), 
                         U = 0.75, alpha = 0.25)
pc.u.theta <- 0.75
pc.alpha.theta <- 0.5
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
theta_dens <- rbind(theta_dens,
                    data.frame(Theta = theta_prior$theta, 
                         Density = exp(theta_prior$log_density), 
                         U = 0.75, alpha = 0.5))

pc.u.theta <- 0.75
pc.alpha.theta <- 0.75
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
theta_dens <- rbind(theta_dens,
                    data.frame(Theta = theta_prior$theta, 
                               Density = exp(theta_prior$log_density), 
                               U = 0.75, alpha = 0.75))

pc.u.theta <- 0.75
pc.alpha.theta <- 0.9
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
theta_dens <- rbind(theta_dens,
                    data.frame(Theta = theta_prior$theta, 
                               Density = exp(theta_prior$log_density), 
                               U = 0.75, alpha = 0.9))


pc.u.theta <- 0.6
pc.alpha.theta <- 0.25
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
theta_dens <- rbind(theta_dens,
                    data.frame(Theta = theta_prior$theta, 
                         Density = exp(theta_prior$log_density), 
                         U = 0.6, alpha = 0.25))
pc.u.theta <- 0.6
pc.alpha.theta <- 0.5
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
theta_dens <- rbind(theta_dens,
                    data.frame(Theta = theta_prior$theta, 
                               Density = exp(theta_prior$log_density), 
                               U = 0.6, alpha = 0.5))

pc.u.theta <- 0.6
pc.alpha.theta <- 0.75
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
theta_dens <- rbind(theta_dens,
                    data.frame(Theta = theta_prior$theta, 
                               Density = exp(theta_prior$log_density), 
                               U = 0.6, alpha = 0.75))

pc.u.theta <- 0.6
pc.alpha.theta <- 0.9
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
theta_dens <- rbind(theta_dens,
                    data.frame(Theta = theta_prior$theta, 
                               Density = exp(theta_prior$log_density), 
                               U = 0.6, alpha = 0.9))

theta_dens %>% 
  mutate(alpha = as.character(alpha),
         U = paste("U =", U)) %>%
  filter(alpha != "0.9") %>%
  ggplot(aes(x = Theta, y = Density, colour = alpha)) + 
  geom_line() + xlab("Theta") + ylab("Density")  + facet_wrap(facets = vars(U))

ggsave(filename = "../Results/priors.pdf", height = 5, width = 8)

#### Compare prior and posterior for theta ####
pc.u.theta <- 0.75
pc.alpha.theta <- 0.75
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
inla.mmarginal(data.frame(x = theta_prior$theta, 
                          y = exp(theta_prior$log_density)))
theta_post <- 
    inla.tmarginal(fun = SUMMER::expit, 
                   marginal = smoothed_direct_adaptive_bym2$fit$marginals.hyperpar$`Theta3 for region.struct`)
theta_dens <- rbind(data.frame(theta = theta_prior$theta, 
                               density = exp(theta_prior$log_density), 
                               Which = "Prior"),
                    data.frame(theta = theta_post[, 1], 
                               density = theta_post[, 2], 
                               Which = "Post"))
theta_dens %>% ggplot(aes(x = theta, y = density, colour = Which)) + 
  geom_line() + xlab("Theta") + ylab("Density") 
ggsave(filename = "../Results/theta_prior_post_comp.pdf",
       height = 5, width = 8)

#### Parameter Summaries ####
survey_years <- c(1992, 2000, 2005, 2008, 2010, 2015)

smoothed_summaries <- 
  data.frame(Parameter = c("mu", "beta", "tau", "phi", "theta",
                           paste0("nu_", survey_years)),
             Mean = NA, SD = NA, `2.5% Quantile` = NA, Median = NA, 
             `97.5% Quantile` = NA, Mode = NA)
smoothed_summaries[1:2, 2:7] <- 
  smoothed_direct_bym2$fit$summary.fixed[, 1:6] %>% unname()
smoothed_summaries[3:4, 2:7] <- 
  smoothed_direct_bym2$fit$summary.hyperpar[1:2, 1:6] %>% unname() 
smoothed_summaries[6:11, 2:7] <- 
  smoothed_direct_bym2$fit$summary.random$survey.id[, 2:7] %>% unname() 


tau_post <- 
  inla.tmarginal(fun = exp, 
                 marginal = smoothed_direct_adaptive_bym2$fit$marginals.hyperpar$`Theta1 for region.struct`)
tau_sum <- inla.zmarginal(tau_post)
tau_mode <- inla.mmarginal(tau_post)
phi_post <- 
  inla.tmarginal(fun = SUMMER::expit, 
                 marginal = smoothed_direct_adaptive_bym2$fit$marginals.hyperpar$`Theta2 for region.struct`)
phi_sum <- inla.zmarginal(phi_post)
phi_mode <- inla.mmarginal(phi_post)
theta_sum <- inla.zmarginal(theta_post)
theta_mode <- inla.mmarginal(theta_post)

proposed_summaries <- 
  data.frame(Parameter = c("mu", "beta", "tau", "phi", "theta",
                           paste0("nu_", survey_years)),
             Mean = NA, SD = NA, `2.5% Quantile` = NA, Median = NA, 
             `97.5% Quantile` = NA, Mode = NA)
proposed_summaries[1:2, 2:7] <- 
  smoothed_direct_adaptive_bym2$fit$summary.fixed[, 1:6] %>% unname()
proposed_summaries[3, 2:7] <- 
  c(tau_sum[c(1:3, 5, 7)], tau_mode)
proposed_summaries[4, 2:7] <- 
  c(phi_sum[c(1:3, 5, 7)], phi_mode)
proposed_summaries[5, 2:7] <- 
  c(theta_sum[c(1:3, 5, 7)], theta_mode)
proposed_summaries[6:11, 2:7] <- 
  smoothed_direct_adaptive_bym2$fit$summary.random$survey.id[, 2:7] %>% unname() 

summaries <- rbind(cbind(Model = "Smoothed Direct", smoothed_summaries),
                   cbind(Model = "Proposed", proposed_summaries))

print(xtable(summaries[, 1:7]), include.rownames = FALSE)

#### Compare smoothed direct estimates ####
out_combined <- out_combined %>% mutate(years = as.numeric(years))
out_bym2 <- out_combined %>% filter(prior == "bym2")
out_adaptive_bym2 <- out_combined %>% filter(prior == "adaptive bym2")
out_combined_for_plot <- out_combined %>% 
    mutate(est = median, model = prior, years = as.numeric(region)) %>%
    select(years, est, upper, lower, model) %>% 
    rbind(data.frame(years = igme$years, est = igme$est, upper = igme$upper,
                     lower = igme$lower, model = "IGME"),
          data.frame(years = as.numeric(meta_analysis_estimates$years), 
                     est = meta_analysis_estimates$est, upper = NA,
                     lower = NA, model = "Survey Meta")) %>%
    mutate(years = if_else(model == "bym2", years - 0.25, years),
           years = if_else(model == "adaptive bym2", years + 0.25, years),)

plot_direct <- 
  ggplot(direct_estimates, 
         aes(x = as.numeric(years), y = mean, color = surveyYears)) + 
  geom_point() + geom_line() + ylab("U5MR") + xlab("Year") + 
  labs(color = "Survey Year") + ylim(0.015, 0.325)
ggsave(filename = "../Results/direct_estimates.pdf", height = 5, width = 8)

plot_direct_lines <- 
  ggplot(direct_estimates, 
         aes(x = as.numeric(years), y = mean, color = surveyYears)) + 
  geom_point() + geom_line() + ylab("U5MR") + xlab("Year") + 
  labs(color = "Survey Year") +
  annotate("rect", xmin = 1992.5, xmax = 1999.5, ymin = -Inf, ymax = Inf,
           alpha = 0.25) 
ggsave(filename = "../Results/direct_estimates_lines.pdf", 
       height = 5, width = 8)

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
  labs(color = "Estimate") + ylim(0.015, 0.325)
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
  mutate(Model = model, 
         Model = ifelse(Model == "bym2", "Smoothed Direct", Model),
         Model = ifelse(Model == "adaptive bym2", "Proposed", Model)) %>%
ggplot(aes(x = years, y = est, color = Model)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) + 
    geom_point() + 
    ylab("U5MR") + xlab("Year") + 
    scale_color_manual(values = c("Smoothed Direct" = "grey40", 
                                  "Proposed" = "black", "IGME" = "blue", 
                                  "Survey Meta" = "red")) 
ggsave(filename = "../Results/smoothed_direct_comp_1_rwanda.pdf", 
       height = 5, width = 8)

out_combined_for_plot %>%
  mutate(Model = model, 
         Model = ifelse(Model == "bym2", "Smoothed Direct", Model),
         Model = ifelse(Model == "adaptive bym2", "Proposed", Model)) %>%
  ggplot(aes(x = years, y = est, color = Model)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) + 
  geom_point() + 
  ylab("U5MR") + xlab("Year") + 
  scale_color_manual(values = c("Smoothed Direct" = "grey40", 
                                "Proposed" = "black", "IGME" = "blue", 
                                "Survey Meta" = "red")) +
  annotate("rect", xmin = -Inf, xmax = 1993, ymin = -Inf, ymax = Inf, 
           alpha = 0.5) + 
  annotate("rect", xmin = 1995, xmax = Inf, ymin = -Inf, ymax = Inf, 
           alpha = 0.5) +
  annotate("rect", xmin = 1993, xmax = 1995, ymin = -Inf, ymax = 0.25, 
         alpha = 0.5) + 
  annotate("rect", xmin = 1993, xmax = 1995, ymin = 0.325, ymax = Inf, 
           alpha = 0.5)
ggsave(filename = "../Results/smoothed_direct_comp_2_rwanda.pdf", 
       height = 5, width = 8)

out_combined_for_plot %>%
  mutate(Model = model, 
         Model = ifelse(Model == "bym2", "Smoothed Direct", Model),
         Model = ifelse(Model == "adaptive bym2", "Proposed", Model)) %>%
  ggplot(aes(x = years, y = est, color = Model)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) + 
  geom_point() + 
  ylab("U5MR") + xlab("Year") + 
  scale_color_manual(values = c("Smoothed Direct" = "grey40", 
                                "Proposed" = "black", "IGME" = "blue", 
                                "Survey Meta" = "red")) +
  annotate("rect", xmin = -Inf, xmax = 2014, ymin = -Inf, ymax = Inf, 
           alpha = 0.5) + 
  annotate("rect", xmin = 2020, xmax = Inf, ymin = -Inf, ymax = Inf, 
           alpha = 0.5) +
  annotate("rect", xmin = 2014, xmax = 2020, ymin = -Inf, ymax = 0.015, 
           alpha = 0.5) + 
  annotate("rect", xmin = 2014, xmax = 2020, ymin = 0.09, ymax = Inf, 
           alpha = 0.5)
ggsave(filename = "../Results/smoothed_direct_comp_3_rwanda.pdf", 
       height = 5, width = 8)

    

