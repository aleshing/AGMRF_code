library(SUMMER)
library(tidyverse)
library(rgdal)
library(sf)
library(ggpubr)
library(xtable)
library(INLA)

#### Load results ####
load(file = "../Data/generated_data/direct_estimates.RData")
load(file = "../Data/generated_data/structure_matrices.RData")
load(file = "../Data/generated_data/structure_matrices_general.RData")
load(file = "../Results/smoothed_direct_estimates.RData")
load(file = "../Results/smoothed_direct_separate.RData")
load(file = "../Results/smoothed_direct_country_int.RData")
load(file = "../Results/smoothed_direct_estimates_general.RData")
load(file = "../Results/smoothed_direct_adaptive_country_int.RData")
load(file = "../Data/generated_data/admin_info.RData")
source("../../helper_files/calc_theta_prior.R")

#### Master data frame ####
out_master <- out_separate %>%
  mutate(Estimate = median, 
         Model = "Smoothed Direct Country-Specific") %>%
  select(region, Estimate, lower, upper, Model) %>%
  rbind(out_country_int %>%
          mutate(Estimate = median, 
                 Model = "Smoothed Direct Country-Intercept") %>%
          select(region, Estimate, lower, upper, Model)) %>%
  rbind(out_combined %>%
          mutate(Estimate = median, 
                 Model = ifelse(prior == "bym2",
                                "Smoothed Direct",
                                "Proposed")) %>%
          select(region, Estimate, lower, upper, Model)) %>%
  rbind(out_adaptive_country_int %>%
          mutate(Estimate = median, 
                 Model = "Proposed Country-Intercept") %>%
          select(region, Estimate, lower, upper, Model)) %>%
  rbind(out_general %>%
          mutate(Estimate = median, 
                 Model = "Proposed General") %>%
          select(region, Estimate, lower, upper, Model))

#### Smoothed Direct Separate Parameter Estimates ####
num_countries <- length(unique(country_admin_table$Country))
phi_marginals <- vector(mode = "list", length = num_countries)

summaries_sd_separate <- 
  data.frame(Model = "Smoothed Direct Country-Specific",
             Country = rep(unique(country_admin_table$Country), each = 4),
             Parameter = rep(c("mu", "tau", "phi", "theta"), num_countries),
             Mean = NA, SD = NA, `2.5% Quantile` = NA, Median = NA, 
             `97.5% Quantile` = NA, Mode = NA)
for(i in 1:num_countries){
  inds_fixed <- (i - 1) * 4 + 1
  inds_hyperpar <- ((i - 1) * 4 + 2):((i - 1) * 4 + 3)
  summaries_sd_separate[inds_fixed, 4:9] <- 
    smoothed_direct_separate[[i]]$fit$summary.fixed[, 1:6] %>% unname()
  summaries_sd_separate[inds_hyperpar, 4:9] <- 
    smoothed_direct_separate[[i]]$fit$summary.hyperpar[1:2, 1:6] %>% unname() 
  phi_marginals[[i]] <- 
    smoothed_direct_separate[[i]]$fit$marginals.hyperpar[[2]]
}

plot(phi_marginals[[1]][, 1], phi_marginals[[1]][, 2], type = "l")
for(i in 2:num_countries){
  lines(phi_marginals[[i]][, 1], phi_marginals[[1]][, 2])
}

#### Smoothed Direct Country Intercept Parameter Estimates ####
summaries_sd_country_int <- 
  data.frame(Model = "Smoothed Direct Country-Intercept",
             Country = c(unique(country_admin_table$Country), rep(NA, 3)),
             Parameter = c(rep("mu", num_countries), "tau", "phi", "theta"),
             Mean = NA, SD = NA, `2.5% Quantile` = NA, Median = NA, 
             `97.5% Quantile` = NA, Mode = NA)
summaries_sd_country_int[1:6, 4:9] <- 
  fit_country_int$summary.fixed[, 1:6] %>% unname()
summaries_sd_country_int[7:8, 4:9] <- 
  fit_country_int$summary.hyperpar[1:2, 1:6] %>% unname()

#### Smoothed Direct Parameter Estimates ####
summaries_sd <- 
  data.frame(Model = "Smoothed Direct",
             Country = NA,
             Parameter = c("mu", "tau", "phi", "theta"),
             Mean = NA, SD = NA, `2.5% Quantile` = NA, Median = NA, 
             `97.5% Quantile` = NA, Mode = NA)
summaries_sd[1, 4:9] <- 
  smoothed_direct_bym2$fit$summary.fixed[, 1:6] %>% unname()
summaries_sd[2:3, 4:9] <- 
  smoothed_direct_bym2$fit$summary.hyperpar[1:2, 1:6] %>% unname()

#### Proposed Parameter Estimates ####
summaries_proposed <- 
  data.frame(Model = "Proposed",
             Country = NA,
             Parameter = c("mu", "tau", "phi", "theta"),
             Mean = NA, SD = NA, `2.5% Quantile` = NA, Median = NA, 
             `97.5% Quantile` = NA, Mode = NA)
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

theta_post <- 
  inla.tmarginal(fun = SUMMER::expit, 
                 marginal = smoothed_direct_adaptive_bym2$fit$marginals.hyperpar$`Theta3 for region.struct`)
theta_sum <- inla.zmarginal(theta_post)
theta_mode <- inla.mmarginal(theta_post)

summaries_proposed[1, 4:9] <- 
  smoothed_direct_adaptive_bym2$fit$summary.fixed[, 1:6] %>% unname()
summaries_proposed[2, 4:9] <- 
  c(tau_sum[c(1:3, 5, 7)], tau_mode)
summaries_proposed[3, 4:9] <- 
  c(phi_sum[c(1:3, 5, 7)], phi_mode)
summaries_proposed[4, 4:9] <- 
  c(theta_sum[c(1:3, 5, 7)], theta_mode)

# Compare theta prior and post
pc.u.theta <- 0.75
pc.alpha.theta <- 0.75
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
inla.mmarginal(data.frame(x = theta_prior$theta, 
                          y = exp(theta_prior$log_density)))
inla.qmarginal(p = c(0.025, 0.5, 0.975), 
               data.frame(x = theta_prior$theta, 
                          y = exp(theta_prior$log_density)))
theta_dens <- rbind(data.frame(theta = theta_prior$theta, 
                               density = exp(theta_prior$log_density), 
                               Which = "Prior"),
                    data.frame(theta = theta_post[, 1], 
                               density = theta_post[, 2], 
                               Which = "Posterior"))
theta_dens %>% ggplot(aes(x = theta, y = density, colour = Which)) + 
  geom_line() + xlab("Theta") + ylab("Density") 
ggsave(filename = "../Results/theta_prior_post_no_country-int.pdf",
       height = 5, width = 8)

# Compare smoothed direct and proposed
plot(smoothed_direct_bym2$fit$marginals.hyperpar[[2]][, 1], 
     smoothed_direct_bym2$fit$marginals.hyperpar[[2]][, 2], type = "l")
lines(phi_post[, 1], phi_post[, 2], col = "red")

#### Proposed Country Intercept Parameter Estimates ####
summaries_proposed_country_intercept <- 
  data.frame(Model = "Proposed Country-Intercept",
             Country = c(unique(country_admin_table$Country), rep(NA, 3)),
             Parameter = c(rep("mu", num_countries), "tau", "phi", "theta"),
             Mean = NA, SD = NA, `2.5% Quantile` = NA, Median = NA, 
             `97.5% Quantile` = NA, Mode = NA)
tau_post <- 
  inla.tmarginal(fun = exp, 
                 marginal = fit_adaptive_country_int$marginals.hyperpar$`Theta1 for region.struct`)
tau_sum <- inla.zmarginal(tau_post)
tau_mode <- inla.mmarginal(tau_post)

phi_post <- 
  inla.tmarginal(fun = SUMMER::expit, 
                 marginal = fit_adaptive_country_int$marginals.hyperpar$`Theta2 for region.struct`)
phi_sum <- inla.zmarginal(phi_post)
phi_mode <- inla.mmarginal(phi_post)

theta_post <- 
  inla.tmarginal(fun = SUMMER::expit, 
                 marginal = fit_adaptive_country_int$marginals.hyperpar$`Theta3 for region.struct`)
theta_sum <- inla.zmarginal(theta_post)
theta_mode <- inla.mmarginal(theta_post)

summaries_proposed_country_intercept[1:6, 4:9] <- 
  fit_adaptive_country_int$summary.fixed[, 1:6] %>% unname()
summaries_proposed_country_intercept[7, 4:9] <- 
  c(tau_sum[c(1:3, 5, 7)], tau_mode)
summaries_proposed_country_intercept[8, 4:9] <- 
  c(phi_sum[c(1:3, 5, 7)], phi_mode)
summaries_proposed_country_intercept[9, 4:9] <- 
  c(theta_sum[c(1:3, 5, 7)], theta_mode)

# Compare theta prior and post
pc.u.theta <- 0.75
pc.alpha.theta <- 0.75
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
inla.mmarginal(data.frame(x = theta_prior$theta, 
                          y = exp(theta_prior$log_density)))
inla.qmarginal(p = c(0.025, 0.5, 0.975), 
               data.frame(x = theta_prior$theta, 
                          y = exp(theta_prior$log_density)))
theta_dens <- rbind(data.frame(theta = theta_prior$theta, 
                               density = exp(theta_prior$log_density), 
                               Which = "Prior"),
                    data.frame(theta = theta_post[, 1], 
                               density = theta_post[, 2], 
                               Which = "Posterior"))
theta_dens %>% ggplot(aes(x = theta, y = density, colour = Which)) + 
  geom_line() + xlab("Theta") + ylab("Density") 
ggsave(filename = "../Results/theta_prior_post_country-int.pdf",
       height = 5, width = 8)

# Compare smoothed direct and proposed
plot(fit_country_int$marginals.hyperpar[[2]][, 1], 
     fit_country_int$marginals.hyperpar[[2]][, 2], type = "l")
lines(phi_post[, 1], phi_post[, 2], col = "red")

#### Proposed Parameter Estimates General ####
summaries_general <- 
  data.frame(Model = "Proposed General",
             Country = c(rep(NA, 3), unique(country_admin_table$Country)),
             Parameter = c("mu", "phi", "theta", "tau", 
                           rep("psi", num_countries - 1)),
             Mean = NA, SD = NA, `2.5% Quantile` = NA, Median = NA, 
             `97.5% Quantile` = NA, Mode = NA)
tau_post <- 
  inla.tmarginal(fun = exp, 
                 marginal = smoothed_direct_adaptive_bym2_general$fit$marginals.hyperpar$`Theta1 for region.struct`)
tau_sum <- inla.zmarginal(tau_post)
tau_mode <- inla.mmarginal(tau_post)

phi_post <- 
  inla.tmarginal(fun = SUMMER::expit, 
                 marginal = smoothed_direct_adaptive_bym2_general$fit$marginals.hyperpar$`Theta2 for region.struct`)
phi_sum <- inla.zmarginal(phi_post)
phi_mode <- inla.mmarginal(phi_post)

theta_post <- 
  inla.tmarginal(fun = SUMMER::expit, 
                 marginal = smoothed_direct_adaptive_bym2_general$fit$marginals.hyperpar$`Theta3 for region.struct`)
theta_sum <- inla.zmarginal(theta_post)
theta_mode <- inla.mmarginal(theta_post)

summaries_general[1, 4:9] <- 
  smoothed_direct_adaptive_bym2_general$fit$summary.fixed[, 1:6] %>% unname()
summaries_general[2, 4:9] <- 
  c(phi_sum[c(1:3, 5, 7)], phi_mode)
summaries_general[3, 4:9] <- 
  c(theta_sum[c(1:3, 5, 7)], theta_mode)
summaries_general[4, 4:9] <- 
  c(tau_sum[c(1:3, 5, 7)], tau_mode)

for(m in 2:num_countries){
  par_name <- paste0("Theta", m + 2," for region.struct")
  psi_post <- 
    inla.tmarginal(fun = exp, 
                   marginal = smoothed_direct_adaptive_bym2_general$fit$marginals.hyperpar[[par_name]])
  psi_sum <- inla.zmarginal(psi_post)
  psi_mode <- inla.mmarginal(psi_post)
  summaries_general[m + 3, 4:9] <- 
    c(psi_sum[c(1:3, 5, 7)], psi_mode)
}

# Compare theta prior and post
pc.u.theta <- 0.75
pc.alpha.theta <- 0.75
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
inla.mmarginal(data.frame(x = theta_prior$theta, 
                          y = exp(theta_prior$log_density)))
inla.qmarginal(p = c(0.025, 0.5, 0.975), 
               data.frame(x = theta_prior$theta, 
                          y = exp(theta_prior$log_density)))
theta_dens <- rbind(data.frame(theta = theta_prior$theta, 
                               density = exp(theta_prior$log_density), 
                               Which = "Prior"),
                    data.frame(theta = theta_post[, 1], 
                               density = theta_post[, 2], 
                               Which = "Posterior"))
theta_dens %>% ggplot(aes(x = theta, y = density, colour = Which)) + 
  geom_line() + xlab("Theta") + ylab("Density") 
ggsave(filename = "../Results/theta_prior_post_general.pdf",
       height = 5, width = 8)

# Compare smoothed direct and proposed general
plot(smoothed_direct_bym2$fit$marginals.hyperpar[[2]][, 1], 
     smoothed_direct_bym2$fit$marginals.hyperpar[[2]][, 2], type = "l")
lines(phi_post[, 1], phi_post[, 2], col = "red")


#### Bind together summaries ####
summaries_master <- rbind(summaries_sd_separate, summaries_sd, 
                          summaries_sd_country_int, summaries_proposed,
                          summaries_proposed_country_intercept, 
                          summaries_general)

print(xtable(summaries_general), include.rownames = FALSE)


#### Do some data munging ####
direct_estimates_2010_2014_country <- direct_estimates %>% 
  bind_rows(.id = "Country") %>% filter(years == "10-14") %>%
  filter(grepl("All", region))

direct_estimates_2010_2014 <- direct_estimates %>% 
  bind_rows(.id = "Country") %>% filter(years == "10-14") %>%
  filter(!grepl("All", region)) %>% group_by(Country) %>%
  mutate(country_mean = mean(mean), 
         country_meta = SUMMER::expit(sum(logit.est * logit.prec) /
                                        sum(logit.prec))) %>% ungroup %>%
  merge(direct_estimates_2010_2014_country %>% 
          mutate(country_direct = mean) %>%
          select(Country, country_direct)) %>%
  select(Country, region, mean, lower, upper, country_direct)

num_regions <- nrow(direct_estimates_2010_2014)
country_admin_table$borders_other_countries <- NA
for(i in 1:num_regions){
  temp_country <- country_admin_table$Country[i]
  temp_countries <- 
    unique(country_admin_table$Country[which(amat_combined[i, ] == 1)])
  if(sum(temp_country != temp_countries) > 0){
    country_admin_table$borders_other_countries[i] <- "Yes"
  }
  else{
    country_admin_table$borders_other_countries[i] <- "No"
  }
}

out_master <- out_master %>%
  merge(country_admin_table %>% 
          mutate(region = Admin1) %>%
          select(region, borders_other_countries), 
        by = "region")

#### Get country outlines ####
countries <- unique(geo_combined$CNTRYNAMEE)
geos <- vector(mode = "list", length = length(countries))
sf::sf_use_s2(FALSE)
for(i in 1:length(countries)){
  geos[[i]] <- geo_combined %>%
    filter(CNTRYNAMEE == countries[i]) %>%
    st_union() %>% st_as_sf() %>%
    mutate(country = countries[i])
}
names(geos) <- countries
geo_country <- bind_rows(geos)

#### Map results ####

# Smoothed direct and proposed general model
geo <- 
  rbind(cbind(geo_combined, 
              value = (out_master %>% 
                         filter(Model == "Smoothed Direct"))$Estimate, 
              Model = "Smoothed Direct"),
        cbind(geo_combined, 
              value = (out_master %>% 
                         filter(Model == "Proposed General"))$Estimate, 
              Model = "Proposed General")) %>%
  mutate(U5MR = 
           cut(value, c(-Inf, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 
                        0.075, 0.08, 0.085, 0.09, 0.095, 0.1, Inf),
               labels = c("<0.04", "0.04 to <0.045", "0.045 to <0.05",
                          "0.05 to <0.055", "0.055 to <0.06",
                          "0.06 to <0.065", "0.065 to <0.07",
                          "0.07 to <0.075", "0.075 to <0.08",
                          "0.08 to <0.85", "0.085 to <0.09", 
                          "0.09 to <0.95", "0.095 to <0.1", ">= 0.1")))
plot_ests <- ggplot(geo) + geom_sf(aes(fill = U5MR)) +  
  scale_fill_viridis_d(option = 'E', direction = -1) + facet_wrap(~Model)  +
  geom_sf(data = geo_country, col = "black", fill = NA, lwd = 1)
ggsave(filename = 
         "../Results/smoothed_direct_comp_no_country_int_1_general.pdf")

plot_ests_rb <- geo %>% filter(CNTRYNAMEE %in% c("Rwanda", "Burundi")) %>%
  ggplot() + geom_sf(aes(fill = U5MR)) +  
  scale_fill_viridis_d(option = 'E', direction = -1) + facet_wrap(~Model)  +
  geom_sf(data = geo_country %>% filter(country %in% c("Rwanda", "Burundi")),
          col = "black", fill = NA, lwd = 1)
ggsave(filename = 
         "../Results/smoothed_direct_comp_no_country_int_rb_1_general.pdf")

# Proposed general and proposed model
geo <- rbind(cbind(geo_combined, 
                   value = (out_master %>% 
                              filter(Model == "Proposed General"))$Estimate, 
                   Model = "Proposed General"),
             cbind(geo_combined, 
                   value = (out_master %>% 
                              filter(Model == "Proposed"))$Estimate, 
                   Model = "Proposed")) %>%
  mutate(U5MR = 
           cut(value, c(-Inf, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 
                        0.075, 0.08, 0.085, 0.09, 0.095, 0.1, Inf),
               labels = c("<0.04", "0.04 to <0.045", "0.045 to <0.05",
                          "0.05 to <0.055", "0.055 to <0.06",
                          "0.06 to <0.065", "0.065 to <0.07",
                          "0.07 to <0.075", "0.075 to <0.08",
                          "0.08 to <0.85", "0.085 to <0.09", 
                          "0.09 to <0.95", "0.095 to <0.1", ">= 0.1")))
plot_ests <- ggplot(geo) + geom_sf(aes(fill = U5MR)) +  
  scale_fill_viridis_d(option = 'E', direction = -1) + facet_wrap(~Model)  +
  geom_sf(data = geo_country, col = "black", fill = NA, lwd = 1)
ggsave(filename = "../Results/smoothed_direct_comp_proposed_general_1.pdf")

plot_ests_rb <- geo %>% filter(CNTRYNAMEE %in% c("Rwanda", "Burundi")) %>%
  ggplot() + geom_sf(aes(fill = U5MR)) +  
  scale_fill_viridis_d(option = 'E', direction = -1) + facet_wrap(~Model)  +
  geom_sf(data = geo_country %>% filter(country %in% c("Rwanda", "Burundi")),
          col = "black", fill = NA, lwd = 1)
ggsave(filename = "../Results/smoothed_direct_comp_proposed_general_rb_1.pdf")

#### Non-map results ####
dodge <- position_dodge(1)

# Smoothed direct and proposed country-intercept model
out_master %>%
  mutate(Country = word(region, sep = "-")) %>%
  merge(direct_estimates_2010_2014 %>% 
          mutate(direct = mean) %>%
          select(c(region, country_direct, direct))) %>%
  mutate(region = word(region, sep = "-", start = 2)) %>%
  filter(Model %in% c("Smoothed Direct Country-Specific", 
                      "Smoothed Direct Country-Intercept",
                      "Proposed Country-Intercept",
                      "Proposed General",
                      "Proposed")) %>%
  ggplot(aes(x = Estimate, y = region, col = Model)) + 
  geom_vline(aes(xintercept = country_direct), col = "red") +
  geom_point(aes(shape = borders_other_countries), 
             position = dodge, cex = 2) +
  geom_linerange(aes(xmin = lower, xmax = upper, col = Model), 
                 position = dodge) +
  xlab("U5MR") + ylab("Region") + 
  geom_point(aes(x = direct), col = "black") + 
  labs(shape = "Borders Another Country") + 
  facet_wrap(facets = vars(Country), scales = "free") 
ggsave(filename = "../Results/smoothed_direct_comp_country_int_2_general.pdf", 
       height = 12, width = 12)
