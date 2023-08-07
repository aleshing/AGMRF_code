library(dplyr)
library(tidyverse)

#### Load in results ####
load(file = paste0("../Results/model_comp"))

#### Plot results ####
results <- results %>% 
  mutate(prec_setting_descrip = factor(prec_setting, 
                                       levels = c("3", "1", "2")))
levels(results$prec_setting_descrip) <- c("V: 1 / 75", 
                                          "V: 1 / 150", 
                                          "V: 1 / 300")

results %>% group_by(prec_setting, par_setting, model) %>% 
  summarise(rmse_logit_risk = mean(rmse_logit_risk),
            rmse_logit_obs = mean(rmse_logit_obs),
            rmse_risk = mean(rmse_risk),
            rmse_obs = mean(rmse_obs),
            cpo = mean(cpo),
            dic = mean(dic)) %>% as.data.frame()

results %>% filter(re_settings == "Random Effect: Equal Precisions") %>%
  ggplot(aes(x = model, y = rmse_logit_risk)) + geom_boxplot() + 
  facet_wrap(~ trend_settings + prec_setting_descrip, ncol = 3) +
  xlab("Model") + ylab("RMSE") + theme_gray(base_size = 16)
ggsave(filename = "../Plots/rmse_logit_risk_plots_equal.pdf", 
       width = 10, height = 7)
results %>% filter(re_settings == "Random Effect: Equal Precisions") %>%
  ggplot(aes(x = model, y = cpo)) + geom_boxplot() + 
  facet_wrap(~ trend_settings + prec_setting_descrip, ncol = 3) +
  xlab("Model") + ylab("LS") + theme_gray(base_size = 16)
ggsave(filename = "../Plots/cpo_plots_equal.pdf", width = 10, height = 7)
results %>%  filter(re_settings == "Random Effect: Equal Precisions") %>%
  ggplot(aes(x = model, y = dic)) + geom_boxplot() + 
  facet_wrap(~ trend_settings + prec_setting_descrip, ncol = 3) +
  xlab("Model") + ylab("DIC") + theme_gray(base_size = 16)
ggsave(filename = "../Plots/dic_plots_equal.pdf", width = 10, height = 7)

results %>% filter(re_settings == "Random Effect: Unequal Precisions") %>%
  ggplot(aes(x = model, y = rmse_logit_risk)) + geom_boxplot() + 
  facet_wrap(~ trend_settings + prec_setting_descrip, ncol = 3) +
  xlab("Model") + ylab("RMSE") + theme_gray(base_size = 16)
ggsave(filename = "../Plots/rmse_logit_risk_plots_unequal.pdf", 
       width = 10, height = 7)
results %>% filter(re_settings == "Random Effect: Unequal Precisions") %>%
  ggplot(aes(x = model, y = cpo)) + geom_boxplot() + 
  facet_wrap(~ trend_settings + prec_setting_descrip, ncol = 3) +
  xlab("Model") + ylab("LS") + theme_gray(base_size = 16)
ggsave(filename = "../Plots/cpo_plots_unequal.pdf", width = 10, height = 7)
results %>%  filter(re_settings == "Random Effect: Unequal Precisions") %>%
  ggplot(aes(x = model, y = dic)) + geom_boxplot() + 
  facet_wrap(~ trend_settings + prec_setting_descrip, ncol = 3) +
  xlab("Model") + ylab("DIC") + theme_gray(base_size = 16)
ggsave(filename = "../Plots/dic_plots_unequal.pdf", width = 10, height = 7)

#### Plot results as differences ####
results_wider <- results %>% select(-c(rmse_risk, rmse_obs, rmse_logit_obs)) %>%
  pivot_wider(names_from = model, values_from = c(rmse_logit_risk, cpo, dic))

results_wider %>% filter(re_settings == "Random Effect: Equal Precisions") %>%
  ggplot(aes(y = `rmse_logit_risk_Smoothed Direct` - 
               `rmse_logit_risk_Proposed`)) + geom_boxplot() + 
  facet_wrap(~ trend_settings + prec_setting_descrip, ncol = 3) +
  ylab("RMSE Difference") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  theme_gray(base_size = 16)
ggsave(filename = "../Plots/rmse_logit_risk_plots_equal_diff.pdf", 
       width = 10, height = 7)
results_wider %>% filter(re_settings == "Random Effect: Equal Precisions") %>%
  ggplot(aes(y = `cpo_Smoothed Direct` - `cpo_Proposed`)) + geom_boxplot() + 
  facet_wrap(~ trend_settings + prec_setting_descrip, ncol = 3) +
  ylab("LS Difference") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  theme_gray(base_size = 16)
ggsave(filename = "../Plots/cpo_plots_equal_diff.pdf", width = 10, height = 7)
results_wider %>%  filter(re_settings == "Random Effect: Equal Precisions") %>%
  ggplot(aes(y = `dic_Smoothed Direct` - `dic_Proposed`)) + geom_boxplot() + 
  facet_wrap(~ trend_settings + prec_setting_descrip, ncol = 3) +
  ylab("DIC Difference") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  theme_gray(base_size = 16)
ggsave(filename = "../Plots/dic_plots_equal_diff.pdf", width = 10, height = 7)

results_wider %>% filter(re_settings == "Random Effect: Unequal Precisions") %>%
  ggplot(aes(y = `rmse_logit_risk_Smoothed Direct` -
               `rmse_logit_risk_Proposed`)) + geom_boxplot() + 
  facet_wrap(~ trend_settings + prec_setting_descrip, ncol = 3) +
  ylab("RMSE Difference") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  theme_gray(base_size = 16)
ggsave(filename = "../Plots/rmse_logit_risk_plots_unequal_diff.pdf", 
       width = 10, height = 7)
results_wider %>% filter(re_settings == "Random Effect: Unequal Precisions") %>%
  ggplot(aes(y = `cpo_Smoothed Direct` - `cpo_Proposed`)) + geom_boxplot() + 
  facet_wrap(~ trend_settings + prec_setting_descrip, ncol = 3) +
  ylab("LS Difference") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  theme_gray(base_size = 16)
ggsave(filename = "../Plots/cpo_plots_unequal_diff.pdf", width = 10, height = 7)
results_wider %>%  
  filter(re_settings == "Random Effect: Unequal Precisions") %>%
  ggplot(aes(y = `dic_Smoothed Direct` - `dic_Proposed`)) + geom_boxplot() + 
  facet_wrap(~ trend_settings + prec_setting_descrip, ncol = 3) +
  ylab("DIC Difference") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  theme_gray(base_size = 16)
ggsave(filename = "../Plots/dic_plots_unequal_diff.pdf", width = 10, height = 7)

#### Trend plots #### 
trend_dat <- data.frame(time = rep(1:30, 3), mu = NA, 
                        trend_name = 
                          rep(c("Constant", "Level Change", "Triangle"),
                              each = 30))
trend_dat$mu[1:30] <- -2
trend_dat$mu[31:60] <- -2
trend_dat$mu[(30 + 9):(30 + 15)] <- -1
trend_dat$mu[61:90] <- c(rep(-2, 8), 
                         seq(from = -2, to = -0.75, length.out = 5)[2:4], 
                         -0.75, seq(from = -2, to = -0.75, length.out = 5)[4:2], 
                         rep(-2, 15))

trend_dat %>% ggplot(aes(x = time, y = mu)) + geom_line() +
  facet_wrap(vars(trend_name), ncol = 3) + ylab(expression(mu)) + xlab("Time") 
ggsave(filename = "../Plots/trends.pdf", width = 7, height = 3)
