library(SUMMER)
library(readstata13)
library(dplyr)
library(rgdal)

#### Load helper functions and data ####
load(file = "../Data/generated_data/direct_estimates.RData")
load(file = "../Data/generated_data/admin_info.RData")
source("../../helper_files/calc_theta_prior.R")
source("../../helper_files/adaptive_bym2_rgeneric.R")
source("../../helper_files/smoothDirect_extra.R")

#### Get smoothed direct estimates, bym2 ####
direct_estimates_2010_2014 <- direct_estimates %>% 
  bind_rows() %>% filter(years == "10-14") %>%
  filter(!grepl("All", region))
amat <- amat_combined
smoothed_direct_bym2 <- smoothDirect_extra(data = direct_estimates_2010_2014,
                                           Amat = amat, 
                                           time.model = NULL,
                                           year_label = "10-14",
                                           year_range = c(2010, 2014), 
                                           is.yearly = FALSE, m = 5, 
                                           verbose = FALSE)
out_bym2 <- getSmoothed(smoothed_direct_bym2)

#### Get smoothed direct estimates, adaptive bym2 ####
# Create structure matrices
num_region <- nrow(country_admin_table)
W_within <- amat
W_between <- amat
for(i in 1:num_region){
  for(j in i:num_region){
    if(amat[i, j] == 1){
      if(country_admin_table[i, "Country"] == 
         country_admin_table[j, "Country"]){
        W_between[i, j] = 0
        W_between[j, i] = 0
      }
      else{
        W_within[i, j] = 0
        W_within[j, i] = 0
      }
    }
  }
}

D_within <- diag(rowSums(W_within))
D_between <- diag(rowSums(W_between))
R_1 <- D_within - W_within
R_2 <- D_between - W_between
scaled_Q <- INLA:::inla.scale.model.bym.internal(R_1 + R_2, 
                                                 adjust.for.con.comp = TRUE)$Q
gv <- scaled_Q[1, 1] / (R_1 + R_2)[1, 1]
R_1_star <- gv * (D_within - W_within)
R_2_star <- gv * (D_between - W_between)
vals <- (1:num_region)[-num_region]
R_1_star_hat <- R_1_star[vals, vals]
R_2_star_hat <- R_2_star[vals, vals]
eps <- eigen(solve(R_1_star_hat + R_2_star_hat) %*% R_2_star_hat)$values
gamma_tilde <- c(1 / eigen(R_1_star + R_2_star)$values[1:(num_region - 1)], 0)
save(R_1_star, R_2_star, eps, gamma_tilde,
     file = "../Data/generated_data/structure_matrices.RData")

# Specify PC prior hyperparameters
pc.u <- 1
pc.alpha <- 0.01
pc.u.phi <- 0.5
pc.alpha.phi <- 2 / 3
pc.u.theta <- 0.75
pc.alpha.theta <- 0.75
theta_prior <- calc_theta_prior(eps, pc.u.theta, pc.alpha.theta)
plot(theta_prior$theta, exp(theta_prior$log_density), type = "l")

# Fit model 
smoothed_direct_adaptive_bym2 <- smoothed_direct_bym2
adaptive_bym2_model <- 
  INLA::inla.rgeneric.define(model = inla.rgeneric.adaptive.bym2.model, 
                             n = num_region, R_1_star = R_1_star, 
                             R_2_star = R_2_star, gamma_tilde = gamma_tilde, 
                             eps = eps, U_prec = pc.u, alpha_prec = pc.alpha, 
                             U_phi = pc.u.phi, alpha_phi = pc.alpha.phi,
                             U_theta = pc.u.theta, 
                             alpha_theta = pc.alpha.theta)
constr <- list(A = matrix(c(rep(0, num_region), rep(1, num_region)), 
                          nrow = 1, ncol = 2 * num_region), e = 0)
mod <- logit.est ~ f(region.struct, model = adaptive_bym2_model, 
                     diagonal = 1e-06, extraconstr = constr, n = 2 * num_region)
options <- list(dic = TRUE, mlik = TRUE, cpo = TRUE, 
                openmp.strategy = "default")
control.inla <- list(strategy = "adaptive", int.strategy = "auto")
fit <- INLA::inla(mod, family = "gaussian", control.compute = options, 
                  data = smoothed_direct_adaptive_bym2$newdata, 
                  control.predictor = list(compute = TRUE), 
                  control.family = 
                    list(hyper =  list(prec = list(initial = log(1), 
                                                   fixed = TRUE))), 
                  scale = smoothed_direct_adaptive_bym2$newdata$logit.prec, 
                  lincomb = smoothed_direct_adaptive_bym2$lincombs.fit, 
                  control.inla = control.inla, verbose = FALSE)
smoothed_direct_adaptive_bym2$model <- mod
smoothed_direct_adaptive_bym2$fit <- fit
out_adaptive_bym2 <- getSmoothed(smoothed_direct_adaptive_bym2)

out_combined <- rbind(cbind(out_bym2, prior = "bym2"),
                      cbind(out_adaptive_bym2, prior = "adaptive bym2"))
save(out_combined, smoothed_direct_bym2, smoothed_direct_adaptive_bym2, 
     file = "../Results/smoothed_direct_estimates.RData")

#### Get smoothed direct estimates, bym2 with country intercept ####
hyperpc2 <- list(prec = list(prior = "pc.prec", 
                             param = c(pc.u, pc.alpha)), 
                 phi = list(prior = "pc", 
                            param = c(pc.u.phi, pc.alpha.phi)))
mod <- logit.est ~ -1 + intercept_Burundi + intercept_Ethiopia + 
  intercept_Kenya + intercept_Rwanda + intercept_Tanzania + intercept_Uganda +
  f(region.struct, graph = amat, model = "bym2", hyper = hyperpc2, 
    scale.model = TRUE, adjust.for.con.comp = TRUE)
options <- list(dic = TRUE, mlik = TRUE, cpo = TRUE, 
                openmp.strategy = "default", return.marginals.predictor = TRUE)
dat_country_int <- smoothed_direct_bym2$newdata %>% arrange(region) %>%
  filter(!is.na(years)) %>% cbind(Country = country_admin_table$Country) %>%
  mutate(intercept_Burundi = as.numeric(Country == "Burundi"), 
         intercept_Ethiopia = as.numeric(Country == "Ethiopia"),
         intercept_Kenya = as.numeric(Country == "Kenya"),
         intercept_Rwanda = as.numeric(Country == "Rwanda"),
         intercept_Tanzania = as.numeric(Country == "Tanzania"),
         intercept_Uganda = as.numeric(Country == "Uganda"))
fit_country_int <- INLA::inla(mod, family = "gaussian", 
                              control.compute = options, 
                              data = dat_country_int, 
                              control.predictor = list(compute = TRUE), 
                              control.family = 
                                list(hyper =  list(prec = list(initial = log(1), 
                                                               fixed = TRUE))), 
                              scale = dat_country_int$logit.prec, 
                              control.inla = control.inla, verbose = FALSE)

out_country_int <- dat_country_int %>% select(region, years, Country) %>%
  mutate(median = NA, lower = NA, upper = NA, logit.median = NA,
         logit.lower = NA, logit.upper = NA)
for (i in 1:nrow(dat_country_int)) {
  tmp.logit <- 
    INLA::inla.rmarginal(1e+05, 
                         fit_country_int$marginals.fitted.values[[i]])
  tmp <- expit(tmp.logit)
  out_country_int$median[i] <- median(tmp)
  out_country_int$lower[i] <- quantile(tmp, probs = 0.025)
  out_country_int$upper[i] <- quantile(tmp, probs = 0.975)
  out_country_int$logit.median[i] <- median(tmp.logit)
  out_country_int$logit.lower[i] <- quantile(tmp.logit, probs = 0.025)
  out_country_int$logit.upper[i] <- quantile(tmp.logit, probs = 0.975)
}

save(out_country_int, fit_country_int, 
     file = "../Results/smoothed_direct_country_int.RData")


#### Get smoothed direct estimates, adaptive bym2 with country intercept ####
mod_country_int <- logit.est ~ -1 + intercept_Burundi + intercept_Ethiopia + 
  intercept_Kenya + intercept_Rwanda + intercept_Tanzania + intercept_Uganda +
  f(region.struct, model = adaptive_bym2_model, 
    diagonal = 1e-06, extraconstr = constr, n = 2 * num_region)
fit_adaptive_country_int <- 
  INLA::inla(mod_country_int, family = "gaussian", control.compute = options, 
             data = dat_country_int, control.predictor = list(compute = TRUE), 
             control.family = 
               list(hyper =  list(prec = list(initial = log(1), fixed = TRUE))), 
             scale = dat_country_int$logit.prec, 
             control.inla = control.inla, verbose = FALSE)
out_adaptive_country_int <- dat_country_int %>% 
  select(region, years, Country) %>%
  mutate(median = NA, lower = NA, upper = NA, logit.median = NA,
         logit.lower = NA, logit.upper = NA)
for (i in 1:nrow(dat_country_int)) {
  tmp.logit <- 
    INLA::inla.rmarginal(1e+05, 
                         fit_adaptive_country_int$marginals.fitted.values[[i]])
  tmp <- expit(tmp.logit)
  out_adaptive_country_int$median[i] <- median(tmp)
  out_adaptive_country_int$lower[i] <- quantile(tmp, probs = 0.025)
  out_adaptive_country_int$upper[i] <- quantile(tmp, probs = 0.975)
  out_adaptive_country_int$logit.median[i] <- median(tmp.logit)
  out_adaptive_country_int$logit.lower[i] <- quantile(tmp.logit, probs = 0.025)
  out_adaptive_country_int$logit.upper[i] <- quantile(tmp.logit, probs = 0.975)
}

save(out_adaptive_country_int, fit_adaptive_country_int, 
     file = "../Results/smoothed_direct_adaptive_country_int.RData")
