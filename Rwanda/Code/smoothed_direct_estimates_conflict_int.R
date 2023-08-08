library(SUMMER)
library(readstata13)
library(dplyr)

#### Load helper functions and data ####
load(file = "../Data/generated_data/direct_estimates.RData")
source("../../helper_files/calc_theta_prior.R")
source("../../helper_files/adaptive_bym2_rgeneric.R")
source("../../helper_files/smoothDirect_extra.R")

num_surveys <- 6
survey_years <- c(1992, 2000, 2005, 2008, 2010, 2015)
years <- levels(births_list[[1]]$time)
num_years <- length(years)

#### Get smoothed direct estimates, bym2 ####
direct_estimates_bym2format <- direct_estimates %>%
  mutate(region = years, years = 0, region_num = as.numeric(region) - 1984)
amat <- matrix(0, nrow = num_years, ncol = num_years)
for(i in 1:num_years){
  if(i < num_years){
    amat[i, i + 1] <- 1
    amat[i + 1, i] <- 1
  }
}
colnames(amat) <- years
rownames(amat) <- years

smoothed_direct_bym2 <- smoothDirect_extra(data = direct_estimates_bym2format, 
                                           Amat = amat, 
                                           time.model = NULL,
                                           year_label = "0", 
                                           year_range = c(0), 
                                           is.yearly = FALSE, m = 1, 
                                           survey.effect = TRUE,
                                           verbose = FALSE)
out_bym2 <- getSmoothed(smoothed_direct_bym2) %>% mutate(years = region)

pc.u <- 1
pc.alpha <- 0.01
pc.u.phi <- 0.5
pc.alpha.phi <- 2/3
hyperpc2 <- list(prec = list(prior = "pc.prec", 
                             param = c(pc.u, pc.alpha)), 
                 phi = list(prior = "pc", 
                            param = c(pc.u.phi, pc.alpha.phi)))
hyperpc1 <- list(prec = list(prior = "pc.prec", 
                             param = c(pc.u, pc.alpha)))
mod <- logit.est ~ -1 + intercept_conflict + intercept_nonconflict  + time +
  f(region.struct, graph = amat, model = "bym2", hyper = hyperpc2,
    scale.model = TRUE, adjust.for.con.comp = TRUE) +
  f(survey.id,  model = "iid", hyper = hyperpc1)
options <- list(dic = TRUE, mlik = TRUE, cpo = TRUE, 
                openmp.strategy = "default", return.marginals.predictor = TRUE)
control.inla <- list(strategy = "adaptive", int.strategy = "auto")

dat <- smoothed_direct_bym2$newdata %>% arrange(region) %>%
  select(region, years, logit.est, logit.prec, region.struct, survey.id) %>% 
  rbind(data.frame(region = as.character(2016:2019),
                   years = NA,
                   logit.est = NA,
                   logit.prec = NA, 
                   region.struct = 2016:2019 - 1984,
                   survey.id = NA)) %>%
  mutate(conflict = as.numeric(region) %in% 1993:1999) %>%
  mutate(intercept_conflict = as.numeric(conflict == TRUE),
         intercept_nonconflict = as.numeric(conflict == FALSE),
         time = region.struct)

fit_conflict_int <- INLA::inla(mod, family = "gaussian", 
                               control.compute = options, 
                               data = dat, 
                               control.predictor = list(compute = TRUE), 
                               control.family = 
                                 list(hyper =  
                                        list(prec = list(initial = log(1),
                                                         fixed = TRUE))), 
                               scale = dat$logit.prec, 
                               control.inla = control.inla, verbose = FALSE)

out_conflict_int <- dat %>% select(region, years, conflict) %>%
  mutate(median = NA, lower = NA, upper = NA, logit.median = NA,
         logit.lower = NA, logit.upper = NA)
for (i in 1:nrow(dat)) {
  tmp.logit <- 
    INLA::inla.rmarginal(1e+05, 
                         fit_conflict_int$marginals.fitted.values[[i]])
  tmp <- expit(tmp.logit)
  out_conflict_int$median[i] <- median(tmp)
  out_conflict_int$lower[i] <- quantile(tmp, probs = 0.05)
  out_conflict_int$upper[i] <- quantile(tmp, probs = 0.95)
  out_conflict_int$logit.median[i] <- median(tmp.logit)
  out_conflict_int$logit.lower[i] <- quantile(tmp.logit, probs = 0.05)
  out_conflict_int$logit.upper[i] <- quantile(tmp.logit, probs = 0.95)
}

out_conflict_int <- out_conflict_int %>%
  filter(is.na(years))

save(out_conflict_int, fit_conflict_int, 
     file = "../Results/smoothed_direct_conflict_int.RData")

#### Get smoothed direct estimates, adaptive bym2 ####
# Create structure matrices
conflict_years <- 1993:1999 - 1984
conflict_years_long <- rep(0, num_years)
conflict_years_long[conflict_years] <- 1

R_conflict <- matrix(0, num_years, num_years)
R_nonconflict <- matrix(0, num_years, num_years)
for(i in 1:num_years){
  if(i == 1){
    if(conflict_years_long[i] || conflict_years_long[i + 1]){
      R_conflict[i, i] <- 1
    }
    else{
      R_nonconflict[i, i] <- 1
    }
  }
  else if(i == num_years){
    if(conflict_years_long[i] || conflict_years_long[i - 1]){
      R_conflict[i, i] <- 1
    }
    else{
      R_nonconflict[i, i] <- 1
    }
  }
  else{
    if(conflict_years_long[i]|| (conflict_years_long[i - 1] && 
                                 conflict_years_long[i + 1])){
      R_conflict[i, i] <- 2
      
    }
    else if(conflict_years_long[i - 1] || conflict_years_long[i + 1]){
      R_conflict[i, i] <- 1
      R_nonconflict[i, i] <- 1
    }
    else{
      R_nonconflict[i, i] <- 2
    }
  }
  
  for(j in 1:num_years){
    if(abs(i - j) == 1){
      if(conflict_years_long[i] || conflict_years_long[j]){
        R_conflict[i, j] <- -1
      }
      else{
        R_nonconflict[i, j] <- -1
      }
    }
  }
}

R_1 <- R_nonconflict
R_2 <- R_conflict
scaled_Q <- INLA:::inla.scale.model.bym.internal(R_1 + R_2, 
                                                 adjust.for.con.comp = TRUE)$Q
gv <- scaled_Q[1, 1] / (R_1 + R_2)[1, 1]
R_1_star <- gv * R_1
R_2_star <- gv * R_2
vals <- (1:num_years)[-num_years]
R_1_star_hat <- R_1_star[vals, vals]
R_2_star_hat <- R_2_star[vals, vals]
eps <- eigen(solve(R_1_star_hat + R_2_star_hat) %*% R_2_star_hat)$values
gamma_tilde <- c(1 / eigen(R_1_star + R_2_star)$values[1:(num_years - 1)], 0)
save(R_1_star, R_2_star, eps, gamma_tilde,
     file = "../Data/generated_data/structure_matrices.RData")

# Specify PC prior hyperparameters
pc.u <- 1
pc.alpha <- 0.01
pc.u.phi <- 0.5
pc.alpha.phi <- 2 / 3
pc.u.theta <- 0.75
pc.alpha.theta <- 0.75

# Fit model 
smoothed_direct_adaptive_bym2 <- smoothed_direct_bym2
adaptive_bym2_model <- 
  INLA::inla.rgeneric.define(model = inla.rgeneric.adaptive.bym2.model, 
                             n = num_years, R_1_star = R_1_star, 
                             R_2_star = R_2_star, gamma_tilde = gamma_tilde, 
                             eps = eps, U_prec = pc.u, alpha_prec = pc.alpha, 
                             U_phi = pc.u.phi, alpha_phi = pc.alpha.phi,
                             U_theta = pc.u.theta, 
                             alpha_theta = pc.alpha.theta)
constr <- list(A = matrix(c(rep(0, num_years), rep(1, num_years)), 
                          nrow = 1, ncol = 2 * num_years), e = 0)
mod <- logit.est ~ -1 + intercept_conflict + intercept_nonconflict  + time +
  f(region.struct, model = adaptive_bym2_model, 
    diagonal = 1e-06, extraconstr = constr, n = 2 * num_years) +
  f(survey.id,  model = "iid", hyper = hyperpc1)
options <- list(dic = TRUE, mlik = TRUE, cpo = TRUE, 
                openmp.strategy = "default", return.marginals.predictor = TRUE)
control.inla <- list(strategy = "adaptive", int.strategy = "auto")
fit_adaptive_conflict_int <- 
  INLA::inla(mod, family = "gaussian", control.compute = options, 
                  data = dat, 
                  control.predictor = list(compute = TRUE), 
                  control.family = 
                    list(hyper =  list(prec = list(initial = log(1), 
                                                   fixed = TRUE))), 
                  scale = dat$logit.prec, 
                  control.inla = control.inla, verbose = FALSE)

out_adaptive_conflict_int <- dat %>% select(region, years, conflict) %>%
  mutate(median = NA, lower = NA, upper = NA, logit.median = NA,
         logit.lower = NA, logit.upper = NA)
for (i in 1:nrow(dat)) {
  tmp.logit <- 
    INLA::inla.rmarginal(1e+05, 
                         fit_adaptive_conflict_int$marginals.fitted.values[[i]])
  tmp <- expit(tmp.logit)
  out_adaptive_conflict_int$median[i] <- median(tmp)
  out_adaptive_conflict_int$lower[i] <- quantile(tmp, probs = 0.05)
  out_adaptive_conflict_int$upper[i] <- quantile(tmp, probs = 0.95)
  out_adaptive_conflict_int$logit.median[i] <- median(tmp.logit)
  out_adaptive_conflict_int$logit.lower[i] <- quantile(tmp.logit, probs = 0.05)
  out_adaptive_conflict_int$logit.upper[i] <- quantile(tmp.logit, probs = 0.95)
}

out_adaptive_conflict_int <- out_adaptive_conflict_int %>%
  filter(is.na(years))

save(out_adaptive_conflict_int, fit_adaptive_conflict_int, 
     file = "../Results/smoothed_direct_adaptive_conflict_int.RData")
