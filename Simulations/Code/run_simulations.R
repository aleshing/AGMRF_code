library(SUMMER)
library(dplyr)

#### Load data and helper functions ####
load("../Data/structure_matrices.RData")
source("../../helper_files/adaptive_bym2_rgeneric.R")
source("../../helper_files/bym2_rgeneric.R")
source("../../helper_files/smoothDirect_extra.R")

#### Set simulation parameters ####
num_reps <- 200
num_years <- 30
num_par_settings <- 6
num_prec_settings <- 3

#### Specify PC prior hyperparameters ####
pc.u <- 1
pc.alpha <- 0.01
pc.u.phi <- 0.5
pc.alpha.phi <- 2 / 3
pc.u.theta <- 0.75
pc.alpha.theta <- 0.75

#### Get SUMMER internals ####
load("../Data/dat_prec_1_par_1_rep_1")
dat_bym2format <- dat %>%
    mutate(region = years, years = 0, region_num = as.numeric(region))
amat <- matrix(0, nrow = num_years, ncol = num_years)
for(i in 1:num_years){
    if(i < num_years){
        amat[i, i + 1] <- 1
        amat[i + 1, i] <- 1
    }
}
colnames(amat) <- 1:num_years
rownames(amat) <- 1:num_years

smoothed_direct_bym2 <- smoothDirect_extra(data = dat_bym2format, 
                                           Amat = amat, time.model = NULL,
                                           year_label = "0", year_range = c(0), 
                                           is.yearly = FALSE, m = 1, 
                                           verbose = FALSE)

#### Set up INLA details #### 
# AGMRF
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
mod <- logit.est ~ f(region.struct, model = adaptive_bym2_model, 
                     diagonal = 1e-06, extraconstr = constr, n = 2 * num_years)
options <- list(dic = TRUE, mlik = TRUE, cpo = TRUE, 
                openmp.strategy = "default")
control.inla <- list(strategy = "adaptive", int.strategy = "auto")
temp <- smoothed_direct_bym2

# BYM2
bym2_model <- 
  INLA::inla.rgeneric.define(model = inla.rgeneric.bym2.model, 
                             n = num_years, Q_star = R_1_star + R_2_star, 
                             gamma_tilde = gamma_tilde, 
                             U_prec = pc.u, alpha_prec = pc.alpha, 
                             U_phi = pc.u.phi, alpha_phi = pc.alpha.phi)
mod.bym2 <- logit.est ~ f(region.struct, model = bym2_model, 
                          diagonal = 1e-06, extraconstr = constr, 
                          n = 2 * num_years)

#### Fit models ####
# 
for(i in 1:num_prec_settings){
    for(j in 1:num_par_settings){
        for(k in 1:num_reps){
            set.seed(37 * i * j * k)
            print(paste(i, j, k))
            load(paste0("../Data/dat_prec_", i, "_par_", j, "_rep_", k))
            temp$newdata[1:num_years, 
                          c("logit.est", "var.est", "logit.prec")] <- 
                dat[, c("logit.est", "var.est", "logit.prec")]
            fit <- INLA::inla(mod, family = "gaussian", 
                              control.compute = options, data = temp$newdata, 
                              control.predictor = list(compute = TRUE), 
                              control.family = 
                                  list(hyper = list(prec = 
                                                        list(initial = log(1), 
                                                             fixed = TRUE))), 
                              scale = temp$newdata$logit.prec, 
                              lincomb = temp$lincombs.fit, 
                              control.inla = control.inla, verbose = FALSE)
            save(fit, 
                 file = paste0("../Results/Fits_proposed/fit_prec_", i, 
                               "_par_", j, "_rep_", k))
            
            fit_bym2 <- INLA::inla(mod.bym2, family = "gaussian", 
                                   control.compute = options, 
                                   data = temp$newdata, 
                                   control.predictor = list(compute = TRUE), 
                                   control.family = 
                                     list(hyper = list(prec = 
                                                         list(initial = log(1), 
                                                              fixed = TRUE))), 
                                   scale = temp$newdata$logit.prec, 
                                   lincomb = temp$lincombs.fit, 
                                   control.inla = control.inla, verbose = FALSE)
            
            # dat_bym2format <- dat %>%
            #     mutate(region = years, years = 0, 
            #            region_num = as.numeric(region))
            # smoothed_direct_fit <- 
            #     smoothDirect_extra(data = dat_bym2format, 
            #                        Amat = amat, time.model = NULL,
            #                        year_label = "0", year_range = c(0), 
            #                        is.yearly = FALSE, m = 1, 
            #                        verbose = FALSE)
            save(fit_bym2, 
                 file = paste0("../Results/Fits_smoothed_direct/fit_prec_", i, 
                               "_par_", j, "_rep_", k))
        }
    }
}