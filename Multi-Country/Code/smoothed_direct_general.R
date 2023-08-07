library(SUMMER)
library(readstata13)
library(dplyr)
library(rgdal)

#### Load helper functions and data ####
load(file = "../Data/generated_data/direct_estimates.RData")
load(file = "../Data/generated_data/admin_info.RData")
source("../../helper_files/calc_theta_prior.R")
source("../../helper_files/adaptive_bym2_rgeneric.R")
source("../../helper_files/multicountry_rgeneric.R")
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

#### Get smoothed direct estimates, more general adaptive bym2 ####
# Create structure matrices
num_region <- nrow(country_admin_table)
countries <- unique(country_admin_table$Country)
num_countries <- length(countries)
W_list <- vector(mode = "list", length = num_countries + 1)
for(m in 1:(num_countries + 1)){
  W_list[[m]] <- amat
}
for(i in 1:num_region){
  for(j in i:num_region){
    if(amat[i, j] == 1){
      
      if(country_admin_table[i, "Country"] == 
         country_admin_table[j, "Country"]){
        W_list[[num_countries + 1]][i, j] = 0
        W_list[[num_countries + 1]][j, i] = 0
        for(m in 1:num_countries){
          if(country_admin_table[i, "Country"] != countries[m]){
            W_list[[m]][i, j] = 0
            W_list[[m]][j, i] = 0
          }
        }
      }
      else{
        for(m in 1:num_countries){
          W_list[[m]][i, j] = 0
          W_list[[m]][j, i] = 0
        }
      }
    }
  }
}

D_list <- vector(mode = "list", length = num_countries + 1)
R_list <- vector(mode = "list", length = num_countries + 1)
for(m in 1:(num_countries + 1)){
  D_list[[m]] <- diag(rowSums(W_list[[m]]))
  R_list[[m]] <- D_list[[m]] - W_list[[m]]
}

R_sum <- R_list[[num_countries + 1]]
for(m in 1:num_countries){
  R_sum <- R_sum + R_list[[m]]
}
scaled_Q <- INLA:::inla.scale.model.bym.internal(R_sum, 
                                                 adjust.for.con.comp = TRUE)$Q
gv <- scaled_Q[1, 1] / (R_sum)[1, 1]

R_star_list <- vector(mode = "list", length = num_countries + 1)
for(m in 1:(num_countries + 1)){
  R_star_list[[m]] <- gv * (D_list[[m]] - W_list[[m]])
}

vals <- (1:num_region)[-num_region]
R_1_star <- gv * (R_sum - R_list[[num_countries + 1]])
R_2_star <- R_star_list[[num_countries + 1]]
R_1_star_hat <- R_1_star[vals, vals]
R_2_star_hat <- R_2_star[vals, vals]
eps <- eigen(solve(R_1_star_hat + R_2_star_hat) %*% R_2_star_hat)$values
gamma_tilde <- c(1 / eigen(R_1_star + R_2_star)$values[1:(num_region - 1)], 0)
save(R_star_list, eps, gamma_tilde,
     file = "../Data/generated_data/structure_matrices_general.RData")

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
  INLA::inla.rgeneric.define(model = inla.rgeneric.multicountry.model, 
                             n = num_region, 
                             M = num_countries,
                             R_star_list = R_star_list,
                             gamma_tilde = gamma_tilde, 
                             eps = eps, U_prec = pc.u, alpha_prec = pc.alpha, 
                             U_phi = pc.u.phi, alpha_phi = pc.alpha.phi,
                             U_theta = pc.u.theta, 
                             alpha_theta = pc.alpha.theta,
                             psi_prior_setting = 1)
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

smoothed_direct_adaptive_bym2_general <- smoothed_direct_adaptive_bym2

out_general <- cbind(out_adaptive_bym2, prior = "adaptive bym2 general")
save(out_general, smoothed_direct_adaptive_bym2_general, 
     file = "../Results/smoothed_direct_estimates_general.RData")
