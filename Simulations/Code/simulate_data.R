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
prec_settings <- data.frame(rbind(rep(150, num_years),
                                  rep(300, num_years),
                                  rep(75, num_years)))
colnames(prec_settings) <- 1:num_years

#### Create structure matrices####
conflict_years <- 9:15
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
R_1_star_hat <- R_1_star[2:num_years, 2:num_years]
R_2_star_hat <- R_2_star[2:num_years, 2:num_years]
eps <- eigen(solve(R_1_star_hat + R_2_star_hat) %*% R_2_star_hat)$values
gamma_tilde <- c(1 / eigen(R_1_star + R_2_star)$values[1:(num_years - 1)], 0)
save(R_1_star, R_2_star, eps, gamma_tilde,
     file = "../Data/structure_matrices.RData")

#### Simulate data ####
data_template <- data.frame(region = "All", years = 1:num_years, logit.est = 0,
                            var.est = 0, logit.prec = 0, region_num = 0,
                            survey = NA, mean = NA, lower = NA, upper = NA,
                            risk = 0)
for(i in 1:num_prec_settings){
    for(j in 1:num_par_settings){
        for(k in 1:num_reps){
            set.seed(37 * i * j * k)
            print(paste(i, j, k))
            dat <- data_template
            dat$logit.prec <- unlist(prec_settings[i, ])
            dat$var.est <- 1 / dat$logit.prec
            
            for(l in 1:num_years){
                if(l %in% conflict_years){
                  if(!is.na(par_settings$mu_2[j])){
                    trend <- par_settings$mu_2[j]
                  }
                  else{
                    trend <- triangle[l]
                  }
                  trend <-
                    dat$risk[l] <- trend + 
                        rnorm(1, mean = 0, sd = 1 / sqrt(par_settings$tau_2[j]))
                        
                }
                else{
                  if(!is.na(par_settings$mu_1[j])){
                    trend <- par_settings$mu_1[j]
                  }
                  else{
                    trend <- triangle[l]
                  }
                    dat$risk[l] <- trend + 
                        rnorm(1, mean = 0, sd = 1 / sqrt(par_settings$tau_1[j]))
                    
                }
            }
            dat$logit.est <- rnorm(num_years, mean = dat$risk,
                                   sd = sqrt(dat$var.est))
            
            save(dat, 
                 file = paste0("../Data/dat_prec_", i, "_par_", j, "_rep_", k))
        }
    }
}
