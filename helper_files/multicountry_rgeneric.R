inla.rgeneric.multicountry.model = function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", 
            "quit"),
    theta = NULL)
{
  # Assume we are passed the following inputs
  # n: the number of regions
  # M: the number of countries 
  # R_star_list: a list of the the scaled structure matrices for the structured 
  # component, should be length M + 1
  # gamma_tilde: the inverse eigenvalues of Q_star = sum_{m=1}^{M+1} R_m_star
  # eps: the eigenvalues of 
  # (sum_{m=1}^{M+1} \hat{R}_m_star)^{-1}\hat{R}_{M+1}_star, where \hat{} 
  # denotes removing the first row and column of the matrix
  # U_prec, alpha_prec: U and alpha parameters for precision PC prior
  # U_phi, alpha_phi: U and alpha parameters for mixing parameter PC prior
  # U_theta, alpha_theta: U and alpha parameters for theta parameter PC prior
  # psi_prior_setting: If psi_prior_setting = 0, use U_prec, alpha_prec. Else, 
  # if psi_prior_setting = 1, use U = 1, alpha = 0.5
  
  # for reference and potential storage for objects to
  # cache, this is the environment of this function
  # which holds arguments passed as `...` in
  # `inla.rgeneric.define()`.
  envir = parent.env(environment())
  
  interpret.theta = function() {
    pars <- list(prec = exp(theta[1L]),
                 phi = 1 / (1 + exp(-theta[2L])),
                 thet = 1 / (1 + exp(-theta[3L])))
    if(M == 1){
      return(pars)
    }
    else{
      for(m in 2:M){
        # Is not having an L here OK?
        pars[[2 + m]] <- exp(theta[2 + m])
        names(pars)[2 + m] <- paste0("psi_", m)
      }
      return(pars)
    }
  }
  
  graph = function(){ return (Q()) }
  
  Q = function() { 
    p = interpret.theta()
    D = (1 / (1 - p$phi)) * diag(n) 
    R_sum <- R_star_list[[1]] + p$thet * R_star_list[[M + 1]]
    if(M > 1){
      for(m in 2:M){
        R_sum <- R_sum + p[[2 + m]] * R_star_list[[m]]
      }
    }
    QQ = rbind(cbind(p$prec * D, -sqrt(p$phi * p$prec) * D),
               cbind(-sqrt(p$phi * p$prec) * D, 
                     R_sum + p$phi * D))
    return (inla.as.sparse(QQ))
  }
  
  mu = function() { return(numeric(0)) } 
  
  log.norm.const = function() { return (numeric(0)) } 
  
  log.prior = function() { 
    p = interpret.theta()
    
    # Construct prior for the precision parameter
    lambda_prec = -log(alpha_prec) / U_prec
    prec_prior = log(lambda_prec) - log(2) - theta[1L] / 2 -
      lambda_prec / sqrt(p$prec)
    
    # Construct prior for psi
    if(psi_prior_setting == 0){
      psi_prior = 0
      for(m in 2:M){
        psi_prior = psi_prior + 
          log(lambda_prec) - log(2) - theta[2 + m] / 2 -
          lambda_prec / sqrt(p[[2 + m]])
      }
    }
    else{
      lambda_prec = -log(0.5) 
      psi_prior = 0
      for(m in 2:M){
        psi_prior = psi_prior + 
          log(lambda_prec) - log(2) - theta[2 + m] / 2 -
          lambda_prec / sqrt(p[[2 + m]])
      }
    }
    
    # Construct prior for the mixing parameter phi
    dU_phi = sqrt(U_phi * sum(gamma_tilde - 1) -
                    sum(log(1 + U_phi * (gamma_tilde - 1))))
    lambda_phi = -log(1 - alpha_phi) / dU_phi
    d2_phi = p$phi * sum(gamma_tilde - 1) - 
      sum(log(1 + p$phi * (gamma_tilde - 1)))
    phi_prior = log(lambda_phi) + 2 * theta[2L] -
      2 * log1p(exp(theta[2L])) - log(2) - log(d2_phi) / 2 -
      lambda_phi * sqrt(d2_phi) +
      log(abs(sum( (gamma_tilde - 1) ^ 2 /
                     (1 + exp(theta[2L]) * gamma_tilde))))
    
    # Construct prior for the theta parameter
    dU_theta <- sqrt(sum(1 / (1 + (U_theta - 1) * eps)) + 
                       sum(log(1 + (U_theta - 1) * eps)) - (n - 1))
    lambda_theta <- -log(alpha_theta) / dU_theta
    d2_theta <- sum(1 / (1 + (p$thet - 1) * eps)) + 
      sum(log(1 + (p$thet - 1) * eps)) - (n - 1)
    theta_prior <- log(lambda_theta) - 3 * log1p(exp(theta[3L])) - log(2) -
      log(d2_theta) / 2 + 
      log(sum((eps ^ 2) / ((1 + (p$thet - 1) * eps) ^ 2))) -
      lambda_theta * sqrt(d2_theta) + theta[3L]
    
    return(prec_prior + phi_prior + theta_prior + psi_prior)
    
  } 
  
  initial = function() { 
    init <- c(4, 0, 0)
    if(M == 1){
      return(init) 
    }
    else{
      for(m in 2:M){
        init[2 + m] <- 4
      }
      return(init)
    }
    
  }
  
  quit = function() { return (invisible()) }
  
  if (!length(theta)) theta = initial()
  val = do.call(match.arg(cmd), args = list())
  return (val) 
}