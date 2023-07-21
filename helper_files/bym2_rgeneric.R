inla.rgeneric.bym2.model = function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", 
            "quit"),
    theta = NULL)
{
    # Assume we are passed the following inputs
    # n: the number of subdivisions (i.e. time points or regions)
    # Q_star: the scaled structure matrix for the structured component
    # gamma_tilde: the inverse eigenvalues of Q_star
    # U_prec, alpha_prec: U and alpha parameters for precision PC prior
    # U_phi, alpha_phi: U and alpha parameters for mixing parameter PC prior
    
    # for reference and potential storage for objects to
    # cache, this is the environment of this function
    # which holds arguments passed as `...` in
    # `inla.rgeneric.define()`.
    envir = parent.env(environment())
    
    interpret.theta = function() {
        return(list(prec = exp(theta[1L]),
                    phi = 1 / (1 + exp(-theta[2L]))))
    }
    
    graph = function(){ return (Q()) }
    
    Q = function() { 
        p = interpret.theta()
        D = (1 / (1 - p$phi)) * diag(n) 
        QQ = rbind(cbind(p$prec * D, -sqrt(p$phi * p$prec) * D),
                   cbind(-sqrt(p$phi * p$prec) * D, Q_star + p$phi * D))
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
        
        # Construct prior for the mixing parameter phi
        dU = sqrt(U_phi * sum(gamma_tilde - 1) -
                      sum(log(1 + U_phi * (gamma_tilde - 1))))
        lambda_phi = -log(1 - alpha_phi) / dU
        d2 = p$phi * sum(gamma_tilde - 1) - 
            sum(log(1 + p$phi * (gamma_tilde - 1)))
        phi_prior = log(lambda_phi) + 2 * theta[2L] -
            2 * log1p(exp(theta[2L])) - log(2) - log(d2) / 2 -
            lambda_phi * sqrt(d2) +
            log(abs(sum( (gamma_tilde - 1) ^ 2 /
                             (1 + exp(theta[2L]) * gamma_tilde))))
        
        return(prec_prior + phi_prior)
            
    } 
    
    initial = function() { return(c(4, 0)) }
    
    quit = function() { return (invisible()) }
    
    if (!length(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())
    return (val) 
}