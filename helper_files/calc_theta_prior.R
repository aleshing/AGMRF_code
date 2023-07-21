calc_theta_prior_point <- function(theta, eigen, lambda){
    n <- length(eigen)
    temp1 <- 1 + (theta - 1) * eigen
    d2 <- sum(1 / temp1) - n + sum(log(temp1))
    jacob <- sum((eigen ^ 2) / (temp1 ^ 2)) 
    
    log(lambda) - log(2) - log(d2) /2 - lambda * sqrt(d2) + log(1 - theta) + 
      log(jacob)
}

calc_theta_prior <- function(eigen, U, alpha){
    n <- length(eigen)
    temp2 <- 1 + (U - 1) * eigen
    dU <- sqrt(sum(1 / temp2) - n + sum(log(temp2)))
    lambda <- -log(alpha) / dU
    
    thetas <- seq(0.0001, 0.9999, by = 0.0001)
    log_prior <- rep(0, length(thetas))
    for(i in 1:length(thetas)){
        log_prior[i] <- calc_theta_prior_point(thetas[i], eigen, lambda)
    }
    return(data.frame(theta = thetas, log_density = log_prior))
}
