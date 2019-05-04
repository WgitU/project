#########################################
#############Function parallel###########
#########################################
update_beta <- function(Z_t, Gam_t, beta_t, Y_t, K, beta_0=rep(0,K), Sigma2_beta=1, sig_beta=0.1) {
  G <- dim(Y_t)[1]
  N <- dim(Y_t)[2]
  beta_star <- matrix(rnorm(G*K, beta_t, sd = sig_beta),G,K)
  
  
  beta_t <- foreach(g=1:G,.combine=rbind) %dorng% {
    Gam_g <- diag(K)*Gam_t[g]
    Gam_g[1,1] <- 1
    part1 <- 0
    for (i in 1:N) {
      part1 = part1 + Z_t[i,]%*%Gam_g%*%(beta_star[g,]-beta_t[g,])*Y_t[g,i] - exp(Z_t[i,]%*%Gam_g%*%beta_star[g,]) + exp(Z_t[i,]%*%Gam_g%*%beta_t[g,])
    }
    part2 <- (beta_star[g,]%*%beta_star[g,] - beta_t[g,]%*%beta_t[g,] - 2*(beta_star[g,]-beta_t[g,])%*%beta_0)/Sigma2_beta/2
    log_r <- part1 -part2
    r <- exp(log_r)
    
    if (runif(1) < r) {
      beta_star[g,]
    } else {
      beta_t[g,]
    }
  }
  
  return(beta_t)
}


update_alpha <- function(Z_t, Omega_t, alpha_t, U, K, alpha_0=rep(0,K), Sigma2_alpha=1, sig_alpha=0.1) {
  N <- dim(U)[2]
  J <- dim(U)[1]
  
  alpha_star <- matrix(rnorm(J*K, alpha_t, sd = sig_alpha),J,K)
  
  alpha_t <- foreach(j=1:J,.combine=rbind) %dorng% {
    omega_g <- diag(K)*Omega_t[j]
    omega_g[1,1] <- 1
    part1 <- 0
    for (i in 1:N) {
      part1 = part1 + Z_t[i,]%*%omega_g%*%(alpha_star[j,]-alpha_t[j,])*U[j,i] - exp(Z_t[i,]%*%omega_g%*%alpha_star[j,]) + exp(Z_t[i,]%*%omega_g%*%alpha_t[j,])
    }
    part2 <- (alpha_star[j,]%*%alpha_star[j,] - alpha_t[j,]%*%alpha_t[j,] - 2*(alpha_star[j,]-alpha_t[j,])%*%alpha_0)/Sigma2_alpha/2
    log_r <- part1 - part2
    r <- exp(log_r)
    
    if (runif(1) < r) {
      alpha_star[j,]
    } else {
      alpha_t[j,] 
    }
  }
  
  return(alpha_t)
}


update_gam <- function(Z_t, beta_t, Y_t, q_gam=0.5) {
  G <- dim(Y_t)[1]
  N <- dim(Y_t)[2]
  
  Gam <- foreach(g=1:G,.combine=c) %dorng% {
    pr1 <- 0
    pr2 <- 0
    for (i in 1:N) {
      pr1 <- pr1 + Z_t[i,]%*%beta_t[g,]*Y_t[g,i] - exp(Z_t[i,]%*%beta_t[g,])
      pr2 <- pr2 + beta_t[g,1]*Y_t[g,i] - exp(beta_t[g,1])
    }
    pr1 <- pr1 + log(q_gam)
    pr2 <- pr2 + log(1-q_gam)
    pr <- 1/(1+exp(pr2-pr1))
    if (runif(1) < pr) {
      1
    } else {
      0
    }
  }
  
  return(Gam)
}


update_omega <- function(Z_t, alpha_t, U, q_omega=0.5) {
  N <- dim(U)[2]
  J <- dim(U)[1]
  
  Omega <- foreach(j=1:J,.combine=c) %dorng% {
    pr1 <- 0
    pr2 <- 0
    for (i in 1:N) {
      pr1 <- pr1 + Z_t[i,]%*%alpha_t[j,]*U[j,i] - exp(Z_t[i,]%*%alpha_t[j,])
      pr2 <- pr2 + alpha_t[j,1]*U[j,i] - exp(alpha_t[j,1])
    }
    pr1 <- pr1 + log(q_omega)
    pr2 <- pr2 + log(1-q_omega)
    pr <- 1/(1+exp(pr2-pr1))
    if (runif(1) < pr) {
      1
    } else {
      0
    }
  }
  
  return(Omega)
}


update_z <- function(Z_t, Gam_t, beta_t, Y_t, Omega_t, alpha_t, U, sig_z=0.1) {
  G <- dim(Y_t)[1]
  N <- dim(Y_t)[2]
  J <- dim(U)[1]
  K <- dim(Z_t)[2]
  
  z_o_star <- matrix(rnorm(N*(K-1), Z_t[,2:K], sd = sig_z),N,K-1)
  z_star <- cbind(rep(1,N),z_o_star)
  
  Z_t_new <- foreach(i=1:N,.combine=rbind) %dorng% {
    part1 <- 0
    for (g in 1:G) {
      Gam_g <- diag(K)*Gam_t[g]
      Gam_g[1,1] <- 1
      part1 <- part1 + (z_star[i,] - Z_t[i,])%*%Gam_g%*%beta_t[g,]*Y_t[g,i] - exp(z_star[i,]%*%Gam_g%*%beta_t[g,]) + exp(Z_t[i,]%*%Gam_g%*%beta_t[g,])
    }
    part2 <- 0
    for (j in 1:J) {
      omega_g <- diag(K)*Omega_t[j]
      omega_g[1,1] <- 1
      part2 <- part2 + (z_star[i,] - Z_t[i,])%*%omega_g%*%alpha_t[j,]*U[j,i] - exp(z_star[i,]%*%omega_g%*%alpha_t[j,]) + exp(Z_t[i,]%*%omega_g%*%alpha_t[j,])
    }
    part3 <- (z_star[i,]%*%z_star[i,] - Z_t[i,]%*%Z_t[i,])/2
    
    log_r <- part1 + part2 - part3
    r <- exp(log_r)
    if (runif(1) < r) {
      z_star[i,]
    } else {
      Z_t[i,]
    }
  }
  
  return(Z_t_new)
}


update_lambda <- function(Y_t, X, lambda_t, sig_lam = c(1,0.1), lambda_0=c(0,-1), Sigma2_lambda=c(1,0.1)) {
  ##library(BayesLogit)
  G <- dim(Y_t)[1]
  N <- dim(Y_t)[2]
  
  lambda_t_new <- foreach(g=1:G,.combine=rbind) %dorng% {
    lambda1_g_star <- rnorm(1,lambda_t[g,2], sd = sig_lam[2])
    if (lambda1_g_star < 0) {
      lambda0_g_star <- rnorm(1,lambda_t[g,1], sd = sig_lam[1])
      
      log2_y_g <- log2(Y_t[g, ] + 1)
      
      ind_x0_yn0 <- (X[g, ] == 0 && Y_t[g, ] > 0)
      ind_n0 = (X[g, ] > 0)
      
      part1 <- prod(pnorm(lambda0_g_star + lambda1_g_star * log2_y_g[ind_x0_yn0])/
                      pnorm(lambda_t[g,1] + lambda_t[g,2] * log2_y_g[ind_x0_yn0]))
      part2 <- prod((1 - pnorm(lambda0_g_star + lambda1_g_star * log2_y_g[ind_n0]))/
                      (1 - pnorm(lambda_t[g,1] + lambda_t[g,2] * log2_y_g[ind_n0])))
      part3 <- exp(-((lambda0_g_star - lambda_t[g,1])*(lambda0_g_star + lambda_t[g,1] - 2*lambda_0[1])/Sigma2_lambda[1] +
                       (lambda1_g_star - lambda_t[g,2])*(lambda1_g_star + lambda_t[g,2] - 2*lambda_0[2])/Sigma2_lambda[2])/2)
      r <- part1 *part2 *part3
      
      if (runif(1) < r) {
        lambda_g <- c(lambda0_g_star, lambda1_g_star)
      } else {
        lambda_t[g,]
      }
    } else {
      lambda_t[g,]
    }
  }
  
  return(lambda_t_new)
}


update_y <- function(Y_t, X, ind_zero, lambda_t, Z_t, Gam_t, beta_t, radius=1) {
  G <- dim(Y_t)[1]
  N <- dim(Y_t)[2]
  y_star <- Y_t
  tmp <- sample(seq(-radius, radius, by=1), sum(ind_zero), replace=TRUE)
  y_star[ind_zero] <- Y_t[ind_zero] + tmp
  
  r <- matrix(0, G, N)
  gam_beta <- NULL
  for (g in 1:G) {
    gam_beta_tmp <- beta_t[g,]
    gam_beta_tmp[2:K] <- gam_beta_tmp[2:K]*Gam_t[g]
    gam_beta <- cbind(gam_beta, gam_beta_tmp)
  }
  theta <- t(Z_t%*%gam_beta)
  
  lambda_tmp <- matrix(lambda_t[,2], G, N)
  
  ind_tmp <- (y_star >= 0 & ind_zero == T & Y_t > 0)
  ind_0tmp <- (y_star >= 0 & ind_zero == T & Y_t == 0)
  part1 <- theta[ind_tmp]*(y_star[ind_tmp] - Y_t[ind_tmp])
  part2 <- lfactorial(Y_t[ind_tmp]) - lfactorial(y_star[ind_tmp])
  
  lambda0_tmp <- matrix(lambda_t[,1], G, N)
  lambda1_tmp <- matrix(lambda_t[,2], G, N)
  part3 <- pnorm(lambda0_tmp[ind_tmp] + lambda1_tmp[ind_tmp] * log2(y_star[ind_tmp]+1)) /
    pnorm(lambda0_tmp[ind_tmp] + lambda1_tmp[ind_tmp] * log2(Y_t[ind_tmp]+1))
  
  r[ind_tmp] <- exp(part1 + part2) * part3
  
  part4 <- theta[ind_0tmp]*(y_star[ind_0tmp] - Y_t[ind_0tmp])
  part5 <- lfactorial(Y_t[ind_0tmp]) - lfactorial(y_star[ind_0tmp])
  r[ind_0tmp] <- exp(part4 + part5)
  
  tmp_r <- matrix(runif(G*N), G, N)
  Y_t[tmp_r <= r] <- y_star[tmp_r <= r]
  
  return(Y_t)
}


BIC_k <- function(X, U, Y, Gam, Omega, alpha, beta, Z, lambda, K) {
  G <- dim(X)[1]
  N <- dim(X)[2]
  J <- dim(U)[1]
  
  omega_alp <- NULL
  for (j in 1:J) {
    tmp <- alpha[j,]
    tmp[2:K] <- tmp[2:K]*Omega[j]
    omega_alp <- cbind(omega_alp, tmp)
  }
  para_u <- t(Z%*%omega_alp)
  log_u <- sum(para_u * U - exp(para_u) - lfactorial(U))
  
  Gam_beta <- NULL
  for (g in 1:G) {
    tmp <- beta[g,]
    tmp[2:K] <- tmp[2:K]*Gam[g]
    Gam_beta <- cbind(Gam_beta, tmp)
  }
  para_y <- t(Z%*%Gam_beta)
  
  ind_x <- (X > 0)
  ind_x0_y <- (X == 0 & Y >0)
  ind_0 <- (X == 0 & Y == 0)
  
  para_y[ind_0] <- -Inf
  
  S <- 1000
  
  sample_y <- array(rpois(S*G*N, exp(para_y)), dim = c(G,N,S))
  
  pr1 <- pnorm(lambda[,1] + lambda[,2]*log2(sample_y+1))
  
  pr2 <- apply(pr1, c(1,2), mean) 
  
  log_x <- sum(log(1-pr2[ind_x])) + sum(log(pr2[ind_x0_y])) 
  
  BIC_res <- -2*(log_x+log_u) + log(N*(G+J))*(G*(K+3)+J*(K+1)+N*(K-1))
  
  return(BIC_res)
}