#########################################
############ 2 Clusters #################
#########################################
library(doSNOW)
library(doRNG)
set.seed(20190327)
setwd("../data")
U <-read.table("DNAm_pre.txt",header =T)
U <- U[,-1:-2]
X <- read.csv("RNA_500.csv")
X <- X[,-1]
X <- as.matrix(X)

G <- dim(X)[1]
N <- dim(X)[2]
J <- dim(U)[1]
K <- 4 #K cluster

Z_o_t <- matrix(rnorm(N*(K-1),0,1),N,K-1)
Z_t <- cbind(rep(1,N),Z_o_t)

beta_t <- matrix(rnorm(G*K,0,1),G,K)

alpha_t <- matrix(rnorm(J*K,0,1),J,K)

Y_t <- as.matrix(X)
ind_zero <- (X==0)
Y_t[ind_zero] <- sample(c(1,2), sum(ind_zero), replace=TRUE)

lambda_t <- matrix(c(rep(0,G),rep(-1,G)),G,2)

numCores <- 5
cl <- makeCluster(numCores)
registerDoSNOW(cl)
Gam_t <- update_gam(Z_t, beta_t, Y_t)
Omega_t <- update_omega(Z_t, alpha_t, U)
stopCluster(cl)

#number of iterations
num_iter <- 5000

#collect posterior samples after burn-in
collection_Gam <- NULL
collection_Omega <- NULL
collection_z <- array(NA, dim=c(num_iter/2, N, K)) 
collection_y <- array(NA, dim=c(num_iter/2, G, N))
collection_lambda <- array(NA, dim=c(num_iter/2, G, 2))
collection_alpha <- array(NA, dim=c(num_iter/2, J, K)) 
collection_beta <- array(NA, dim=c(num_iter/2, G, K)) 

for(t in 1:num_iter){
  Y_t <- update_y(Y_t, X, ind_zero, lambda_t, Z_t, Gam_t, beta_t)
  
  numCores <- 5
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)
  
  beta_t <- update_beta(Z_t, Gam_t, beta_t, Y_t, K)
  alpha_t <- update_alpha(Z_t, Omega_t, alpha_t, U, K)
  Z_t <- update_z(Z_t, Gam_t, beta_t, Y_t, Omega_t, alpha_t, U)
  lambda_t <- update_lambda(Y_t, X, lambda_t)
  Gam_t <- update_gam(Z_t, beta_t, Y_t)
  Omega_t <- update_omega(Z_t, alpha_t, U)
  
  stopCluster(cl)
  
  if(t > num_iter/2){
    collection_Gam <- rbind(collection_Gam, Gam_t)
    collection_Omega <- rbind(collection_Omega, Omega_t)
    collection_z[t-num_iter/2,,] <- Z_t
    collection_y[t-num_iter/2,,] <- Y_t
    collection_lambda[t-num_iter/2,,] <- lambda_t
    collection_alpha[t-num_iter/2,,] <- alpha_t
    collection_beta[t-num_iter/2,,] <- beta_t
  }
}

## Inference

lambda_est <- apply(collection_lambda, c(2,3), mean)
alpha_est <- apply(collection_alpha, c(2,3), mean)
beta_est <- apply(collection_beta, c(2,3), mean)
Z_est <- apply(collection_z, c(2,3), mean)

getMode <- function(x){
  return(as.numeric(names(table(x))[table(x) == max(table(x))][1]))
}

Gam_est <- apply(collection_Gam, 2, getMode)
Omega_est <- apply(collection_Omega, 2, getMode)
Y_est <- apply(collection_y,c(2,3),getMode)

BIC_k(X, U, Y_est, Gam_est, Omega_est, alpha_est, beta_est, Z_est, lambda_est, K)