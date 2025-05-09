library(glmnet)
library(iterators)
library(foreach)
library(doParallel)
library(fastclime)
library(MASS)
library(Matrix)
library(mvtnorm)
generate_model_1 <- function(p, n1, n2) {
  Omega_I <- matrix(0, nrow = p, ncol = p)
  diag(Omega_I) <- 1
  D <- matrix(0, nrow = p, ncol = p)
  for (i in 1:(p - 1)) {
    D[i, i + 1] <- 0.3
    D[i + 1, i] <- 0.3
  }
  Omega_II <- Omega_I + D
  data_I <- MASS::mvrnorm(n1, mu = rep(0, p), Sigma = solve(Omega_I))
  data_II <- MASS::mvrnorm(n2, mu = rep(0, p), Sigma = solve(Omega_II))
  return(list(Omega_I = Omega_I, Omega_II = Omega_II, data_I = data_I, data_II = data_II))
}
generate_model_2 <- function(p, n1, n2) {
  Omega_I <- diag(1, p)
  D=matrix(0,p,p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      D[i,j]=0.3*rbinom(1,1,5/p)
    }
  }
  D=D+t(D)
  Omega_II <- Omega_I + D
  data_I <- MASS::mvrnorm(n1, mu = rep(0, p), Sigma = solve(Omega_I))
  data_II <- MASS::mvrnorm(n2, mu = rep(0, p), Sigma = solve(Omega_II))
  return(list(Omega_I = Omega_I, Omega_II = Omega_II, data_I = data_I, data_II = data_II))
}
cross_validation <- function(X, Y, t) {
  lambda1_values <- 10^seq(-5, 2, length = 50 )
  ratio_values <- 10^seq(-2, 2, length = 20)
  # ratio_values <- seq(0.1, 2, length = 10)
  results <- data.frame(ratio = numeric(), lambda1 = numeric(), mse = numeric())
  for (ratio in ratio_values){
    Z=c(X[,t]/sqrt(2*nrow(X)),Y[,t]/sqrt(2*nrow(Y)))
    Z=as.matrix(Z)
    X_upper <- cbind(X[,-t]/sqrt(2*nrow(X)), Y[,-t]/sqrt(2*nrow(Y)))
    X_lower <- cbind(matrix(0, nrow = nrow(Y[,-t]), ncol = ncol(X[,-t])), ratio * Y[,-t]/sqrt(2*nrow(Y)))
    tilde_X <- rbind(X_upper, X_lower)
    cv_fit <- cv.glmnet(tilde_X, Z, alpha = 1, lambda = lambda1_values, family="gaussian",type.measure = "mse", nfolds = 5)
    mse <- min(cv_fit$cvm)
    best_lambda1 <- cv_fit$lambda.min
    results <- rbind(results, data.frame(ratio = ratio, lambda1 = best_lambda1, mse = mse))
  }
  best_result <- results[which.min(results$mse), ]
  return(best_result)
}
DDME <- function(X_data,Y_data,M) {
n1=dim(X_data)[1]
n2=dim(Y_data)[1]
p=dim(X_data)[2]
upper_index = ceiling(n1/M)
lower_index = floor(n1/M)
truncation_index = n1 - lower_index * M
# M Number of machines
#p Number of variables
#n1 and n2 are Number of samples x and y
nworkers <- detectCores() - 1  
cl <- makeCluster(nworkers)
registerDoParallel(cl)
############### our method
local_Result1 = foreach(i = 1:M, .packages = c("glmnet","fastclime"),.export = c("cross_validation")) %dopar%{
    if(i <= truncation_index) {
      X_temp=matrix(0,upper_index,p)
      Y_temp=matrix(0,upper_index,p)
      X_temp=X[(1+upper_index*(i-1)):(upper_index*i), ]
      Y_temp=Y[(1+upper_index*(i-1)):(upper_index*i), ]
    }else{
      X_temp=matrix(0,lower_index,p)
      Y_temp=matrix(0,lower_index,p)
      X_temp=X_data[(truncation_index * upper_index + 1 + lower_index * ( i - 1 - truncation_index)):(truncation_index * upper_index + lower_index * ( i - truncation_index)), ]
      Y_temp=Y_data[(truncation_index * upper_index + 1 + lower_index * ( i - 1 - truncation_index)):(truncation_index * upper_index + lower_index * ( i - truncation_index)), ]
    }
Beta1=matrix(0,p,p)###beta1_debias
Beta2=matrix(0,p,p)###beta2_debias
beta1=matrix(0,p,p)###beta1_non_debias
beta2=matrix(0,p,p)###beta2_non_debias
    r=rep(0,p)
    for (t in 1:p){
      best_lambdas <- cross_validation(X_temp, Y_temp, t)
      Z=c(X_temp[,t]/sqrt(2*nrow(X_temp)),Y_temp[,t]/sqrt(2*nrow(Y_temp)))
      Z=as.matrix(Z)
      X_upper <- cbind(X_temp[,-t]/sqrt(2*nrow(X_temp)), Y_temp[,-t]/sqrt(2*nrow(Y_temp)))
      X_lower <- cbind(matrix(0, nrow = nrow(Y_temp[,-t]), ncol = ncol(X_temp[,-t])),  (best_lambdas$ratio) * Y_temp[,-t]/sqrt(2*nrow(Y_temp)))
      tilde_X <- rbind(X_upper, X_lower)
      beta_estimates <- glmnet(tilde_X, Z, family="gaussian", alpha=1, lambda = best_lambdas$lambda1)
      beta_1=beta_estimates$beta[1:(p-1)]
      beta_d=beta_estimates$beta[p:(2*p-2)]*(best_lambdas$ratio)
      k1=t(X_temp[,-t])%*%(X_temp[,t]-X_temp[,-t]%*%beta_1)/nrow(X_temp)+t(Y_temp[,-t])%*%(Y_temp[,t]-Y_temp[,-t]%*%(beta_1+beta_d))/nrow(Y_temp)
      k2=t(Y_temp[,-t])%*%(Y_temp[,t]-Y_temp[,-t]%*%(beta_1+beta_d))/nrow(Y_temp)
      Sigma1=cov(X_temp[,-t])
      Sigma2=cov(Y_temp[,-t])
      Clime_X_temp=fastclime(Sigma1,lambda.min=0.01)
      M1=Clime_X_temp$icovlist
      index1=0;
      for( y in 2:length(M1)){
        if (norm(M1[[y-1]]%*%Sigma1-diag(1,p-1,p-1),type = "M") <= norm(M1[[y]]%*%Sigma1-diag(1,p-1,p-1),type = "M")){
          index1=y-1;
        }
        else {index1=y}
      }
      Clime_Y_temp=fastclime(Sigma2,lambda.min=0.01)
      M2=Clime_Y_temp$icovlist
      index2=0;
      for( y in 2:length(M2)){
        if (norm(M2[[y-1]]%*% Sigma2-diag(1,p-1,p-1),type = "M") <= norm(M2[[y]]%*% Sigma2-diag(1,p-1,p-1),type = "M")){
          index2=y-1;
        }
        else {index2=y}
      }
      beta1[t,-t] = beta_1
      beta2[t,-t] = beta_d+beta_1
      Beta1[t,-t] = beta_1+M1[[index1]]%*%(k1-k2)
      Beta2[t,-t] = beta_d+beta_1+M2[[index2]]%*%k2
      r[t]=(sum((X_temp[,t]-X_temp[,-t]%*%Beta1[t,-t])^2)+sum((Y_temp[,t]-Y_temp[,-t]%*%Beta2[t,-t])^2))/(n1/M+n2/M)
    }
    list(Beta1 = Beta1, Beta2 = Beta2, beta1 = beta1, beta2 = beta2, r=r)
  }
stopCluster(cl)

Beta1_sum=matrix(0,p,p)###sum_Beta1
Beta2_sum=matrix(0,p,p)###sum_Beta2
Beta1_agg=matrix(0,p,p)###agg_Beta1
Betad_agg=matrix(0,p,p)###agg_Betad
for (m in 1:M){
    Beta1_sum=Beta1_sum+local_Result1[[m]]$Beta1
    Beta2_sum=Beta2_sum+local_Result1[[m]]$Beta2
  } 
pi1=2*sqrt(log(p)/n1)
pi2=M*(log(p)/n2)
pi_agg=max(pi1,pi2)
Beta1_agg=Beta1_sum/M
Betad_agg=(Beta2_sum-Beta1_sum)/M
Beta1_agg=Beta1_agg*(!abs(Beta1_agg) <= pi_agg)
Betad_agg=Betad_agg*(!abs(Betad_agg) <= pi_agg)
cl <- makePSOCKcluster(nworkers)
clusterSetRNGStream(cl)
registerDoParallel(cl)
local_Result2 = foreach(i = 1:M, .packages = c("glmnet")) %dopar%{
  if(i <= truncation_index) {
    X_temp=matrix(0,upper_index,p)
    Y_temp=matrix(0,upper_index,p)
    X_temp=X_data[(1+upper_index*(i-1)):(upper_index*i), ]
    Y_temp=Y_data[(1+upper_index*(i-1)):(upper_index*i), ]
  }else{
    X_temp=matrix(0,lower_index,p)
    Y_temp=matrix(0,lower_index,p)
    X_temp=X_data[(truncation_index * upper_index + 1 + lower_index * ( i - 1 - truncation_index)):(truncation_index * upper_index + lower_index * ( i - truncation_index)), ]
    Y_temp=Y_data[(truncation_index * upper_index + 1 + lower_index * ( i - 1 - truncation_index)):(truncation_index * upper_index + lower_index * ( i - truncation_index)), ]
  }
  r=rep(0,p)
  for (t in 1:p){
    r[t]=(sum((X_temp[,t]-X_temp[,-t]%*%Beta1_agg[t,-t])^2)+sum((Y_temp[,t]-Y_temp[,-t]%*%(Beta1_agg[t,-t]+Betad_agg[t,-t]))^2))/(n1/M+n2/M)
    }
  list(r = r)
}
stopCluster(cl)
r_sum=rep(0,p)
for (m in 1:M){
  r_sum=r_sum+local_Result2[[m]]$r
} 
r_agg=rep(0,p)
r_agg=r_sum/M
D_est=matrix(0,p,p)
for (t in 1:p) {
  D_est[t, ] <- Betad_agg[t, ] / r_agg[t]
}
D_est=-(D_est+t(D_est))/2
return(D_est=D_est)
}




