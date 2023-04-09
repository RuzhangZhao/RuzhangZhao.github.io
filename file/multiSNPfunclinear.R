pacman::p_load(inline,data.table,dplyr,speedglm,Rfast,feather,locfit,bigmemory,stringr,glmnet,PRROC)
library(inline)
library(data.table)
library(dplyr)
library(speedglm)
library(Rfast)
library(feather)
library(locfit)
library(bigmemory)
library(stringr)
library(glmnet)
library(PRROC)
library(ROCR)
library(expm)
U_func<-function(
  UKBB_pop,
  beta,
  theta_UKBB_GPC,
  study_info){
  N_Pop<-nrow(UKBB_pop)
  expit_beta<-(UKBB_pop[,-1]%*%beta)
  u1<-crossprod((UKBB_pop[,1] -expit_beta),UKBB_pop[,-1])
  u2<-c()
  u2_part1<-c(expit_beta)%*%UKBB_pop[,var_SNP]
  u2_part2<-sapply(1:N_SNP, function(snp_id){
    u2_id<-UKBB_pop[,c("V",var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)
    c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
  })
  u2<-u2_part1 - u2_part2
  u<-c(u1,u2)*(1/N_Pop)
}

utcu_C<-function(u,C){
  c(u%*%C%*%u)
}
utcu<-function(UKBB_pop,beta,theta_UKBB_GPC,study_info,C){
  u<-U_func(UKBB_pop,beta,theta_UKBB_GPC,study_info)
  utcu_C(u,C)
}


grad_U_wrt_beta_func<-function(
  UKBB_pop,
  beta){
  N_Pop<-nrow(UKBB_pop)
  U1_beta_gradient<-crossprod(UKBB_pop[,-1],UKBB_pop[,-1])
  U2_beta_gradient<-crossprod(UKBB_pop[,var_SNP],UKBB_pop[,-1])
  rbind(U1_beta_gradient,U2_beta_gradient)*(1/N_Pop)
}

U3_func<-function(
  UKBB_pop,
  theta_UKBB_GPC){
  expit_PC<-(UKBB_pop[,c("V",var_GPC)]%*%theta_UKBB_GPC)
  # UKBB_pop[,1] is Y(phenotype)
  crossprod((UKBB_pop[,"Y"]-expit_PC),UKBB_pop[,c("V",var_GPC)])
}

### Generate coefficients of GPCs. 
inv_grad_U3_wrt_theta1_func<-function(
  UKBB_pop,
  theta_UKBB_GPC){
  N_Pop<-nrow(UKBB_pop)
  mat<--(1/N_Pop)*crossprod(UKBB_pop[,c("V",var_GPC)],UKBB_pop[,c("V",var_GPC)])
  #-solve(mat,tol=1e-60)
  -ginv(mat)
}

grad_U1_wrt_theta_func<-function(
  len_U1,
  len_theta){
  matrix(0,nrow = len_U1,ncol = len_theta)
}

grad_U2_wrt_theta_func<-function(
  UKBB_pop,
  theta_UKBB_GPC,
  study_info){
  N_Pop<-nrow(UKBB_pop)
  u2_theta<-sapply(1:N_SNP, function(snp_id){
    dexpit_id<-c(UKBB_pop[,paste0("SNP",snp_id)])
    u2id_theta1_gradient<-dexpit_id%*%UKBB_pop[,c("V",var_GPC)]
    u2id_gammaid_gradient<-dexpit_id%*%c(UKBB_pop[,paste0("SNP",snp_id)])
    c(u2id_theta1_gradient,u2id_gammaid_gradient)
  })
  #u2_theta is (1+len_GPC+1)*N_SNP
  N_snp_id<-nrow(u2_theta)
  # Gradient is N_equation*N_variables
  -(1/N_Pop)*t(rbind(u2_theta[-N_snp_id,],diag(u2_theta[N_snp_id,])))
}

var_U_beta_theta_func<-function(
  UKBB_pop,
  beta,
  theta_UKBB_GPC,
  study_info){
  N_Pop<-nrow(UKBB_pop)
  expit_beta<-(UKBB_pop[,-1]%*%beta)
  var_11<-crossprod(UKBB_pop[,-1]*c(expit_beta-UKBB_pop[,1]),UKBB_pop[,-1]*c(expit_beta-UKBB_pop[,1]))
  u2_theta_coef<-sapply(1:N_SNP, function(snp_id){
    expit_id<-c((UKBB_pop[,c("V",var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)))
    expit_beta-expit_id
  }) #col is SNP #row is sample 
  var_22<-crossprod(u2_theta_coef*UKBB_pop[,var_SNP],u2_theta_coef*UKBB_pop[,var_SNP])
  var_12<-crossprod(UKBB_pop[,-1]*c(expit_beta-UKBB_pop[,1]),u2_theta_coef*UKBB_pop[,var_SNP])
  (1/N_Pop)*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
}

var_theta1_hat_func<-function(
  UKBB_pop,
  theta_UKBB_GPC,
  study_info){
  N_Pop<-nrow(UKBB_pop)
  expit_PC<-(as.matrix(UKBB_pop[,c("V",var_GPC)])%*%theta_UKBB_GPC)
  mat_inside<-(1/N_Pop)*crossprod(c(UKBB_pop[,"Y"]-expit_PC)*UKBB_pop[,c("V",var_GPC)],c(UKBB_pop[,"Y"]-expit_PC)*UKBB_pop[,c("V",var_GPC)])
  mat_outside<-inv_grad_U3_wrt_theta1_func(
    UKBB_pop=UKBB_pop,
    theta_UKBB_GPC = theta_UKBB_GPC)
  mat_outside%*%mat_inside%*%t(mat_outside)
}
var_theta2_hat_vec_func<-function(study_info){
  var_vec<-sapply(1:N_SNP, function(i){
    study_info[[i]]$Covariance
  })
}

cov_U_with_theta_hat_func<-function(
  UKBB_pop,
  beta,
  theta_UKBB_GPC,
  study_info){
  N_Pop<-nrow(UKBB_pop)
  expit_beta<-(UKBB_pop[,-1]%*%beta)
  expit_PC<-(as.matrix(UKBB_pop[,c("V",var_GPC)])%*%theta_UKBB_GPC)
  cov_U1_theta_hat<-(1/N_Pop)*crossprod(UKBB_pop[,-1]*c(expit_beta-UKBB_pop[,1]),UKBB_pop[,c("V",var_GPC)]*c(UKBB_pop[,1]-expit_PC))
  u2_theta_coef<-sapply(1:N_SNP, function(snp_id){
    expit_id<-c((UKBB_pop[,c("V",var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)))
    expit_beta-expit_id
  }) #col is SNP #row is sample 
  cov_U2_theta_hat<-(1/N_Pop)*crossprod(u2_theta_coef*UKBB_pop[,var_SNP],UKBB_pop[,c("V",var_GPC)]*c(UKBB_pop[,1]-expit_PC))
  rbind(cov_U1_theta_hat,cov_U2_theta_hat)
}

final_var_U_beta_theta_hat_func<-function(
  UKBB_pop,
  beta,
  theta_UKBB_GPC,
  study_info,
  len_U1,
  len_theta){
  N_Pop<-nrow(UKBB_pop)
  var_1st_U_beta_theta<-var_U_beta_theta_func(
    UKBB_pop=UKBB_pop,
    theta_UKBB_GPC = theta_UKBB_GPC,
    study_info = study_info,
    beta = beta
  )
  var_theta2_vec<-var_theta2_hat_vec_func(study_info = study_info)
  U_theta_gradient<-rbind(grad_U1_wrt_theta_func(len_U1 = len_U1,len_theta = len_theta),
    grad_U2_wrt_theta_func(UKBB_pop,theta_UKBB_GPC,study_info))
    var_theta1<-var_theta1_hat_func(
      UKBB_pop = UKBB_pop,
      theta_UKBB_GPC = theta_UKBB_GPC,
      study_info = study_info)
    var_grad_times_theta_hat_fromPC<-U_theta_gradient[,1:(1+len_GPC)]%*%var_theta1%*%t(U_theta_gradient[,1:(1+len_GPC)])
    var_grad_times_theta_hat_fromSNP<-U_theta_gradient[,-(1:(1+len_GPC))]%*%(var_theta2_vec*t(U_theta_gradient[,-(1:(1+len_GPC))]))*N_Pop
    var_2nd_grad_times_theta_hat<-var_grad_times_theta_hat_fromPC+var_grad_times_theta_hat_fromSNP## There 
    mat_outside<-inv_grad_U3_wrt_theta1_func(
      UKBB_pop=UKBB_pop,
      theta_UKBB_GPC = theta_UKBB_GPC)
    cov_U<-cov_U_with_theta_hat_func(UKBB_pop = UKBB_pop,beta = beta,theta_UKBB_GPC = theta_UKBB_GPC,study_info = study_info)
    cov_3rd_between_1st_2nd<-cov_U%*%mat_outside%*%t(U_theta_gradient[,1:(1+len_GPC)])
    
    res<-(var_1st_U_beta_theta+var_2nd_grad_times_theta_hat+cov_3rd_between_1st_2nd+t(cov_3rd_between_1st_2nd))
  res
}


