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
    u2_id<-UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]%*%c(study_info[[snp_id]]$Coeff)
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
    
    u2id_gammaid_gradient<-dexpit_id%*%c(UKBB_pop[,paste0("SNP",snp_id)])
    c(u2id_gammaid_gradient)
  })
  #u2_theta is (1+len_GPC+1)*N_SNP
  N_snp_id<-nrow(u2_theta)
  # Gradient is N_equation*N_variables
  -(1/N_Pop)*t(rbind(diag(u2_theta)))
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
    expit_id<-c((UKBB_pop[,c(paste0("SNP",snp_id))]*c(study_info[[snp_id]]$Coeff)))
    expit_beta-expit_id
  }) #col is SNP #row is sample 
  var_22<-crossprod(u2_theta_coef*UKBB_pop[,var_SNP],u2_theta_coef*UKBB_pop[,var_SNP])
  var_12<-crossprod(UKBB_pop[,-1]*c(expit_beta-UKBB_pop[,1]),u2_theta_coef*UKBB_pop[,var_SNP])
  (1/N_Pop)*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
}

var_theta2_hat_vec_func<-function(study_info){
  var_vec<-sapply(1:N_SNP, function(i){
    study_info[[i]]$Covariance
  })
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
 
  #var_grad_times_theta_hat_fromPC<-U_theta_gradient[,1:(1+len_GPC)]%*%var_theta1%*%t(U_theta_gradient[,1:(1+len_GPC)])
  var_grad_times_theta_hat_fromSNP<-U_theta_gradient[,-(1:(1+len_GPC))]%*%(var_theta2_vec*t(U_theta_gradient[,-(1:(1+len_GPC))]))*N_Pop
  var_2nd_grad_times_theta_hat<-var_grad_times_theta_hat_fromSNP## There 
 
  res<-(var_1st_U_beta_theta+var_2nd_grad_times_theta_hat)
  res
}


