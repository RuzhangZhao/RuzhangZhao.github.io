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
# compute the variance of ridge regression
var_ridge<-function(
  UKBB_pop,
  beta,
  A_penalty,
  lambda){
  A_penalty_with_lambda<-A_penalty*lambda
  dexpit_beta<-dexpit(UKBB_pop[,-1]%*%as.vector(beta))
  Hessian_<-crossprod(UKBB_pop[,-1]*c(dexpit_beta),UKBB_pop[,-1])
  diag(solve(Hessian_+A_penalty_with_lambda))
}

U_func<-function(
  UKBB_pop,
  beta,
  theta_UKBB_GPC,
  study_info){
  N_Pop<-nrow(UKBB_pop)
  expit_beta<-expit(UKBB_pop[,-1]%*%beta)
  u1<-crossprod((UKBB_pop[,1] -expit_beta),UKBB_pop[,-1])
  u2<-c()
  u2_part1<-c(expit_beta)%*%UKBB_pop[,var_SNP]
  u2_part2<-sapply(1:N_SNP, function(snp_id){
    u2_id<-expit(UKBB_pop[,c("V",var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff))
    c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
  })
  u2<-u2_part1 - u2_part2
  u<-c(u1,u2)*(1/N_Pop)
}
dexpit<-function(x){
  expit(x)*(1-expit(x))
}

grad_U_wrt_beta_func<-function(
  UKBB_pop,
  beta){
  N_Pop<-nrow(UKBB_pop)
  dexpit_beta<-dexpit(UKBB_pop[,-1]%*%beta)
  U1_beta_gradient<-crossprod(UKBB_pop[,-1]*c(dexpit_beta),UKBB_pop[,-1])*(-1)
  U2_beta_gradient<-crossprod(UKBB_pop[,var_SNP]*c(dexpit_beta),UKBB_pop[,-1])
  rbind(U1_beta_gradient,U2_beta_gradient)*(1/N_Pop)
}

U3_func<-function(UKBB_pop,theta_UKBB_GPC){
  expit_PC<-expit(UKBB_pop[,c("V",var_GPC)]%*%theta_UKBB_GPC)
  # UKBB_pop[,1] is Y(phenotype)
  crossprod((UKBB_pop[,"Y"]-expit_PC),UKBB_pop[,c("V",var_GPC)])
}

### Generate coefficients of GPCs. 
inv_grad_U3_wrt_theta1_func<-function(
  UKBB_pop,
  theta_UKBB_GPC){
  N_Pop<-nrow(UKBB_pop)
  dexpit_PC<-dexpit(as.matrix(UKBB_pop[,c("V",var_GPC)])%*%theta_UKBB_GPC)
  mat<--(1/N_Pop)*crossprod(UKBB_pop[,c("V",var_GPC)]*c(dexpit_PC),UKBB_pop[,c("V",var_GPC)])
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
    dexpit_id<-c(UKBB_pop[,paste0("SNP",snp_id)])*c(dexpit(UKBB_pop[,c("V",var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)))
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
  expit_beta<-expit(UKBB_pop[,-1]%*%beta)
  var_11<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-expit_beta),UKBB_pop[,-1]*c(UKBB_pop[,1]-expit_beta))
  u2_theta_coef<-sapply(1:N_SNP, function(snp_id){
    expit_id<-c(expit(UKBB_pop[,c("V",var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)))
    expit_beta-expit_id
  }) #col is SNP #row is sample 
  var_22<-crossprod(u2_theta_coef*UKBB_pop[,var_SNP],u2_theta_coef*UKBB_pop[,var_SNP])
  var_12<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-expit_beta),u2_theta_coef*UKBB_pop[,var_SNP])
  (1/N_Pop)*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
}

var_theta1_hat_func<-function(
  UKBB_pop,
  theta_UKBB_GPC,
  study_info){
  N_Pop<-nrow(UKBB_pop)
  expit_PC<-expit(as.matrix(UKBB_pop[,c("V",var_GPC)])%*%theta_UKBB_GPC)
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
  expit_beta<-expit(UKBB_pop[,-1]%*%beta)
  expit_PC<-expit(as.matrix(UKBB_pop[,c("V",var_GPC)])%*%theta_UKBB_GPC)
  cov_U1_theta_hat<-(1/N_Pop)*crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-expit_beta),UKBB_pop[,c("V",var_GPC)]*c(UKBB_pop[,1]-expit_PC))
  u2_theta_coef<-sapply(1:N_SNP, function(snp_id){
    expit_id<-c(expit(UKBB_pop[,c("V",var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)))
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



utcu_C<-function(u,C){
  c(u%*%C%*%u)
}
utcu<-function(beta,C,lambda,A_penalty){
  u<-U_func(UKBB_pop,beta,theta_UKBB_GPC,study_info)
  # we only need to consider the optimization 
  # problem related with beta 
  utcu_C(u,C)+lambda*c(c(beta)%*%(A_penalty%*%c(beta)))
}

gradient_utcu<-function(beta,C_current){
  N_Pop<-nrow(UKBB_pop)
  u<-U_func(UKBB_pop,beta,theta_UKBB_GPC,study_info)
  R_mat<-as.matrix(UKBB_pop[,-1])
  # R = [1,Z,X] includes "V" as intercept,
  # SNP as "Z" for matrix which we want to penalize
  # genetic PC & risk factors
  Z_mat<-as.matrix(UKBB_pop[,var_SNP])
  # SNP as "Z" for matrix which we want to penalize
  # Since we have the GWAS data for each SNP 
  # we need to have score function for each
  p_vec<-c(expit(R_mat%*%c(beta)))
  # p = expit(R\beta), the estimated log prob
  M_vec<-p_vec*(1-p_vec)
  # M_vec should be diagonal matrix, just use it as 
  # vector here 
  g0<-M_vec*cbind(-R_mat,Z_mat)%*%C_current%*%u
  # g0 is the gradient vector without 
  # considering the desgin matrix part
  g<-c(g0)%*%R_mat/N_Pop
  g
}

Hessian_utcu<-function(beta,C_current){
  N_Pop<-nrow(UKBB_pop)
  u<-U_func(UKBB_pop,beta,theta_UKBB_GPC,study_info)
  R_mat<-as.matrix(UKBB_pop[,-1])
  # R = [1,Z,X] includes "V" as intercept,
  # SNP as "Z" for matrix which we want to penalize
  # genetic PC & risk factors
  Z_mat<-as.matrix(UKBB_pop[,var_SNP])
  # SNP as "Z" for matrix which we want to penalize
  # Since we have the GWAS data for each SNP 
  # we need to have score function for each
  p_vec<-c(expit(R_mat%*%c(beta)))
  # p = expit(R\beta), the estimated log prob
  M_vec<-p_vec*(1-p_vec)
  # M_vec should be diagonal matrix, just use it as 
  # vector here 
  #H0_part1<-M_vec*cbind(-R_mat,Z_mat)%*%C_current%*%t(M_vec*cbind(-R_mat,Z_mat))
  H_part1_left<-t(M_vec*cbind(-R_mat,Z_mat))%*%R_mat
  H_part1<-(t(H_part1_left)%*%H_part1_left)/N_Pop
  #H_part1<-crossprod(R_mat,H0_part1%*%R_mat)
  # The first part of Hessian matrix is U'CU'
  # where both U are taken as derivative
  
  M_long_vec<-M_vec*(1-2*p_vec)
  H0_part2<-c(M_long_vec*cbind(-R_mat,Z_mat)%*%C_current%*%u)
  H_part2<-crossprod(R_mat,H0_part2*R_mat)/N_Pop
  # The second part of Hessian matrix is UCU''
  # where only one U is taken as second derivate
  # We need to discuss whether the second part should be included 
  
  H<-H_part1+H_part2
  # After getting H0, the computation of H is easy
  # Here, we first use H0_part1. 
  #det(H)
  # Notice the objective function is based on UTCU/N
  # where N is sample size
  H
}
tilde_final_var_func<-function(
  UKBB_pop,
  beta,
  theta_UKBB_GPC,
  study_info,
  len_U1,
  len_theta,
  A_penalty,
  lambda){
  N_Pop<-nrow(UKBB_pop)
  A_penalty_with_lambda<-A_penalty*lambda
  tilde_U_beta_gradient<-grad_U_wrt_beta_func(UKBB_pop = UKBB_pop,beta = c(beta))
  
  tilde_var_U_beta_theta<-
    var_U_beta_theta_func(
      UKBB_pop=UKBB_pop,
      beta = beta,
      theta_UKBB_GPC = theta_UKBB_GPC,
      study_info = study_info)
  
  tilde_U_beta_square<-crossprod(tilde_U_beta_gradient,C)%*%tilde_var_U_beta_theta%*%C%*%tilde_U_beta_gradient
  tilde_U_beta_theta<-U_func(
    UKBB_pop = UKBB_pop,
    beta = beta,
    theta_UKBB_GPC =theta_UKBB_GPC,
    study_info = study_info)
  
  tilde_left<-crossprod(tilde_U_beta_gradient,C)%*%tilde_U_beta_theta
  tilde_U_beta_square2<-tilde_left%*%t(tilde_left)
  
  
  ### STop here
  tilde_U_beta_theta_times_Abeta<-crossprod(tilde_U_beta_gradient,C)%*%tilde_U_beta_theta%*%t(A_penalty_with_lambda%*%beta)*sqrt(N_Pop)
  
  tilde_Abeta_square<-(A_penalty_with_lambda%*%beta)%*%t(A_penalty_with_lambda%*%beta)*N_Pop
  
  tilde_var_1st_U_beta_theta<-tilde_U_beta_square - tilde_U_beta_square2
  #tilde_var_1st_U_beta_theta<-tilde_U_beta_square+tilde_U_beta_theta_times_Abeta+t(tilde_U_beta_theta_times_Abeta)+tilde_Abeta_square
  
  ####### This part about var of theta hat is the same 
  
  var_theta2_vec<-var_theta2_hat_vec_func(study_info = study_info)
  U_theta_gradient<-rbind(grad_U1_wrt_theta_func(len_U1 = len_U1,len_theta = len_theta),
    grad_U2_wrt_theta_func(UKBB_pop,theta_UKBB_GPC,study_info))

  var_grad_times_theta_hat_fromSNP<-U_theta_gradient[,-(1:(1+len_GPC))]%*%(var_theta2_vec*t(U_theta_gradient[,-(1:(1+len_GPC))]))*N_Pop

    var_theta1<-var_theta1_hat_func(
      UKBB_pop = UKBB_pop,
      theta_UKBB_GPC = theta_UKBB_GPC,
      study_info = study_info)
    var_theta2_vec<-var_theta2_hat_vec_func(study_info = study_info)
    U_theta_gradient<-rbind(grad_U1_wrt_theta_func(len_U1 = len_U1,len_theta = len_theta),
      grad_U2_wrt_theta_func(UKBB_pop,theta_UKBB_GPC,study_info))
    
    var_grad_times_theta_hat_fromPC<-U_theta_gradient[,1:(1+len_GPC)]%*%var_theta1%*%t(U_theta_gradient[,1:(1+len_GPC)])
    var_grad_times_theta_hat_fromSNP<-U_theta_gradient[,-(1:(1+len_GPC))]%*%(var_theta2_vec*t(U_theta_gradient[,-(1:(1+len_GPC))]))*N_Pop
    
    var_2nd_grad_times_theta_hat<-var_grad_times_theta_hat_fromPC+var_grad_times_theta_hat_fromSNP
    #######
    
    tilde_var_2nd<-crossprod(tilde_U_beta_gradient,C)%*%var_2nd_grad_times_theta_hat%*%C%*%tilde_U_beta_gradient
    
    cov_U<-cov_U_with_theta_hat_func(UKBB_pop = UKBB_pop,beta = beta,theta_UKBB_GPC = theta_UKBB_GPC,study_info = study_info)
    mat_outside<-inv_grad_U3_wrt_theta1_func(
      UKBB_pop=UKBB_pop,
      theta_UKBB_GPC = theta_UKBB_GPC)
    tilde_cov_3rd_between_1st_2nd<-crossprod(tilde_U_beta_gradient,C)%*%cov_U%*%mat_outside%*%t(U_theta_gradient[,1:(1+len_GPC)])%*%C%*%tilde_U_beta_gradient
    
    #print(max(tilde_var_1st_U_beta_theta))
    #print(max(tilde_var_2nd))
    #print(max(tilde_cov_3rd_between_1st_2nd))
    res<-tilde_var_1st_U_beta_theta+tilde_var_2nd+tilde_cov_3rd_between_1st_2nd+t(tilde_cov_3rd_between_1st_2nd)
  
  
  res
}



library(reticulate)
library(glue)

py_run_string(glue(
  "
import torch
import torch.nn as nn
import numpy as np
torch.set_default_tensor_type('torch.FloatTensor')
class GMMNet(nn.Module):
    def __init__(self,n_feature):
        super(GMMNet, self).__init__()
        self.fc = nn.Linear(n_feature,1)
        self.sig = nn.Sigmoid()
    def forward(self,UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_):
        #UKBB_pop = torch.FloatTensor(UKBB_pop)
        #beta = torch.FloatTensor(beta)
        e1 = self.sig(self.fc(UKBB_pop[:,1:UKBB_pop.shape[1]]))
        e1 = e1.squeeze()
        u1 = torch.matmul((UKBB_pop[:,0] - e1),UKBB_pop[:,1:UKBB_pop.shape[1]])
        index_varSNP = [colname_UKBB.index(var_SNP[i]) for i in range(len(var_SNP))]
        u2_part1 = torch.matmul(e1,UKBB_pop[:,index_varSNP])
        if var_GPC is None:
          index_varGPC = []
        else:
          index_varGPC = [colname_UKBB.index(var_GPC[i]) for i in range(len(var_GPC))]
        u2_part2 = torch.empty(len(var_SNP))
        for snp_id in range(1,(len(var_SNP)+1) ):#range(100,101):#

          index_snp_id = [0]+index_varGPC+[colname_UKBB.index('SNP'+str(snp_id))]
          u2_id = torch.matmul(UKBB_pop[:,index_snp_id],torch.FloatTensor(np.append(theta_UKBB_GPC,study_info[snp_id-1]['Coeff'])))
          u2_part2[snp_id-1] =torch.dot(u2_id,UKBB_pop[:,colname_UKBB.index('SNP'+str(snp_id))])
        Un = torch.cat([u1,u2_part1-u2_part2],dim=0)
        Un = Un.squeeze()
        hatQ = torch.matmul(Un,torch.matmul(C_,Un))
        return hatQ/UKBB_pop.shape[0]**2



def torchoptim(UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_,initial_val,D,lam=1,lr=0.01,EPOCH=1000):
    UKBB_pop = torch.FloatTensor(UKBB_pop)
    net = GMMNet(UKBB_pop.shape[1]-1)
    initial_val_tensor = torch.tensor(np.array(initial_val),dtype=torch.float)
    initial_val_tensor = torch.nn.Parameter(initial_val_tensor.reshape(1,(UKBB_pop.shape[1]-1) ))
    net.fc.weight = initial_val_tensor
    net.fc.bias = torch.nn.Parameter(torch.tensor([0.]))
    net.fc.bias.requires_grad_(False)
    #for k in range(1,no_of_studies+1):
    #    study_info1[str(k)] = torch.FloatTensor(np.array(study_info1[str(k)]))
    C_ = torch.FloatTensor(np.array(C_))
    
    def null_loss(x):
        return x

    loss_metric = null_loss
    ## Initialize Optimizer and Learning Rate Scheduler
    learning_rate = lr
    optimizer = torch.optim.Adam(net.parameters(),lr=learning_rate)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=100, gamma=0.1)

    loss_collect = [100]
    count_early = 0
    
    for i in range(EPOCH):
        outputs = net(UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_)
        loss = loss_metric(outputs)       
        penalty_val = torch.tensor(0.)
        if lam != 0:
            pen_loc = np.where(np.diag(D)==1)[0]
            
            penalty_val += (torch.norm(net.fc.weight[0,pen_loc]))**2*lam
        loss += penalty_val
        if i == 0:
            print('Initial Loss: '+str(loss.item()), flush = True)
            initial_loss = loss.item()
            eps_early = initial_loss*10**(-5)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        scheduler.step()
        if (i+1)%100 == 0 or i == 0:
            print('Epoch '+str(i+1)+' Loss: '+str(loss.item()))
        print('Epoch '+str(i+1)+' Loss: '+str(loss.item()))
        loss_collect.append(loss.item())
        ## early stopping criteria

        if np.abs(loss_collect[-1] - loss_collect[-2]) < eps_early:
            count_early += 1
        if count_early >= 10:
            break
    outputs = net(UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_)
    penalty_val = torch.tensor(0.)
    if lam != 0:
        pen_loc = np.where(np.diag(D)==1)[0]
        penalty_val += (torch.norm(net.fc.weight[0,pen_loc]))**2*lam
    print('Final Loss: '+str((outputs+penalty_val).item()), flush = True)
               
    return [i.item() for i in net.fc.weight.squeeze()],np.array(loss_collect[1:len(loss_collect)])
")  
  ,local = FALSE, convert = TRUE)



# copy objects from the main python module into the specified R environment
py_main <- import_main(convert = TRUE)
py_main_dict <- py_get_attr(py_main, "__dict__")
py_name<- "torchoptim"
Encoding(py_name) <- "UTF-8"

py_value <- py_main_dict[[py_name]]
py_envir<-globalenv()
assign(py_name, py_value, envir = py_envir) 





library(reticulate)
library(glue)

py_run_string(glue(
  "
import torch
import torch.nn as nn
import numpy as np
torch.set_default_tensor_type('torch.FloatTensor')
class GMMNet(nn.Module):
    def __init__(self,n_feature):
        super(GMMNet, self).__init__()
        self.fc = nn.Linear(n_feature,1)
        self.sig = nn.Sigmoid()
    def forward(self,UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_):
        #UKBB_pop = torch.FloatTensor(UKBB_pop)
        #beta = torch.FloatTensor(beta)
        e1 = self.sig(self.fc(UKBB_pop[:,1:UKBB_pop.shape[1]]))
        e1 = e1.squeeze()
        u1 = torch.matmul((UKBB_pop[:,0] - e1),UKBB_pop[:,1:UKBB_pop.shape[1]])
        index_varSNP = [colname_UKBB.index(var_SNP[i]) for i in range(len(var_SNP))]
        u2_part1 = torch.matmul(e1,UKBB_pop[:,index_varSNP])
        if var_GPC is None:
          index_varGPC = []
        else:
          index_varGPC = [colname_UKBB.index(var_GPC[i]) for i in range(len(var_GPC))]
        
        u2_part2 = torch.empty(len(var_SNP))
        for snp_id in range(1,(len(var_SNP)+1)  ):#range(100,101):#

          index_snp_id = [0]+index_varGPC+[colname_UKBB.index('SNP'+str(snp_id))]
          u2_id = torch.matmul(UKBB_pop[:,index_snp_id],torch.FloatTensor(np.append(theta_UKBB_GPC,study_info[snp_id-1]['Coeff'])))
          u2_id = self.sig(u2_id)
          u2_part2[snp_id-1] =torch.dot(u2_id,UKBB_pop[:,colname_UKBB.index('SNP'+str(snp_id))])
        Un = torch.cat([u1,u2_part1-u2_part2],dim=0)
        Un = Un.squeeze()
        hatQ = torch.matmul(Un,torch.matmul(C_,Un))
        return hatQ/UKBB_pop.shape[0]**2



def torchoptimLBFGS(UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_,initial_val,D,lam=1,lr=1,EPOCH=10,verbose=False):
    UKBB_pop = torch.FloatTensor(UKBB_pop)
    net = GMMNet(UKBB_pop.shape[1]-1)
    initial_val_tensor = torch.tensor(np.array(initial_val),dtype=torch.float)
    initial_val_tensor = torch.nn.Parameter(initial_val_tensor.reshape(1,(UKBB_pop.shape[1]-1) ))
    net.fc.weight = initial_val_tensor
    net.fc.bias = torch.nn.Parameter(torch.tensor([0.]))
    net.fc.bias.requires_grad_(False)
    #for k in range(1,no_of_studies+1):
    #    study_info1[str(k)] = torch.FloatTensor(np.array(study_info1[str(k)]))
    C_ = torch.FloatTensor(np.array(C_))
    
    def null_loss(x):
        return x

    loss_metric = null_loss
    ## Initialize Optimizer and Learning Rate Scheduler
    learning_rate = lr
    optimizer = torch.optim.LBFGS(net.parameters(),max_iter=30)    
    ## scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=100, gamma=0.1)

    loss_collect = [100]
    count_early = 0
    eps_early = 10**(-10)


    def closure():
        outputs = net(UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_)
        penalty_val = torch.tensor(0.)
        if lam != 0:
            pen_loc = np.where(np.diag(D)==1)[0]
            
            penalty_val += (torch.norm(net.fc.weight[0,pen_loc]))**2*lam
        optimizer.zero_grad() 
        loss = loss_metric(outputs)       
        loss += penalty_val
        loss.backward()
        return loss
    #for i in range(EPOCH):
    #    print('here')
    losses = optimizer.step(closure)
    i=0  
    #for i in range(N_LBFGS_STEPS_VALIDATION):
          
        #optimizer.zero_grad()
        #loss.backward()
        #optimizer.step()
        #scheduler.step()
    if (i+1)%100 == 0 or i == 0:
        if verbose:
            print('Epoch '+str(i+1)+' Loss: '+str(losses.item()))
    if verbose:
        print('Epoch '+str(i+1)+' Loss: '+str(losses.item()))
    loss_collect.append(losses.item())
    ## early stopping criteria

    #if np.abs(loss_collect[-1] - loss_collect[-2]) < eps_early:
    #    count_early += 1
    #if count_early >= 100:
    #    break
    
    if verbose:
        outputs = net(UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_)
        penalty_val = torch.tensor(0.)
        if lam != 0:
            pen_loc = np.where(np.diag(D)==1)[0]
            penalty_val += (torch.norm(net.fc.weight[0,pen_loc]))**2*lam
        print('Final Loss: '+str((outputs+penalty_val).item()), flush = True)
               
    return [i.item() for i in net.fc.weight.squeeze()],np.array(loss_collect[1:len(loss_collect)])
")  
  ,local = FALSE, convert = TRUE)



# copy objects from the main python module into the specified R environment
py_main <- import_main(convert = TRUE)
py_main_dict <- py_get_attr(py_main, "__dict__")
py_name<- "torchoptimLBFGS"
Encoding(py_name) <- "UTF-8"

py_value <- py_main_dict[[py_name]]
py_envir<-globalenv()
assign(py_name, py_value, envir = py_envir) 











library(reticulate)
library(glue)

py_run_string(glue(
  "
import torch
import torch.nn as nn
import numpy as np
torch.set_default_tensor_type('torch.FloatTensor')
class GMMNet(nn.Module):
    def __init__(self,n_feature):
        super(GMMNet, self).__init__()
        self.fc = nn.Linear(n_feature,1)
        self.sig = nn.Sigmoid()
    def forward(self,UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_):
        #UKBB_pop = torch.FloatTensor(UKBB_pop)
        #beta = torch.FloatTensor(beta)
        e1 = self.sig(self.fc(UKBB_pop[:,1:UKBB_pop.shape[1]]))
        e1 = e1.squeeze()
        u1 = torch.matmul((UKBB_pop[:,0] - e1),UKBB_pop[:,1:UKBB_pop.shape[1]])
        index_varSNP = [colname_UKBB.index(var_SNP[i]) for i in range(len(var_SNP))]
        u2_part1 = torch.matmul(e1,UKBB_pop[:,index_varSNP])
        if var_GPC is None:
          index_varGPC = []
        else:
          index_varGPC = [colname_UKBB.index(var_GPC[i]) for i in range(len(var_GPC))]
        
        u2_part2 = torch.empty(len(var_SNP))
        for snp_id in range(1,(len(var_SNP)+1) ):#range(100,101):#

          index_snp_id = [0]+index_varGPC+[colname_UKBB.index('SNP'+str(snp_id))]
          u2_id = torch.matmul(UKBB_pop[:,index_snp_id],torch.FloatTensor(np.append(theta_UKBB_GPC,study_info[snp_id-1]['Coeff'])))
          u2_part2[snp_id-1] =torch.dot(u2_id,UKBB_pop[:,colname_UKBB.index('SNP'+str(snp_id))])
        Un = torch.cat([u1,u2_part1-u2_part2],dim=0)
        Un = Un.squeeze()
        hatQ = torch.matmul(Un,torch.matmul(C_,Un))
        return hatQ/UKBB_pop.shape[0]**2



def torchoptimLBFGSv2(UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_,initial_val,D,lam=1,lr=1,EPOCH=10,verbose=False):
    UKBB_pop = torch.FloatTensor(UKBB_pop)
    net = GMMNet(UKBB_pop.shape[1]-1)
    initial_val_tensor = torch.tensor(np.array(initial_val),dtype=torch.float)
    initial_val_tensor = torch.nn.Parameter(initial_val_tensor.reshape(1,(UKBB_pop.shape[1]-1) ))
    net.fc.weight = initial_val_tensor
    net.fc.bias = torch.nn.Parameter(torch.tensor([0.]))
    net.fc.bias.requires_grad_(False)
    #for k in range(1,no_of_studies+1):
    #    study_info1[str(k)] = torch.FloatTensor(np.array(study_info1[str(k)]))
    C_ = torch.FloatTensor(np.array(C_))
    
    def null_loss(x):
        return x

    loss_metric = null_loss
    ## Initialize Optimizer and Learning Rate Scheduler
    learning_rate = lr
    optimizer = torch.optim.LBFGS(net.parameters())    
    ## scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=100, gamma=0.1)

    loss_collect = [100]
    count_early = 0
    eps_early = 10**(-10)


    def closure():
        outputs = net(UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_)
        penalty_val = torch.tensor(0.)
        if lam != 0:
            penalty_val += (torch.matmul(torch.matmul(net.fc.weight,torch.FloatTensor(D)),net.fc.weight.T)[0][0])**2*lam
        optimizer.zero_grad() 
        loss = loss_metric(outputs)       
        loss += penalty_val
        loss.backward()
        return loss

    losses = optimizer.step(closure)
    
    loss_collect.append(losses.item())
               
    return [i.item() for i in net.fc.weight.squeeze()],np.array(loss_collect[1:len(loss_collect)])
")  
  ,local = FALSE, convert = TRUE)



# copy objects from the main python module into the specified R environment
py_main <- import_main(convert = TRUE)
py_main_dict <- py_get_attr(py_main, "__dict__")
py_name<- "torchoptimLBFGSv2"
Encoding(py_name) <- "UTF-8"

py_value <- py_main_dict[[py_name]]
py_envir<-globalenv()
assign(py_name, py_value, envir = py_envir) 


