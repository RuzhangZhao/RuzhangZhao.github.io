library(reticulate)
library(glue)
loadtorchoptim<-function(){
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
    def forward(self,ref_dat,no_of_studies,col_inds,study_info1,C_):
        ref_size = ref_dat.shape[0]
        e1 = self.sig(self.fc(ref_dat))
        e1 = e1.squeeze()
        k=1
        col_ind = col_inds[str(k)]
        col_ind = [ind-1 for ind in col_ind]
        tmp_ref_dat = ref_dat[:,col_ind]
        r_second_U = torch.matmul(tmp_ref_dat,study_info1[str(k)])
        r_second_U = self.sig(r_second_U)
        r_second_U = e1-r_second_U
        r_second_U = r_second_U.reshape(r_second_U.shape[0],1)
        U = r_second_U*tmp_ref_dat
        for k in range(2,no_of_studies+1):
            col_ind = col_inds[str(k)]
            col_ind = [ind-1 for ind in col_ind]
            tmp_ref_dat = ref_dat[:,col_ind]
            r_second_U = torch.matmul(tmp_ref_dat,study_info1[str(k)])
            r_second_U = self.sig(r_second_U)
            r_second_U = e1-r_second_U
            r_second_U = r_second_U.reshape(r_second_U.shape[0],1)
            tmp_U = r_second_U*tmp_ref_dat
            U = torch.cat([U,tmp_U],dim=1)
        Un = torch.mean(U,dim=0)
        Un = Un.squeeze()
        hatQ = torch.matmul(Un,torch.matmul(C_,Un))
        return hatQ*ref_size

def torchoptim(ref_dat,no_of_studies,col_inds,study_info1,C_,initial_val,lam,D=None,penalty='L2',lr=0.01,lr_step_size=100,EPOCH=1000,eps_early_stop = 10**(-4),verbose=True,verbose_step=10):
    ref_dat = torch.FloatTensor(np.array(ref_dat))
    net = GMMNet(ref_dat.shape[1])
    initial_val_tensor = torch.tensor(np.array(initial_val),dtype=torch.float)
    initial_val_tensor = torch.nn.Parameter(initial_val_tensor.reshape(1,ref_dat.shape[1]))
    net.fc.weight = initial_val_tensor
    net.fc.bias = torch.nn.Parameter(torch.tensor([0.]))
    net.fc.bias.requires_grad_(False)
    for k in range(1,no_of_studies+1):
        study_info1[str(k)] = torch.FloatTensor(np.array(study_info1[str(k)]))
    C_ = torch.FloatTensor(np.array(C_))
    def null_loss(x):
        return x

    loss_metric = null_loss
    ## Initialize Optimizer and Learning Rate Scheduler
    learning_rate = lr
    optimizer = torch.optim.Adam(net.parameters(),lr=learning_rate)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=lr_step_size, gamma=0.1)

    loss_collect = [100]
    count_early = 0
    for i in range(EPOCH):
        outputs = net(ref_dat,no_of_studies,col_inds,study_info1,C_)
        loss = loss_metric(outputs)  
        penalty_val = torch.tensor(0.)
        if lam != 0:
            pen_loc = np.where(np.diag(D)==1)[0]
            if penalty == 'L2':
                penalty_val += (torch.norm(net.fc.weight[0,pen_loc]))**2*lam
            elif penalty == 'L1':
                penalty_val += torch.norm(net.fc.weight[0,pen_loc],1)*lam*2
            else:
                D = torch.FloatTensor(np.array(D))
                penalty_val += lam*(torch.matmul(torch.matmul(net.fc.weight,D),net.fc.weight.T).squeeze())
        loss += penalty_val
        if i == 0:
            print('Initial Loss: '+str(loss.item()), flush = True)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        scheduler.step()       
        if verbose:
            if (i+1)%verbose_step == 0 and i > 0:
                print('Epoch '+str(i+1)+' Loss: '+str(loss.item()), flush = True)
        loss_collect.append(loss.item())
        ## early stopping criteria
        if np.abs(loss_collect[-1] - loss_collect[-2]) < eps_early_stop:
            count_early += 1
        if count_early >= 5:
            break
    net_fc_weight = [i.item() for i in net.fc.weight.squeeze()]
    if penalty == 'L1':
        eps = 10**(-6)
        for i in range(len(net_fc_weight)):
            if np.abs(net_fc_weight[i]) < eps:
                net_fc_weight[i] = 0.
        net_fc_weight = [i for i in net_fc_weight]
        final_tensor = torch.tensor(np.array(net_fc_weight),dtype=torch.float)
        final_tensor = torch.nn.Parameter(final_tensor.reshape(1,ref_dat.shape[1]))
        net.fc.weight = final_tensor
    outputs = net(ref_dat,no_of_studies,col_inds,study_info1,C_)
    penalty_val = torch.tensor(0.)
    if lam != 0:
        pen_loc = np.where(np.diag(D)==1)[0]
        if penalty == 'L2':
            penalty_val += (torch.norm(net.fc.weight[0,pen_loc]))**2*lam
        elif penalty == 'L1':
            penalty_val += torch.norm(net.fc.weight[0,pen_loc],1)*lam*2
        else:
            D = torch.FloatTensor(np.array(D))
            penalty_val += lam*(torch.matmul(torch.matmul(net.fc.weight,D),net.fc.weight.T).squeeze())
    print('Final Loss: '+str((outputs+penalty_val).item()), flush = True)
                 
    return net_fc_weight,np.array(loss_collect[1:len(loss_collect)]),outputs.item()
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
}
fastoptim<-function(study_info, 
  ref_dat, 
  col_inds,
  C, 
  initial_val, 
  model,
  missing_covariance_study_indices, 
  learning_rate = 0.01,
  lr_step_size = 100,
  EPOCH = 10000,
  lam = 0,
  D = NULL,
  penalty='L2',
  verbose = TRUE,
  verbose_step = 10,
  return_C = FALSE,
  return_asy_var = TRUE)
{
  EPOCH<-as.integer(EPOCH)
  no_of_studies<-length(study_info)
  sigmoid<-function(vec){
    index_nonneg<-which(vec >= 0)
    index_neg<-which(vec < 0)
    res<-rep(NA,length(vec))
    res[index_nonneg]<-1/(1+exp(-vec[index_nonneg]))
    res[index_neg]<-exp(vec[index_neg])/(1+exp(vec[index_neg]))
    res
  }
  
  ## Calculate new C based on new estimated coef
  C_beta<-function(coef,C){
    
    e1 <- sigmoid(ref_dat %*% coef)
    e3 <- sigmoid(-ref_dat %*% coef)
    if(length(missing_covariance_study_indices) > 0)
    {
      lambda_ref<-lapply(1:no_of_studies, function(k){
        if(k %in% missing_covariance_study_indices){
          col_ind <-  col_inds[[k]]
          if(length(study_info[[k]][[1]]) == 1){
            tmp<-ref_dat[, col_ind] * study_info[[k]][[1]]
          }else{
            tmp<-ref_dat[, col_ind] %*% study_info[[k]][[1]]
          }
          e1_tmp<-sigmoid(tmp)
          e3_tmp<-sigmoid(-tmp)
          w_lambda_ref_logistic_vec <- (e3_tmp^2)*e1 + (e1_tmp^2)*e3
          W_lambda_ref_logistic <- diag(c(w_lambda_ref_logistic_vec))
          lambda_ref_k<-crossprod(ref_dat[,col_ind],(w_lambda_ref_logistic_vec*ref_dat[,col_ind]))/study_info[[k]][[3]]
        }else{
          col_ind <-  col_inds[[k]]
          if(length(study_info[[k]][[1]]) == 1){
            tmp<-ref_dat[, col_ind] * study_info[[k]][[1]]
          }else{
            tmp<-ref_dat[, col_ind] %*% study_info[[k]][[1]]
          }
          temp_weight_logistic <- sigmoid(tmp)*sigmoid(-tmp)
          W_k <- crossprod(ref_dat[,col_ind],temp_weight_logistic*ref_dat[,col_ind])
          lambda_ref_k <- (W_k %*% study_info[[k]][[2]] %*% t(W_k))/(nrow(ref_dat))
        }
        lambda_ref_k
      })
      for(k in missing_covariance_study_indices)
      {
        col_ind <-  col_inds[[k]]
        if(length(study_info[[k]][[1]]) == 1){
          tmp<-ref_dat[, col_ind] * study_info[[k]][[1]]
        }else{
          tmp<-ref_dat[, col_ind] %*% study_info[[k]][[1]]
        }
        e1_tmp<-sigmoid(tmp)
        e3_tmp<-sigmoid(-tmp)
        w_lambda_ref_logistic_vec <- (e3_tmp^2)*e1 + (e1_tmp^2)*e3
        W_lambda_ref_logistic <- diag(c(w_lambda_ref_logistic_vec))
        lambda_ref[[k]] <- crossprod(ref_dat[,col_ind],(w_lambda_ref_logistic_vec*ref_dat[,col_ind]))/study_info[[k]][[3]]
      }
      study_indices <- seq(1,no_of_studies,1)
      non_missing_covariance_study_indices <- study_indices[-which(study_indices %in% missing_covariance_study_indices)]
      for(k in non_missing_covariance_study_indices)
      {
        col_ind <-  col_inds[[k]]
        #Define W_k here
        if(length(study_info[[k]][[1]]) == 1){
          tmp<-ref_dat[, col_ind] * study_info[[k]][[1]]
        }else{
          tmp<-ref_dat[, col_ind] %*% study_info[[k]][[1]]
        }
        temp_weight_logistic <- sigmoid(tmp)*sigmoid(-tmp)
        W_k <- crossprod(ref_dat[,col_ind],temp_weight_logistic*ref_dat[,col_ind])
        lambda_ref[[k]] <- (W_k %*% study_info[[k]][[2]] %*% t(W_k))/(nrow(ref_dat))
      }
      
    }
    ## When all the estimates for var-cov are provided
    if(length(missing_covariance_study_indices) == 0)
    {
      #print("all the estimates for var-cov are provided")
      lambda_ref<-lapply(1:no_of_studies, function(k){
        col_ind <-  col_inds[[k]]
        #Define W_k here
        if(length(study_info[[k]][[1]]) == 1){
          tmp<-ref_dat[, col_ind] * study_info[[k]][[1]]
        }else{
          tmp<-ref_dat[, col_ind] %*% study_info[[k]][[1]]
        }
        temp_weight_logistic <- sigmoid(tmp)*sigmoid(-tmp)
        W_k <- crossprod(ref_dat[,col_ind],temp_weight_logistic*ref_dat[,col_ind])
        (W_k %*% study_info[[k]][[2]] %*% t(W_k))/(nrow(ref_dat))
      })
      
    }
    
    Lambda_ref <- Reduce(magic::adiag,lambda_ref)
    #print(dim(Lambda_ref))
    U<-lapply(1:no_of_studies, function(k){
      col_ind <-  col_inds[[k]]
      if(length(study_info[[k]][[1]]) == 1){
        r_second_U <- sigmoid(ref_dat[,col_ind] * study_info[[k]][[1]])
      }else{
        r_second_U <- sigmoid(ref_dat[,col_ind] %*% study_info[[k]][[1]]) 
      }
      t((e1-r_second_U)*ref_dat[ ,col_ind])
    })
    U<-do.call(rbind,U)
    # Defining delta_hat here...
    Delta_hat <- (U %*% t(U))/(nrow(ref_dat))
    
    # Defining optimal C here...
    inv_C<-Lambda_ref + Delta_hat
    #print("inv_C")
    #print(det(inv_C))
    #print(det(inv_C+(mean(diag(inv_C))*0.01)*diag(1,dim(inv_C)[1])))
    #print(mean(diag(inv_C))*0.01)
    if (abs(det(inv_C))>1e-20){
      #print("original")
      C_beta <- solve(inv_C, tol = 1e-60)
    }else{
      C_beta <- solve(inv_C+(mean(diag(inv_C))*0.01)*diag(1,dim(inv_C)[1]), tol = 1e-60)
    }
    #print("C_beta")
    #print(det(C_beta))
    C_beta
  }
  
  asy_var_opt<-function(coef,C){
    
    e1 <- sigmoid(ref_dat %*% coef)
    e3 <- sigmoid(-ref_dat %*% coef)
    W_Gamma_temp_1 <-e1*e3 
    
    Gamma_hat<-lapply(1:no_of_studies, function(k){
      col_ind <-col_inds[[k]]
      t(W_Gamma_temp_1*ref_dat[ ,col_ind])%*% ref_dat
    })
    Gamma_hat<-do.call(rbind,Gamma_hat)
    Gamma_hat <- Gamma_hat/nrow(ref_dat)
    info <- (crossprod(Gamma_hat, C)%*%Gamma_hat)
    if(det(info) == 0){ 
      asy_var = NULL
    }else{
      asy_var<-solve(info, tol = 1e-60)/nrow(ref_dat)
    }
    asy_var
  }
  
  C_iter<-C
  estimated_coef <- c(initial_val)
  loadtorchoptim()
  study_info1<-list()
  for (i in 1:length(study_info)){
    study_info1[[paste0(i)]]<-c(study_info[[i]][[1]])
  }
  col_inds1<-list()
  for (i in 1:length(col_inds)){
    col_inds1[[paste0(i)]]<-col_inds[[i]]
  }
  res<-torchoptim(
    ref_dat=ref_dat,
    no_of_studies = no_of_studies,
    col_inds = col_inds1,
    study_info1 = study_info1,
    C_ = C_iter,
    initial_val = estimated_coef,
    lam = lam,
    D = D,
    penalty=penalty,
    lr = learning_rate,
    lr_step_size = lr_step_size,
    EPOCH=EPOCH,
    verbose=verbose,
    verbose_step =verbose_step)
  estimated_coef<-res[[1]]
  names(estimated_coef)<-colnames(ref_dat)
  print(paste0("func value before:",res[[2]][1]))
  print(paste0("func value after:",res[[2]][length(res[[2]])]))
  cost_val<-res[[3]]
  newList<-list("estimated_coef" = estimated_coef,
    "cost_val"=cost_val,"no_of_iter"=length(res[[2]])) 
  if(return_C){
    C_new<-C_beta(estimated_coef,C_iter)
    newList$C <- C_new
  }
  if(return_asy_var){
    print("Yes")
    asy_var<-asy_var_opt(estimated_coef,C_iter)
    newList$asy_var<-asy_var
  }
  newList
}

