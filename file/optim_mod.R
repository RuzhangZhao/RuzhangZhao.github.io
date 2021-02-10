#coef<-initial_val
#C<-C_init
sigmoid<-function(vec){
  index_nonneg<-which(vec >= 0)
  index_neg<-which(vec < 0)
  res<-rep(NA,length(vec))
  res[index_nonneg]<-1/(1+exp(-vec[index_nonneg]))
  res[index_neg]<-exp(vec[index_neg])/(1+exp(vec[index_neg]))
  res
}
# Use optim function, Given C, calculate estimated coefficient
useoptim<-function(no_of_studies, 
  study_info, 
  ref_dat, 
  col_inds,
  X_bdiag_list, 
  C, 
  initial_val, 
  threshold, 
  model,
  missing_covariance_study_indices, 
  variable_intercepts, 
  return_C = FALSE,
  return_Hessian = FALSE,
  lambda.gam = 0,
  maxit = 100,
  D = NULL)
{
  ## Define Cost Function U^T*C*U
  UTCU<-function(coef,C){
    e1 <- sigmoid(ref_dat %*% coef)
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
    Un <-Rfast::rowmeans(U)
    hatQ <- t(Un)%*%C%*%Un
    
    if (lambda.gam == 0){
      res<-hatQ[1]*nrow(ref_dat)
    }else{
      res<-hatQ[1]*nrow(ref_dat) + lambda.gam*(coef%*%D%*%coef)[1]
    }
    res
  }
  
  
  
  ## Define Gradient of Cost Function 
  grad<-function(coef,C)
  {
    e1 <- sigmoid(ref_dat%*%coef)
    W_temp_1 <- c(e1*sigmoid(-ref_dat%*%coef))
    W_temp_1_ref_dat<-crossprod(ref_dat,(W_temp_1 * ref_dat))
    Dn_2<- lapply(1:no_of_studies, function(k){
      col_ind<-col_inds[[k]]
      if(length(study_info[[k]][[1]]) == 1){
        r_second <- sigmoid(ref_dat[,col_ind] * study_info[[k]][[1]])
      }else{
        r_second <- sigmoid(ref_dat[,col_ind] %*% study_info[[k]][[1]])
      }
      c((e1 - r_second)%*%ref_dat[,col_ind])
    })
    Dn_2<-unlist(Dn_2)
    Dn_1 <- lapply(1:no_of_studies, function(k){
      W_temp_1_ref_dat[,col_inds[[k]]]
    })
    Dn_1<-do.call(cbind,Dn_1)
    Dn <- Dn_1 %*% (C %*% Dn_2)
    if(lambda.gam == 0){
      res_grad<-Dn/nrow(ref_dat)
    }else{
      res_grad<-Dn/nrow(ref_dat)+lambda.gam*(D%*%coef)
    }
    res_grad
  } 
  
  
  ## Define Hessian Matrix of Cost Function
  Hessian<-function(coef,C){
    ref_dat_size<-nrow(ref_dat)
    e1 <- sigmoid(ref_dat %*% coef)
    e3 <- sigmoid(-ref_dat %*% coef)
    l <- e1*(1-2*e1)*e3
    nan_indices <- which(l %in% NaN == TRUE)
    l[nan_indices] <- 0
    W_temp_1 <- e1*e3
    
    X_abdiag2<-lapply(1:no_of_studies, function(k){
      col_ind<-col_inds[[k]]
      if(length(study_info[[k]][[1]]) == 1){
        r_list_k<- c(e1 - sigmoid(ref_dat[,col_ind] * study_info[[k]][[1]]))
      }else{
        r_list_k<- c(e1 - sigmoid(ref_dat[,col_ind] %*% study_info[[k]][[1]]))
      }
      c(r_list_k%*%X_bdiag_list[[k]])
    })
    X_abdiag2<-unlist(X_abdiag2)
    X_abdiag3<-c(X_abdiag2%*%C)
    col_ind_culcount<-c()
    initial_count<-1
    for(i in 1:no_of_studies){
      ncol_study<-length(col_inds[[i]])
      col_ind_culcount<-c(col_ind_culcount,initial_count)
      initial_count<-initial_count+ncol_study
    }
    
    V<-lapply(1:no_of_studies, function(k){
      k_j<-col_ind_culcount[k]
      ncol_study<-length(col_inds[[k]])
      V_tmp<-c(l*X_bdiag_list[[k]]%*%X_abdiag3[k_j: (k_j + ncol_study -1)])
      k_j<-k_j+ncol_study
      V_tmp
    })
    V<-Reduce('+',V)
    HH<-crossprod(ref_dat,V*ref_dat)
    col_inds_collect<-unlist(col_inds)
    W_star_first_list <- W_temp_1*ref_dat[,col_inds_collect]
    W_star_first_list <-crossprod(W_star_first_list,ref_dat)
    HH<-HH+crossprod(W_star_first_list,C%*%(W_star_first_list))
    HH
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
    print("inv_C")
    print(det(inv_C))
    print(det(inv_C+(10^(-100/dim(inv_C)[1]))*diag(1,dim(inv_C)[1])))
    if (abs(det(inv_C))>1e-20){
      print("original")
      C_beta <- solve(inv_C, tol = 1e-60)
    }else{
      C_beta <- solve(inv_C+(10^(-60/dim(inv_C)[1]))*diag(1,dim(inv_C)[1]), tol = 1e-60)
    }
    print("C_beta")
    print(det(C_beta))
    C_beta
  }
  
  ## Use optim to calculate new estimated coef
  #X_abdiag = Reduce(magic::adiag,X_bdiag_list)
  C_iter<-C
  print(det(C))
  estimated_coef <- c(initial_val)
  #estimated_coef<-c(1.1605597 ,-2.0090653,  3.6657561 ,-1.8231836, -0.6630903 , 0.7562909)
  print("Once")
  res<-optim(par=estimated_coef,
    fn=UTCU,
    C=C_iter,
    gr = grad,
    method = "BFGS")
  estimated_coef<-res$par
  names(estimated_coef)<-colnames(ref_dat)
  C_iter<-C_beta(estimated_coef,C_iter)
  no_of_iter<-c(res$counts[2])
  names(no_of_iter)<-NULL
  cost_val<-res$value
  if (res$convergence == 0){
    convergence <- TRUE
    print("converge")
  }else{
    convergence <- FALSE
    print("Not converge")
  }
  if (return_C & return_Hessian){
    H_<-Hessian(estimated_coef,C_iter)
    newList<-list("estimated_coef" = estimated_coef,
      "cost_val"=cost_val,"no_of_iter"=no_of_iter,
      "convergence"=convergence,"C"=C_iter,"Hessian" = H_)  
  }else if (return_C & (!return_Hessian)){
    newList<-list("estimated_coef" = estimated_coef,
      "cost_val"=cost_val,"no_of_iter"=no_of_iter,
      "convergence"=convergence,"C"=C_iter)
  }else if ((!return_C) & return_Hessian){
    H_<-Hessian(estimated_coef,C_iter)
    newList<-list("estimated_coef" = estimated_coef,
      "cost_val"=cost_val,"no_of_iter"=no_of_iter,
      "convergence"=convergence,"Hessian" = H_)  
  }else{
    newList<-list("estimated_coef" = estimated_coef,
      "cost_val"=cost_val,"no_of_iter"=no_of_iter,
      "convergence"=convergence) 
  }
    newList
}



###############--------------Lambda Tune---------------##############
d_det<-function(H_lD,D){
  tr<-0
  for(pos in which(diag(D)!=0)){
    tr <-tr + det(H_lD[-pos,-pos])
  }
  tr
}

grad_lambda<-function(lambda,coeff,D,H){
  H_lD<-H + lambda*D
  
  -sum(D)/lambda + coeff%*%D%*%coeff + d_det(H_lD,D)/det(H_lD)
}

hessian_lambda<-function(lambda,D,H){
  H_lD<-H + lambda*D
  tr<-0
  
  for(pos in which(diag(D)!=0)){
    tr <-tr + d_det(H_lD = H_lD[-pos,-pos],D = D[-pos,-pos])
  }
  sum(D)/lambda^2 -(d_det(H_lD,D)/det(H_lD))^2+ tr/det(H_lD)
}

