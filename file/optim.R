# Use optim function, Given C, calculate estimated coefficient
useoptim<-function(no_of_studies, 
  study_info, 
  ref_dat, 
  X_rbind, 
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
  D = NULL)
{
  ## Define Cost Function U^T*C*U
  UTCU<-function(coef,C){
    r_first_U <- c()
    r_second_U <- c()
    U <- matrix(NA, nrow(C), nrow(ref_dat))
    
    r = c()
    r_first <- c()
    r_second <- c()
    
    k_j  <- 1
    e1 <- as.vector(1/(1 + exp(-ref_dat %*% coef)))
    e2 <- as.vector((exp(-ref_dat %*% coef) - 1)/(exp(-ref_dat %*% coef) + 1))
    e3 <- as.vector(1/(1 + exp(ref_dat %*% coef)))
    
    #-----------------------------
    #This seems to define some estimating equations
    # I need to write down the estimating equation from linear and logistic regression 
    for(k in 1 : no_of_studies)
    {
      r_first[[k]] <- e1
      col_ind <-  which(colnames(ref_dat) %in% names(study_info[[k]][[1]]) == TRUE)
      r_second[[k]] <-  as.vector(1/(1 + exp(-ref_dat[,col_ind] %*% study_info[[k]][[1]])))
      ncol_study <- length(col_ind)
      k_j <- k_j + ncol_study
      r = c(r, (r_first[[k]] - r_second[[k]]))
    }
    
    study_indices <- seq(1,no_of_studies,1)
    non_missing_covariance_study_indices <- study_indices[-which(study_indices %in% missing_covariance_study_indices)]
    #print(dim(Lambda_ref))
    k_U = 1
    r_first_1 <- as.vector(1/(1 + exp(-ref_dat %*% coef)))
    
    for(k in 1: no_of_studies)
    {
      col_ind <-  which(colnames(ref_dat) %in% names(study_info[[k]][[1]]) == TRUE)
      r_first_U[[k]] <- r_first_1
      r_second_U[[k]] <- as.vector(1/(1 + exp(-ref_dat[,col_ind] %*% study_info[[k]][[1]])))
      U[k_U:(k_U + length(col_ind)-1), ] <- t(ref_dat[ ,col_ind]) %*% diag(r_first_U[[k]]-r_second_U[[k]])
      k_U <- k_U + length(col_ind)
    }
    
    
    Un = rowMeans(U)
    hatQ = t(Un)%*%C%*%Un
    
    if (lambda.gam == 0){
      res<-hatQ[1]*nrow(ref_dat)*nrow(ref_dat)
    }else{
      res<-hatQ[1]*nrow(ref_dat)*nrow(ref_dat) + lambda.gam*coef%*%D%*%coef
    }
    res
  }
  
  ## Define Gradient of Cost Function 
  grad<-function(coef,C)
  {
    r = c()
    r_first <- c()
    r_second <- c()
    Dn_1 <- matrix(NA, ncol(ref_dat), nrow(C))
    Dn_2 <- c()
    
    r_first_U <- c()
    r_second_U <- c()
    U <- matrix(NA, nrow(C), nrow(ref_dat))
    k_j <- 1
    W_temp_1 <- as.vector((1/(1 + exp(-ref_dat %*% coef)))*(1/(1 + exp(ref_dat %*% coef))))
    W <- diag(W_temp_1)
    
    #W <- diag(rep(1/W_temp_1, no_of_studies))
    #print(is.nan(W))
    e1 <- as.vector(1/(1 + exp(-ref_dat %*% coef)))
    e2 <- as.vector((exp(-ref_dat %*% coef) - 1)/(exp(-ref_dat %*% coef) + 1))
    e3 <- as.vector(1/(1 + exp(ref_dat %*% coef)))
    
    
    #-----------------------------
    #This seems to define some estimating equations
    # I need to write down the estimating equation from linear and logistic regression 
    for(k in 1 : no_of_studies)
    {
      r_first[[k]] <- e1
      col_ind <-  which(colnames(ref_dat) %in% names(study_info[[k]][[1]]) == TRUE)
      r_second[[k]] <-  as.vector(1/(1 + exp(-ref_dat[,col_ind] %*% study_info[[k]][[1]])))
      Dn_2 <- c(Dn_2, as.vector(t(ref_dat[,col_ind]) %*% (r_first[[k]] - r_second[[k]])))
      ncol_study <- length(col_ind)
      Dn_1[,k_j: (k_j + ncol_study -1)] <- t(ref_dat) %*% W %*% ref_dat[, col_ind]
      k_j <- k_j + ncol_study
      r = c(r, (r_first[[k]] - r_second[[k]]))
    }
    
    Dn <- Dn_1 %*% C %*% Dn_2
    if(lambda.gam == 0){
      res_grad<-Dn
    }else{
      res_grad<-Dn+lambda.gam*(D%*%coef)
    }
    res_grad
  } 
  
  ## Define Hessian Matrix of Cost Function
  Hessian<-function(coef,C){
    X_abdiag = Reduce(magic::adiag,X_bdiag_list)
    r_first_U <- c()
    r_second_U <- c()
    U <- matrix(NA, nrow(C), nrow(ref_dat))
    Gamma_hat <- matrix(NA, nrow(C), ncol(ref_dat))
    
    r = c()
    r_first <- c()
    r_second <- c()
    W_star_first_list <- list()
    V <- c()
    
    k_j  <- 1
    
    W_temp_1 <- as.vector((1/(1 + exp(-ref_dat %*% coef)))*(1/(1 + exp(ref_dat %*% coef))))
    W <- diag(W_temp_1)
    
    e1 <- as.vector(1/(1 + exp(-ref_dat %*% coef)))
    e2 <- as.vector((exp(-ref_dat %*% coef) - 1)/(exp(-ref_dat %*% coef) + 1))
    e3 <- as.vector(1/(1 + exp(ref_dat %*% coef)))
    l <- e1*e2*e3
    nan_indices <- which(l %in% NaN == TRUE)
    l[nan_indices] <- 0
    L <- diag(l)
    
    for(k in 1 : no_of_studies)
    {
      r_first[[k]] <- e1
      col_ind <-  which(colnames(ref_dat) %in% names(study_info[[k]][[1]]) == TRUE)
      r_second[[k]] <-  as.vector(1/(1 + exp(-ref_dat[,col_ind] %*% study_info[[k]][[1]])))
      ncol_study <- length(col_ind)
      W_star_first_list[[k]] <- W %*% ref_dat[,col_ind]
      k_j <- k_j + ncol_study
      r = c(r, (r_first[[k]] - r_second[[k]]))
    }
    
    
    L = diag(rep(l, no_of_studies))
    
    V = as.vector(t(r) %*% X_abdiag %*% C %*% t(X_abdiag) %*% L)
    #print(class(V))
    W_star_second <- diag(V)
    
    W_star_first <- Reduce(magic::adiag,W_star_first_list) %*% C %*% t(Reduce(magic::adiag,W_star_first_list))
    W_star <- W_star_first + W_star_second
    
    t(X_rbind) %*% W_star %*% X_rbind
  }
  
  ## Calculate new C based on new estimated coef
  C_beta<-function(coef,C){
    r_first_U <- c()
    r_second_U <- c()
    U <- matrix(NA, nrow(C), nrow(ref_dat))
    
    r = c()
    r_first <- c()
    r_second <- c()
    
    k_j  <- 1
    
    W_temp_1 <- as.vector((1/(1 + exp(-ref_dat %*% coef)))*(1/(1 + exp(ref_dat %*% coef))))
    W <- diag(W_temp_1)
    
    #W <- diag(rep(1/W_temp_1, no_of_studies))
    #print(is.nan(W))
    e1 <- as.vector(1/(1 + exp(-ref_dat %*% coef)))
    e2 <- as.vector((exp(-ref_dat %*% coef) - 1)/(exp(-ref_dat %*% coef) + 1))
    e3 <- as.vector(1/(1 + exp(ref_dat %*% coef)))
    
    #-----------------------------
    #This seems to define some estimating equations
    # I need to write down the estimating equation from linear and logistic regression 
    for(k in 1 : no_of_studies)
    {
      r_first[[k]] <- e1
      col_ind <-  which(colnames(ref_dat) %in% names(study_info[[k]][[1]]) == TRUE)
      r_second[[k]] <-  as.vector(1/(1 + exp(-ref_dat[,col_ind] %*% study_info[[k]][[1]])))
      ncol_study <- length(col_ind)
      k_j <- k_j + ncol_study
      r = c(r, (r_first[[k]] - r_second[[k]]))
    }
    
    study_indices <- seq(1,no_of_studies,1)
    non_missing_covariance_study_indices <- study_indices[-which(study_indices %in% missing_covariance_study_indices)]
    #print(non_missing_covariance_study_indices)
    #print(missing_covariance_study_indices)
    
    #Define lambda_ref here
    lambda_ref <- list()
    if(length(missing_covariance_study_indices) > 0)
    {
      for(k in missing_covariance_study_indices)
      {
        col_ind <-  which(colnames(ref_dat) %in% names(study_info[[k]][[1]]) == TRUE)
        #Define W_lambda_ref_logistic
        w_lambda_ref_logistic_vec <- (((1/(1 + exp(ref_dat[,col_ind] %*% study_info[[k]][[1]])))^2)*(1/(1 + exp(-ref_dat %*% coef)))) + (((1/(1 + exp(-ref_dat[,col_ind] %*% study_info[[k]][[1]])))^2)*(1/(1 + exp(ref_dat %*% coef))))
        #print(w_lambda_ref_logistic_vec )
        W_lambda_ref_logistic <- diag(as.vector(w_lambda_ref_logistic_vec))
        #print(class(W_lambda_ref_logistic))
        lambda_ref[[k]] <- (t(ref_dat[,col_ind]) %*% W_lambda_ref_logistic %*% ref_dat[,col_ind])/(study_info[[k]][[3]])
        #print(dim(lambda_ref[[k]]))
      }
      
      for(k in non_missing_covariance_study_indices)
      {
        col_ind <-  which(colnames(ref_dat) %in% names(study_info[[k]][[1]]) == TRUE)
        #Define W_k here
        temp_weight_logistic <- as.vector((1/(1 + exp(-ref_dat[, col_ind] %*% study_info[[k]][[1]])))*(1/(1 + exp(ref_dat[,col_ind] %*% study_info[[k]][[1]]))))
        #temp_weight_logistic <- (exp(ref_dat[,col_ind] %*% study_info[[k]][[1]]))/((1 + exp(ref_dat[,col_ind] %*% study_info[[k]][[1]]))^2)
        W_k <- t(ref_dat[,col_ind]) %*% diag(temp_weight_logistic) %*% ref_dat[,col_ind]
        lambda_ref[[k]] <- (W_k %*% study_info[[k]][[2]] %*% t(W_k))/(nrow(ref_dat))
        #print(dim(lambda_ref[[k]]))
      }
      
    }
    ## When all the estimates for var-cov are provided
    if(length(missing_covariance_study_indices) == 0)
    {
      #print("all the estimates for var-cov are provided")
      for(k in 1 : no_of_studies)
      {
        col_ind <-  which(colnames(ref_dat) %in% names(study_info[[k]][[1]]) == TRUE)
        #Define W_k here
        temp_weight_logistic <- as.vector((1/(1 + exp(-ref_dat[, col_ind] %*% study_info[[k]][[1]])))*(1/(1 + exp(ref_dat[,col_ind] %*% study_info[[k]][[1]]))))
        #temp_weight_logistic <- (exp(ref_dat[,col_ind] %*% study_info[[k]][[1]]))/((1 + exp(ref_dat[,col_ind] %*% study_info[[k]][[1]]))^2)
        W_k <- t(ref_dat[,col_ind]) %*% diag(temp_weight_logistic) %*% ref_dat[,col_ind]
        lambda_ref[[k]] <- (W_k %*% study_info[[k]][[2]] %*% t(W_k))/(nrow(ref_dat))
        #print(dim(lambda_ref[[k]]))
      }
      
    }
    
    Lambda_ref <- Reduce(magic::adiag,lambda_ref)
    #print(dim(Lambda_ref))
    k_U = 1
    r_first_1 <- as.vector(1/(1 + exp(-ref_dat %*% coef)))
    
    for(k in 1: no_of_studies)
    {
      col_ind <-  which(colnames(ref_dat) %in% names(study_info[[k]][[1]]) == TRUE)
      r_first_U[[k]] <- r_first_1
      r_second_U[[k]] <- as.vector(1/(1 + exp(-ref_dat[,col_ind] %*% study_info[[k]][[1]])))
      U[k_U:(k_U + length(col_ind)-1), ] <- t(ref_dat[ ,col_ind]) %*% diag(r_first_U[[k]]-r_second_U[[k]])
      k_U <- k_U + length(col_ind)
    }
    
    # Defining delta_hat here...
    Delta_hat <- (U %*% t(U))/(nrow(ref_dat))
    
    # Defining optimal C here...
    C_beta <- solve(Lambda_ref + Delta_hat, tol = 1e-60)
    
    C_beta
  }
  
  ## Use optim to calculate new estimated coef
  X_abdiag = Reduce(magic::adiag,X_bdiag_list)
  C_iter<-C
  estimated_coef <- as.vector(initial_val)
  #estimated_coef<-c(1.1605597 ,-2.0090653,  3.6657561 ,-1.8231836, -0.6630903 , 0.7562909)
  res<-optim(par=estimated_coef,
    fn=UTCU,
    C=C_iter,
    gr = grad,
    method = "BFGS",hessian = T)
  
  estimated_coef<-res$par
  
  names(estimated_coef)<-colnames(ref_dat)
  C_iter<-C_beta(estimated_coef,C_iter)
  no_of_iter<-c(res$counts[2])
  names(no_of_iter)<-NULL
  cost_val<-res$value
  if (res$convergence == 0){
    convergence <- TRUE
  }else{
    convergence <- FALSE
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

