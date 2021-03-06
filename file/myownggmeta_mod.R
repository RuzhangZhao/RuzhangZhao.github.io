library(parallel)
ggmeta_mod <- function(study_info, ref_dat, 
  model, 
  estimated_intercept = NULL,
  variable_intercepts=FALSE,
  initial_val=NULL,
  tune_lambda=FALSE,
  lambda.gam=NULL,
  D=NULL,
  lambda.gam.first = 0,
  one_step = TRUE,
  control = list(epsilon = 1e-04, 
    maxit = 500, maxit_lam = 1e3,
    lambda_tune_eps=1e-06))
{
  call_ggmeta <- match.call()
  threshold <- control[[1]]
  maxit <- control[[2]]
  maxit_lam <- control[[3]]
  lambda_tune_eps <- control[[4]]
  no_of_studies <- length(study_info)
  
  error_1 <- 0
  error_2 <- 0
  error_3 <- 0
  error_4 <- 0
  
  temp <- c()
  indicator_missing_covariance_sample_size<-0

  missing_study_sample_size <- c()
  missing_covariance_study_indices <- c()

  
  for(i in 1 : no_of_studies)
  {
    if(is.null(study_info[[i]][[3]]) == T){
      missing_study_sample_size<-c(missing_study_sample_size,i) 
    }
      temp <- union(temp,names(study_info[[i]][[1]]))
    if(is.null(study_info[[i]][[2]]) == T &&
        is.null(study_info[[i]][[3]]) == T){
      indicator_missing_covariance_sample_size<-1    
    }
      if(is.null(study_info[[i]][[2]]) == T){
        missing_covariance_study_indices <- 
          c(missing_covariance_study_indices, i)
      }
      
  }
  
  if(indicator_missing_covariance_sample_size == 1){
    print("Error: All the studies should have either
      an estimate for the var-cov matrix or the sample size.
      At least one of them is missing(NULL) in at least one of
      the studies")
    error_1 <-1
  }
  if(indicator_missing_covariance_sample_size == 0)
  {
    for(i in missing_study_sample_size)
    {
      if(length(study_info[[i]][[1]]) != 
          ncol(study_info[[i]][[2]]))
      {
        print("Error: length of the study parameter(effect sizes)
          vector does not match with the dimension of its variance
          covariance matrix")
        error_2 <- 1
      }
      if(sum(names(study_info[[i]][[1]]) != 
          colnames(study_info[[i]][[2]])) > 0)
      {
        print("Error: names of the variables corresponding to the 
          study specific parameter vector is not same(also not in 
          the same order) as in its variance covariance matrix")
        error_2 <- 1
      }
    }
  }
  
  if(ncol(ref_dat) < length(temp))
  {
    print("number of covariates in the reference data does not 
      match with that of the full model")
    error_3 <- 1
  }
  
  if(sum(is.na(match(temp, colnames(ref_dat)))) > 0)
  {
    print("names of covariates in the reference data does not 
      match with that of the study specific covariates")
    error_4 <- 1
  }
  
  ## Needed for calculating initial value
  names_wo_intercept <- c()
  if (is.null(estimated_intercept)){
    for(i in 1:no_of_studies)
    {
      names_wo_intercept <- union(names_wo_intercept,
        names(study_info[[i]][[1]][-1]))
    } 
  }else{
    for(i in 1:no_of_studies)
    {
      names_wo_intercept <- union(names_wo_intercept,
        names(study_info[[i]][[1]]))
    }
  }
  weight_sum <- 0
  sum <- 0
  estimates_in_which_studies_indices <- list()
  for(k in 1: length(names_wo_intercept))
  {
    temp_estimates_in_which <- c()
    for(j in 1: no_of_studies)
    {
      if(names_wo_intercept[k] %in% names(study_info[[j]][[1]]) == T)
        temp_estimates_in_which <- c(temp_estimates_in_which, j)
    }
    estimates_in_which_studies_indices[[k]] <- 
                                temp_estimates_in_which
  }
  ## End of the need for calculating initial value
  
  ## If the given inputs pass all the sanity checks...
  if(!(error_1 == 0 && error_2 == 0 && error_3 == 0 && error_4 == 0)){
    print("Please Fix the Error!")
    return(NULL)
  }
  ## Calculate initial_val if not provided
  if(length(initial_val) == 0)
  {
    initial_val <- c()
    
    # Calculating initial value when the variable_intercepts is FALSE
    if(variable_intercepts == F)
    {
      for(k in 1: length(names_wo_intercept))
      {
        for(j in estimates_in_which_studies_indices[[k]])
        {
          if(is.null(study_info[[j]][[2]]) == F)
          {
            index_cov <- which(names(study_info[[j]][[1]]) %in% names_wo_intercept[k]  == T)
            weight <- 1/study_info[[j]][[2]][index_cov, index_cov]
            weight_sum <- weight_sum + weight
            sum <- sum + study_info[[j]][[1]][index_cov]*weight
          }
          if(is.null(study_info[[j]][[2]]) == T)
          {
            index_cov <- which(names(study_info[[j]][[1]]) %in% names_wo_intercept[k]  == T)
            weight <- study_info[[j]][[3]]
            weight_sum <- weight_sum + weight
            sum <- sum + study_info[[j]][[1]][index_cov]*weight
          }
        }
        initial_val[k+1] <- sum/weight_sum
        weight_sum <- 0
        sum <- 0
      }
      if(is.null(estimated_intercept)){
        for(j in 1:no_of_studies)
        {
          if(is.null(study_info[[j]][[2]]) == F)
          {
            weight <- 1/study_info[[j]][[2]][1, 1]
            weight_sum <- weight_sum + weight
            sum <- sum + study_info[[j]][[1]][1]*weight
          }
          if(is.null(study_info[[j]][[2]]) == T)
          {
            weight <- study_info[[j]][[3]]
            weight_sum <- weight_sum + weight
            sum <- sum + study_info[[j]][[1]][1]*weight
          }
        }
        initial_val[1] <- sum/weight_sum
      }
    }
    if(is.na(initial_val[1])){
      initial_val<-initial_val[-1]
    }
    # End of Calculating initial value when the variable_intercepts is FALSE
  }
  initial_val1<<-initial_val
  print(length(initial_val))
  ## End of initial_val calculation

  
  ## Creating X_rbind and X_bdiag matrices
  
  ## Defining X_rbind here
  #X_rbind <- c()
  #for(k in 1 : no_of_studies)
  #{
  #  X_rbind <- rbind(X_rbind,ref_dat)
  #}
  ## Defining X_abdiag  here
  X_bdiag_list <- list()
  dim_C<- 0
  
  col_inds<-lapply(1:no_of_studies, function(k){
    which(colnames(ref_dat) %in% names(study_info[[k]][[1]]) == TRUE)
  })
  for(k in 1 : no_of_studies)
  {
    col_ind <-  col_inds[[k]]
    dim_C <- dim_C + length(col_ind)
    X_bdiag_list[[k]] <- as.matrix(ref_dat[,col_ind])
  }
  #X_bdiag <- Reduce(magic::adiag,X_bdiag_list)
  ## End of Creating X_rbind and X_bdiag matrices

  ########------End of Preparation-------########
  
  ## When the model is logistic, calls the newoptim function which implements the NR method...
  if(model == "logistic")
  {
    C_init <- diag(dim_C)
    total_iter <- 0
    eps<-1
    if(lambda.gam.first==0){
      D1<-NULL
    }else{
      D1<-D
    }
    output_initial <- useoptim(no_of_studies = no_of_studies, 
      study_info = study_info, 
      ref_dat = ref_dat, 
      col_inds = col_inds,
      X_bdiag_list = X_bdiag_list,
      C = C_init,
      initial_val = initial_val,
      threshold = threshold,
      model = model, 
      missing_covariance_study_indices = missing_covariance_study_indices, 
      variable_intercepts = variable_intercepts,
      return_C = TRUE,
      lambda.gam = lambda.gam.first,
      maxit=maxit,
      D = D1)
    coef_iter <- output_initial$estimated_coef
    # the variance will be considered later
    #asy_var_beta_converged_initial <- output_initial$Asy_var_optim
    total_iter<- total_iter + output_initial$no_of_iter
    C_iter <- output_initial$C
    cost_val_new<-output_initial$cost_val
    if(one_step==FALSE){
      while(total_iter < maxit & eps > threshold)
      {
        print("1")
        output_iter <- 
          useoptim(no_of_studies = no_of_studies, 
            study_info = study_info, 
            ref_dat = ref_dat, 
            col_inds = col_inds,
            X_bdiag_list = X_bdiag_list,
            C = C_iter,
            initial_val = coef_iter,
            threshold = threshold,
            model = model, 
            missing_covariance_study_indices = missing_covariance_study_indices, 
            variable_intercepts = variable_intercepts,
            return_C = TRUE,
            lambda.gam = lambda.gam.first,
            D = D1)
        coef_iter <- output_iter$estimated_coef
        C_iter <- output_iter$C
        cost_val_old<-cost_val_new
        cost_val_new<-output_iter$cost_val
        total_iter<- total_iter + output_iter$no_of_iter
        convergence<-output_iter$convergence
        eps<-abs(cost_val_new - cost_val_old)
      }
      
      output_iter <- 
        useoptim(no_of_studies = no_of_studies, 
          study_info = study_info, 
          ref_dat = ref_dat, 
          col_inds = col_inds,
          X_bdiag_list = X_bdiag_list,
          C = C_iter,
          initial_val = coef_iter,
          threshold = threshold,
          model = model, 
          missing_covariance_study_indices = missing_covariance_study_indices, 
          variable_intercepts = variable_intercepts,
          return_C = FALSE,
          return_Hessian = TRUE,
          lambda.gam = lambda.gam.first,
          D = D1)
      coef_iter <- output_iter$estimated_coef
      total_iter<- total_iter + output_iter$no_of_iter
      convergence<-output_iter$convergence
      res<-eigen(output_iter$Hessian)
    }else{
      print("use this")
      output_iter <- 
        useoptim(no_of_studies = no_of_studies, 
          study_info = study_info, 
          ref_dat = ref_dat, 
          col_inds = col_inds,
          X_bdiag_list = X_bdiag_list,
          C = C_iter,
          initial_val = coef_iter,
          threshold = threshold,
          model = model, 
          missing_covariance_study_indices = missing_covariance_study_indices, 
          variable_intercepts = variable_intercepts,
          return_C = FALSE,
          return_Hessian = TRUE,
          lambda.gam = lambda.gam.first,
          maxit = 100,
          D = D1)
      coef_iter <- output_iter$estimated_coef
      total_iter<- total_iter + output_iter$no_of_iter
      convergence<-output_iter$convergence
      res<-eigen(output_iter$Hessian)
      one_step_output_iter<<-output_iter
    }
    
    if (sum(Re(res$values)[which(abs(Re(res$values)) > threshold)] < 0) > 0 ){
      warning("convergence to Hessian singular result!")
      convergence<-TRUE
    }else{
      print(sum(Re(res$values)[which(abs(Re(res$values)) > threshold)] > 0))
      print("Hessian is good!!!")
      convergence<-TRUE
    }
    
    if (convergence){
        if(!is.null(lambda.gam)){
          lambdais0<-FALSE
          if(length(lambda.gam) == 1){
            if (lambda.gam == 0){
              lambdais0<-TRUE
            }
          }
          if (lambdais0){
              logistic_result<-list("estimated_coef" = coef_iter,
                "convergence"=convergence,
                "no_of_iter"=total_iter,
                "C_final"=C_iter)
          }else{
            print("go here")
            logistic_result <- mclapply(lambda.gam, function(lam){
              coef_iter_lam<-coef_iter
              output_lam <- 
                useoptim(no_of_studies = no_of_studies, 
                  study_info = study_info, 
                  ref_dat = ref_dat, 
                  col_inds = col_inds,
                  X_bdiag_list = X_bdiag_list,
                  C = C_iter,
                  initial_val = coef_iter_lam,
                  threshold = threshold,
                  model = model, 
                  missing_covariance_study_indices = missing_covariance_study_indices, 
                  variable_intercepts = variable_intercepts,
                  return_C = FALSE,
                  return_Hessian = TRUE,
                  lambda.gam = lam,
                  maxit=maxit,
                  D = D)
              coef_iter_lam<-output_lam$estimated_coef
              #output_lam$cost_val<-output_lam$cost_val*nrow(ref_dat)
              #lam<-lam*nrow(ref_dat)
              H<-output_lam$Hessian
              if (det(H + lam*D) > 0){
                REML<-c(output_lam$cost_val+log(det(H + lam*D))-sum(D)*log(lam) -(dim(D)[1]- sum(D))*log(2*pi)  )  
              }else{
                REML<-Inf
              }
  
              
              #print(paste0(lam,",",REML))
              
              
              GCV<-nrow(ref_dat)*(output_lam$cost_val-lam*coef_iter_lam%*%D%*%coef_iter_lam)/
                  (nrow(ref_dat)-sum(diag(crossprod(ref_dat,
                    ref_dat)%*%solve(crossprod(ref_dat,
                      ref_dat)+lam*nrow(ref_dat)*D))))^2
              list("estimated_coef" = coef_iter_lam,
                "marginal_likelihood"=REML, 
                "GCV"=GCV[1], 
                "convergence"=convergence,
                "convergence_lambda_tune"=output_lam$convergence,
                "no_of_iter"=total_iter,
                "no_of_iter_lambda_tune"=output_lam$no_of_iter) 
            },mc.cores = 4)
            
          }
          
        }else{
          logistic_result<-list("estimated_coef" = coef_iter,
            "convergence"=convergence,
            "no_of_iter"=total_iter,
            "C_final"=C_iter)
        }
      
    }else{
      logistic_result<-list("estimated_coef" = coef_iter,
        "convergence"=convergence,
        "no_of_iter"=total_iter,
        "C_final"=C_iter)
    }
    
    return(logistic_result)
    
  }
  
}
library(pracma)
FindIntercept<-function(beta,n1ratio,x){
  
  sigmoid<-function(vec){
    index_nonneg<-which(vec >= 0)
    index_neg<-which(vec < 0)
    res<-rep(NA,length(vec))
    res[index_nonneg]<-1/(1+exp(-vec[index_nonneg]))
    res[index_neg]<-exp(vec[index_neg])/(1+exp(vec[index_neg]))
    res
  }
  f<-function(alpha,beta,n1ratio,x,n){
    if(is.null(dim(x)[1])){
      logitp<-x*beta+alpha 
    }else{
      logitp<-c(as.matrix(x)%*%c(beta))+alpha 
    }
    sum(sigmoid(logitp))/n - n1ratio
  }
  gradf<-function(alpha,beta,x,n){
    if(is.null(dim(x)[1])){
      logitp<-x*beta+alpha 
    }else{
      logitp<-c(as.matrix(x)%*%c(beta))+alpha 
    }
    sum(sigmoid(logitp)*(1-sigmoid(logitp)))/n
  }
  finalf<-function(alpha){
    f(alpha,beta,n1ratio,x,n)
  }
  finalgradf<-function(alpha){
    gradf(alpha,beta,x,n)
  }
  if(is.null(dim(x)[1])){
    n<-length(x)
  }else{
    n<-dim(x)[1]
  }
  res<-newtonRaphson(fun=finalf,x0=log(n1ratio/(1-n1ratio)),dfun = finalgradf)
  return(res$root)
}

