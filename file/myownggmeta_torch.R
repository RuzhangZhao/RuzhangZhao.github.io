library(parallel)
ggmeta <- function(
  study_info, 
  ref_dat,
  model="logistics",
  variable_intercepts=FALSE,
  initial_val=NULL,
  lambda=NULL,
  D=NULL,
  return_asy_var=TRUE,
  initial_run_res = NULL,
  learning_rate = 0.01,
  lr_step_size = 100,
  EPOCH_inital = 1000,
  EPOCH = 200,
  penalty = 'L2',
  verbose=TRUE,
  verbose_step = 1)
{
  call_ggmeta <- match.call()

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
  for(i in 1:no_of_studies)
  {
    names_wo_intercept <- union(names_wo_intercept,
      names(study_info[[i]][[1]]))
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
  
  dim_C<- 0
  
  col_inds<-lapply(1:no_of_studies, function(k){
    which(colnames(ref_dat) %in% names(study_info[[k]][[1]]) == TRUE)
  })
  for(k in 1 : no_of_studies)
  {
    col_ind <-  col_inds[[k]]
    dim_C <- dim_C + length(col_ind)
  }
  
  if(is.null(initial_run_res)){
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
          initial_val[k] <- sum/weight_sum
          weight_sum <- 0
          sum <- 0
        }
      }
      # End of Calculating initial value when the variable_intercepts is FALSE
    }
    initial_val1<<-initial_val
    print(length(initial_val))
    ## End of initial_val calculation
    
    ########------End of Preparation-------########
    
    ## When the model is logistic, calls the newoptim function which implements the NR method...
    C_init <- diag(dim_C)
    total_iter <- 0
    eps<-1
    output_initial <- fastoptim(
      study_info = study_info, 
      ref_dat = ref_dat, 
      col_inds = col_inds,
      C = C_init,
      initial_val = initial_val,
      model = model, 
      missing_covariance_study_indices = missing_covariance_study_indices,
      learning_rate = learning_rate,
      lr_step_size = lr_step_size,
      EPOCH = EPOCH_inital,
      lam = 0,
      D = D,
      penalty = penalty,
      return_C = TRUE,
      return_asy_var = FALSE,
      verbose=verbose,
      verbose_step =verbose_step)
    coef_iter <- output_initial$estimated_coef
    # the variance will be considered later
    #asy_var_beta_converged_initial <- output_initial$Asy_var_optim
    total_iter<- total_iter + output_initial$no_of_iter
    C_iter <- output_initial$C
    cost_val_new<-output_initial$cost_val
    output_iter <- 
      fastoptim(
        study_info = study_info, 
        ref_dat = ref_dat, 
        col_inds = col_inds,
        C = C_iter,
        initial_val = coef_iter,
        model = model, 
        missing_covariance_study_indices = missing_covariance_study_indices,
        learning_rate = learning_rate,
        lr_step_size = lr_step_size,
        EPOCH = EPOCH,
        lam = 0,
        D = D,
        penalty = penalty,
        return_C = FALSE,
        return_asy_var = return_asy_var,
        verbose=verbose,
        verbose_step =verbose_step)
    coef_iter <- output_iter$estimated_coef
    total_iter<- total_iter + output_iter$no_of_iter
    if(!is.null(lambda)){
      lambdais0<-FALSE
      if(length(lambda) == 1){
        if (lambda == 0){
          lambdais0<-TRUE
        }
      }
      if (lambdais0){
        logistic_result<-list("estimated_coef" = coef_iter,
          "no_of_iter"=total_iter,
          "C_final"=C_iter)
        if(return_asy_var){
          logistic_result$asy_var<-output_iter$asy_var
        }
      }else{
        print("go here")
        logistic_result <- lapply(lambda, function(lam){
          coef_iter_lam<-coef_iter
          output_lam <- fastoptim(
            study_info = study_info, 
            ref_dat = ref_dat, 
            col_inds = col_inds,
            C = C_iter,
            initial_val = coef_iter_lam,
            model = model, 
            missing_covariance_study_indices = missing_covariance_study_indices,
            learning_rate = learning_rate,
            lr_step_size = lr_step_size,
            EPOCH = EPOCH,
            lam = lam,
            D = D,
            penalty = penalty,
            return_C = FALSE,
            return_asy_var = return_asy_var,
            verbose=verbose,
            verbose_step =verbose_step)
          coef_iter_lam<-output_lam$estimated_coef
          total_iter_lam<-total_iter+output_lam$no_of_iter
          GCV<-nrow(ref_dat)*(output_lam$cost_val)/
            (nrow(ref_dat)-sum(diag(crossprod(ref_dat,
              ref_dat)%*%solve(crossprod(ref_dat,
                ref_dat)+lam*nrow(ref_dat)*D))))^2
          GCV<-c(GCV)
          newList<-list("estimated_coef" = coef_iter_lam,
            "GCV"=GCV, 
            "no_of_iter"=total_iter_lam) 
          if(return_asy_var){
            newList$asy_var<-output_lam$asy_var
          }
        })
      }
    }else{
      logistic_result<-list("estimated_coef" = coef_iter,
        "no_of_iter"=total_iter,
        "C_final"=C_iter)
      if(return_asy_var){
        logistic_result$asy_var<-output_iter$asy_var
      }
    }
  }else{
    if(is.null(lambda)){
      return("lambda = 0 has been calculated as initial result!")
    }else{
      if(length(lambda) == 1){
        if (lambda == 0){
          return("lambda = 0 has been calculated as initial result!")
        }
      }
    }
    coef_iter<-initial_run_res$estimated_coef
    C_iter<-initial_run_res$C_final
    total_iter<-initial_run_res$no_of_iter
    
    logistic_result <- lapply(lambda, function(lam){
        coef_iter_lam<-coef_iter
        output_lam <- fastoptim(
        study_info = study_info, 
        ref_dat = ref_dat, 
        col_inds = col_inds,
        C = C_iter,
        initial_val = coef_iter_lam,
        model = model, 
        missing_covariance_study_indices = missing_covariance_study_indices,
        learning_rate = learning_rate,
        lr_step_size = lr_step_size,
        EPOCH = EPOCH,
        lam = lam,
        D = D,
        penalty = penalty,
        return_C = FALSE,
        return_asy_var = return_asy_var,
        verbose=verbose,
        verbose_step =verbose_step)
      coef_iter_lam<-output_lam$estimated_coef
      total_iter_lam<-total_iter+output_lam$no_of_iter
      GCV<-nrow(ref_dat)*(output_lam$cost_val)/
        (nrow(ref_dat)-sum(diag(crossprod(ref_dat,
          ref_dat)%*%solve(crossprod(ref_dat,
            ref_dat)+lam*nrow(ref_dat)*D))))^2
      GCV<-c(GCV)
      newList<-list("estimated_coef" = coef_iter_lam,
        "GCV"=GCV, 
        "no_of_iter"=total_iter_lam) 
      print(return_asy_var)
      if(return_asy_var){
        print("Yes")
        newList$asy_var<-output_lam$asy_var
        print(length(output_lam))
        print(output_lam$asy_var)
      }
      newList
    })
  }
  return(logistic_result)
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

