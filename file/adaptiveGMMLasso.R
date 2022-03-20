if(0){
  UKBB_pop<-UKBB_pop_all
  study_info<-study_info_scaled 
}
## UKBB_pop should be scaled with SNPs, 
## but not with phenotype.
## Both SNP and phenotype should be centered. 
library(expm,quietly = T)
library(magic,quietly = T)
adaptiveGMMlasso<-function(UKBB_pop,N_SNP,study_info){
  colnames(UKBB_pop)[1]<-"Y"
  var_SNP<-paste0("SNP",1:(N_SNP))
  len_SNP<-length(var_SNP)
  var_GPC<-NULL
  len_GPC<-length(var_GPC)
  var_nonSNP<-c(var_GPC)
  len_nonSNP<-length(var_nonSNP)
  var_names_full_fit <- c(var_SNP,var_nonSNP)
  colname_UKBB<-colnames(UKBB_pop)
  len_U1 = len_SNP+len_nonSNP
  len_U2 = len_SNP
  len_beta = len_U1
  len_theta = len_SNP+len_GPC
  len_U = len_U1 + len_U2
  N_Pop<-nrow(UKBB_pop)
  
  
  #### adaptive lasso with fast lasso computation 
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop#*var(UKBB_pop[,1])
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11)
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  
  ############## Beta
  pseudo_Xy<-function(
    C_half,UKBB_pop,
    var_SNP,var_GPC,
    N_SNP,theta_UKBB_GPC,study_info){
    N_Pop<-nrow(UKBB_pop)
    pseudo_X<-C_half%*%rbind(t(UKBB_pop[,-1]),t(UKBB_pop[,var_SNP]))%*%UKBB_pop[,-1]
    pseudo_y1<-crossprod(UKBB_pop[,1],UKBB_pop[,-1])
    if(0){
      pseudo_y2<-sapply(1:N_SNP, function(snp_id){
        u2_id<-UKBB_pop[,c("V",var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)
        c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
      })
    }
    pseudo_y22<-sapply(1:N_SNP, function(snp_id){
      u2_id<-UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
      c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
    })
    
    pseudo_y<-c(c(c(pseudo_y1),c(pseudo_y22))%*%C_half)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
  }
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  ridge_fit0<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01)
  #ridge_fit<-glmnet(x= UKBB_pop[,-1],y= UKBB_pop[,1],standardize=F,intercept=F,lambda = lambda_list[100],alpha = 0)
  gamma_adaptivelasso<-1/2
  w_adaptive0<-1/(abs(coef(ridge_fit0)[-1]))^gamma_adaptivelasso
  w_adaptive0[is.infinite(w_adaptive0)]<-max(w_adaptive0[!is.infinite(w_adaptive0)])*5
  #ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list[50],alpha = 0,penalty.factor = w_adaptive)
  ridge_fit<-glmnet(x= UKBB_pop[,-1],y= UKBB_pop[,1],standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01,penalty.factor = w_adaptive0)
  #ridge_fit<-glmnet(x= scale(ref),y= scale(pheno_EUR$T2D,scale = F),standardize=F,intercept=F,lambda = lambda_list[50]/100/mean(w_adaptive),penalty.factor = w_adaptive)
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*10
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  lambda_list<-adaptivelasso_fit$lambda
  UKBB_cor<-cor(UKBB_pop[,-1])
  #nonzero_cor<-UKBB_cor[index_nonzero,index_nonzero]
  if(0){
    for( i in (which(lambda_list == adaptivelasso_fit$lambda.min)) : 2){
      beta_i<-coef(adaptivelasso_fit,s = lambda_list[i])[-1]
      index_nonzero_i<-which(beta_i!=0)
      index_nonzero_i0<-which(abs(beta_i)>1e-3)
      nonzero_cor_i<-UKBB_cor[index_nonzero_i0,index_nonzero_i]
      nonzero_cor_i[nonzero_cor_i == 1] = 0
      lambda_index<-i
      if( sum(nonzero_cor_i > 0.8) == 0 )break
    }
  }
  if(1){
    for( i in length(adaptivelasso_fit$lambda) : 2){
      beta_i<-coef(adaptivelasso_fit,s = lambda_list[i])[-1]
      index_nonzero_i<-which(beta_i!=0)
      index_nonzero_i0<-which(abs(beta_i)>1e-3)
      nonzero_cor_i<-UKBB_cor[index_nonzero_i0,index_nonzero_i]
      nonzero_cor_i[nonzero_cor_i == 1] = 0
      lambda_index<-i
      if( sum(nonzero_cor_i > 0.75) == 0 )break
    }
  }
  beta<-coef(adaptivelasso_fit,s =lambda_list[lambda_index])[-1]
  #index_nonzero_i0<-which(abs(beta)>1e-3)
  
  index_nonzero<-which(beta!=0)
  xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
  Sigsum_half<-cbind(xtx,xtx)%*%C_half
  Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
  Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
  inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
  
  W1 = xtx/var(UKBB_pop[,1])
  W2 = xtx%*%C_22%*%xtx
  W = W1+W2
  W_nonzero = W[index_nonzero,index_nonzero]
  
  final_v<-diag(inv_Sigsum_scaled_nonzero%*%W_nonzero%*%inv_Sigsum_scaled_nonzero)
  
  aa_final<-1-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1)
  pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  pos2<-index_nonzero[which(aa_final<0.05/length(aa_final))]
  #print(pos2)
  #print(beta[index_nonzero])
  #print(final_v)
  #print(index_nonzero)
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos2"=pos2,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive
    )
}

## UKBB_pop should be scaled with SNPs, 
## but not with phenotype.
## Both SNP and phenotype should be centered. 
library(expm,quietly = T)
library(magic,quietly = T)
adaptiveGMMlasso2<-function(UKBB_pop,N_SNP,study_info){
  colnames(UKBB_pop)[1]<-"Y"
  var_SNP<-paste0("SNP",1:(N_SNP))
  len_SNP<-length(var_SNP)
  var_GPC<-NULL
  len_GPC<-length(var_GPC)
  var_nonSNP<-c(var_GPC)
  len_nonSNP<-length(var_nonSNP)
  var_names_full_fit <- c(var_SNP,var_nonSNP)
  colname_UKBB<-colnames(UKBB_pop)
  len_U1 = len_SNP+len_nonSNP
  len_U2 = len_SNP
  len_beta = len_U1
  len_theta = len_SNP+len_GPC
  len_U = len_U1 + len_U2
  N_Pop<-nrow(UKBB_pop)
  
  
  #### adaptive lasso with fast lasso computation 
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop#*var(UKBB_pop[,1])
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11)
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  #C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  
  ############## Beta
  pseudo_Xy<-function(
    C_half,UKBB_pop,
    var_SNP,var_GPC,
    N_SNP,theta_UKBB_GPC,study_info){
    N_Pop<-nrow(UKBB_pop)
    pseudo_X<-C_half%*%rbind(t(UKBB_pop[,-1]),t(UKBB_pop[,var_SNP]))%*%UKBB_pop[,-1]
    pseudo_y1<-crossprod(UKBB_pop[,1],UKBB_pop[,-1])
    if(0){
      pseudo_y2<-sapply(1:N_SNP, function(snp_id){
        u2_id<-UKBB_pop[,c("V",var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)
        c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
      })
    }
    pseudo_y22<-sapply(1:N_SNP, function(snp_id){
      u2_id<-UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
      c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
    })
    
    pseudo_y<-c(c(c(pseudo_y1),c(pseudo_y22))%*%C_half)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
  }
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  
  ridge_fit0<-glmnet(x= UKBB_pop[,-1],y= UKBB_pop[,1],standardize=F,intercept=F,lambda = lambda_list0[50]/10,alpha = 0.01)
  #ridge_fit<-glmnet(x= UKBB_pop[,-1],y= UKBB_pop[,1],standardize=F,intercept=F,lambda = lambda_list[100],alpha = 0)
  gamma_adaptivelasso<-1/2
  w_adaptive0<-1/(abs(coef(ridge_fit0)[-1]))^gamma_adaptivelasso
  w_adaptive0[is.infinite(w_adaptive0)]<-max(w_adaptive0[!is.infinite(w_adaptive0)])*5
  #ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list[50],alpha = 0,penalty.factor = w_adaptive)
  #ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/10,alpha = 0.01,penalty.factor = w_adaptive0)
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/10,alpha = 0.01)
  #ridge_fit<-glmnet(x= scale(ref),y= scale(pheno_EUR$T2D,scale = F),standardize=F,intercept=F,lambda = lambda_list[50]/100/mean(w_adaptive),penalty.factor = w_adaptive)
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*10
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  lambda_list<-adaptivelasso_fit$lambda
  UKBB_cor<-cor(UKBB_pop[,-1])
  beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
  
  
  while(1){
    index_nonzero_i<-which(beta!=0)
    nonzero_cor_i<-UKBB_cor[index_nonzero_i,index_nonzero_i]
    nonzero_cor_i[nonzero_cor_i == 1] = 0
    tmp<-which(abs(nonzero_cor_i)>0.8,arr.ind = T)
    if(nrow(tmp)>0){
      tmp<-cbind(index_nonzero_i[tmp[,2]],index_nonzero_i[tmp[,1]])
      tmp<-tmp[tmp[,1]<tmp[,2],]
      #which(index_nonzero_i[tmp[,2]]%in%index_nonzero_i0[tmp[,1]])
      rm_index<-c()
      if(is.null(dim(tmp)[1])){
        z1<-study_info[[tmp[1]]]$Coeff^2/study_info[[tmp[1]]]$Covariance
        z2<-study_info[[tmp[2]]]$Coeff^2/study_info[[tmp[2]]]$Covariance
        if( z1 > z2){
          rm_index<-c(rm_index,tmp[2])
        }else{
          rm_index<-c(rm_index,tmp[1])
        }
      }else{
        for(i in 1:nrow(tmp)){
          z1<-study_info[[tmp[i,1]]]$Coeff^2/study_info[[tmp[i,1]]]$Covariance
          z2<-study_info[[tmp[i,2]]]$Coeff^2/study_info[[tmp[i,2]]]$Covariance
          if( z1 > z2){
            rm_index<-c(rm_index,tmp[i,2])
          }else{
            rm_index<-c(rm_index,tmp[i,1])
          }
        }
      }
      
      w_adaptive<-w_adaptive
      w_adaptive[rm_index]<-w_adaptive[rm_index]*100
      
      adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
      beta<-coef(adaptivelasso_fit,s = 'lambda.min')[-1]
      
    }else{
      break
    } 
  }
  
  
  index_nonzero<-which(beta!=0)
  xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
  Sigsum_half<-cbind(xtx,xtx)%*%C_half
  Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
  Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
  inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
  
  W1 = xtx/var(UKBB_pop[,1])
  W2 = xtx%*%C_22%*%xtx
  W = W1+W2
  W_nonzero = W[index_nonzero,index_nonzero]
  
  final_v<-diag(inv_Sigsum_scaled_nonzero%*%W_nonzero%*%inv_Sigsum_scaled_nonzero)
  
  aa_final<-1-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1)
  pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  pos2<-index_nonzero[which(aa_final<0.05/length(aa_final))]
  #print(pos)
  #print(beta[index_nonzero])
  #print(final_v)
  #print(index_nonzero)
  print(paste0("original_MSE2:",sum((UKBB_pop[,-1]%*%beta-UKBB_pop[,1])^2)))
  print(paste0("pseudo_MSE2:",sum((pseudo_X%*%beta-pseudo_y)^2)))
  
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos2"=pos2,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive
  )
}


############### 
## 3 is set prior from unbalanced (almost from GWAS)

library(expm,quietly = T)
library(magic,quietly = T)
adaptiveGMMlasso3<-function(UKBB_pop,N_SNP,study_info){
  colnames(UKBB_pop)[1]<-"Y"
  var_SNP<-paste0("SNP",1:(N_SNP))
  len_SNP<-length(var_SNP)
  var_GPC<-NULL
  len_GPC<-length(var_GPC)
  var_nonSNP<-c(var_GPC)
  len_nonSNP<-length(var_nonSNP)
  var_names_full_fit <- c(var_SNP,var_nonSNP)
  colname_UKBB<-colnames(UKBB_pop)
  len_U1 = len_SNP+len_nonSNP
  len_U2 = len_SNP
  len_beta = len_U1
  len_theta = len_SNP+len_GPC
  len_U = len_U1 + len_U2
  N_Pop<-nrow(UKBB_pop)
  
  
  #### adaptive lasso with fast lasso computation 
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop#*var(UKBB_pop[,1])
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11)
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  
  ############## Beta
  pseudo_Xy<-function(
    C_half,UKBB_pop,
    var_SNP,var_GPC,
    N_SNP,theta_UKBB_GPC,study_info){
    N_Pop<-nrow(UKBB_pop)
    pseudo_X<-C_half%*%rbind(t(UKBB_pop[,-1]),t(UKBB_pop[,var_SNP]))%*%UKBB_pop[,-1]
    pseudo_y1<-crossprod(UKBB_pop[,1],UKBB_pop[,-1])
    if(0){
      pseudo_y2<-sapply(1:N_SNP, function(snp_id){
        u2_id<-UKBB_pop[,c("V",var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)
        c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
      })
    }
    pseudo_y22<-sapply(1:N_SNP, function(snp_id){
      u2_id<-UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
      c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
    })
    
    pseudo_y<-c(c(c(pseudo_y1),c(pseudo_y22))%*%C_half)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
  }
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  
  ridge_fit0<-glmnet(x= UKBB_pop[,-1],y= UKBB_pop[,1],standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01)
  #ridge_fit<-glmnet(x= UKBB_pop[,-1],y= UKBB_pop[,1],standardize=F,intercept=F,lambda = lambda_list[100],alpha = 0)
  gamma_adaptivelasso<-1/2
  w_adaptive0<-1/(abs(coef(ridge_fit0)[-1]))^gamma_adaptivelasso
  w_adaptive0[is.infinite(w_adaptive0)]<-max(w_adaptive0[!is.infinite(w_adaptive0)])*5
  #ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list[50],alpha = 0,penalty.factor = w_adaptive)
  #ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/10,alpha = 0.01,penalty.factor = w_adaptive0)
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01,penalty.factor = w_adaptive0)
  #ridge_fit<-glmnet(x= scale(ref),y= scale(pheno_EUR$T2D,scale = F),standardize=F,intercept=F,lambda = lambda_list[50]/100/mean(w_adaptive),penalty.factor = w_adaptive)
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*10
  
  
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  lambda_list<-adaptivelasso_fit$lambda
  UKBB_cor<-cor(UKBB_pop[,-1])
  beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
  
  
  while(1){
    index_nonzero_i<-which(beta!=0)
    nonzero_cor_i<-UKBB_cor[index_nonzero_i,index_nonzero_i]
    nonzero_cor_i[nonzero_cor_i == 1] = 0
    tmp<-which(abs(nonzero_cor_i)>0.75,arr.ind = T)
    if(nrow(tmp)>0){
      tmp<-cbind(index_nonzero_i[tmp[,2]],index_nonzero_i[tmp[,1]])
      tmp<-tmp[tmp[,1]<tmp[,2],]
      #which(index_nonzero_i[tmp[,2]]%in%index_nonzero_i0[tmp[,1]])
      rm_index<-c()
      if(is.null(dim(tmp)[1])){
        z1<-study_info[[tmp[1]]]$Coeff^2/study_info[[tmp[1]]]$Covariance
        z2<-study_info[[tmp[2]]]$Coeff^2/study_info[[tmp[2]]]$Covariance
        if( z1 > z2){
          rm_index<-c(rm_index,tmp[2])
        }else{
          rm_index<-c(rm_index,tmp[1])
        }
      }else{
        for(i in 1:nrow(tmp)){
          z1<-study_info[[tmp[i,1]]]$Coeff^2/study_info[[tmp[i,1]]]$Covariance
          z2<-study_info[[tmp[i,2]]]$Coeff^2/study_info[[tmp[i,2]]]$Covariance
          if( z1 > z2){
            rm_index<-c(rm_index,tmp[i,2])
          }else{
            rm_index<-c(rm_index,tmp[i,1])
          }
        }
      }
      
      w_adaptive<-w_adaptive
      w_adaptive[rm_index]<-w_adaptive[rm_index]*100
      
      adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
      beta<-coef(adaptivelasso_fit,s = 'lambda.min')[-1]
      
    }else{
      break
    } 
  }
  
  
  index_nonzero<-which(beta!=0)
  xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
  Sigsum_half<-cbind(xtx,xtx)%*%C_half
  Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
  Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
  inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
  
  W1 = xtx/var(UKBB_pop[,1])
  W2 = xtx%*%C_22%*%xtx
  W = W1+W2
  W_nonzero = W[index_nonzero,index_nonzero]
  
  final_v<-diag(inv_Sigsum_scaled_nonzero%*%W_nonzero%*%inv_Sigsum_scaled_nonzero)
  
  aa_final<-1-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1)
  pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  pos2<-index_nonzero[which(aa_final<0.05/length(aa_final))]
  #print(pos)
  #print(beta[index_nonzero])
  #print(final_v)
  #print(index_nonzero)
  print(paste0("original_MSE3:",sum((UKBB_pop[,-1]%*%beta-UKBB_pop[,1])^2)))
  print(paste0("pseudo_MSE3:",sum((pseudo_X%*%beta-pseudo_y)^2)))
  
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos2"=pos2,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive
  )
}


adaptiveGMMlasso31<-function(UKBB_pop,N_SNP,study_info,type=1){
  colnames(UKBB_pop)[1]<-"Y"
  var_SNP<-paste0("SNP",1:(N_SNP))
  len_SNP<-length(var_SNP)
  var_GPC<-NULL
  len_GPC<-length(var_GPC)
  var_nonSNP<-c(var_GPC)
  len_nonSNP<-length(var_nonSNP)
  var_names_full_fit <- c(var_SNP,var_nonSNP)
  colname_UKBB<-colnames(UKBB_pop)
  len_U1 = len_SNP+len_nonSNP
  len_U2 = len_SNP
  len_beta = len_U1
  len_theta = len_SNP+len_GPC
  len_U = len_U1 + len_U2
  N_Pop<-nrow(UKBB_pop)
  
  
  
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop#*var(UKBB_pop[,1])
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11)
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  #C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  
  ############## Beta
  pseudo_Xy<-function(
    C_half,UKBB_pop,
    var_SNP,var_GPC,
    N_SNP,theta_UKBB_GPC,study_info){
    N_Pop<-nrow(UKBB_pop)
    pseudo_X<-C_half%*%rbind(t(UKBB_pop[,-1]),t(UKBB_pop[,var_SNP]))%*%UKBB_pop[,-1]
    pseudo_y1<-crossprod(UKBB_pop[,1],UKBB_pop[,-1])
    if(0){
      pseudo_y2<-sapply(1:N_SNP, function(snp_id){
        u2_id<-UKBB_pop[,c("V",var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)
        c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
      })
    }
    pseudo_y22<-sapply(1:N_SNP, function(snp_id){
      u2_id<-UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
      c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
    })
    
    pseudo_y<-c(c(c(pseudo_y1),c(pseudo_y22))%*%C_half)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
  }
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  
  ridge_fit0<-glmnet(x= UKBB_pop[,-1],y= UKBB_pop[,1],standardize=F,intercept=F,lambda = lambda_list0[50]/10,alpha = 0.01)
  #ridge_fit<-glmnet(x= UKBB_pop[,-1],y= UKBB_pop[,1],standardize=F,intercept=F,lambda = lambda_list[100],alpha = 0)
  gamma_adaptivelasso<-1/2
  w_adaptive0<-1/(abs(coef(ridge_fit0)[-1]))^gamma_adaptivelasso
  w_adaptive0[is.infinite(w_adaptive0)]<-max(w_adaptive0[!is.infinite(w_adaptive0)])*5
  #ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list[50],alpha = 0,penalty.factor = w_adaptive)
  #ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/10,alpha = 0.01,penalty.factor = w_adaptive0)
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/10,alpha = 0.01)
  #ridge_fit<-glmnet(x= scale(ref),y= scale(pheno_EUR$T2D,scale = F),standardize=F,intercept=F,lambda = lambda_list[50]/100/mean(w_adaptive),penalty.factor = w_adaptive)
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*10
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  lambda_list<-adaptivelasso_fit$lambda
  UKBB_cor<-cor(UKBB_pop[,-1])
  beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
  
  
  while(1){
    index_nonzero_i<-which(beta!=0)
    nonzero_cor_i<-UKBB_cor[index_nonzero_i,index_nonzero_i]
    nonzero_cor_i[nonzero_cor_i == 1] = 0
    tmp<-which(abs(nonzero_cor_i)>0.8,arr.ind = T)
    if(nrow(tmp)>0){
      tmp<-cbind(index_nonzero_i[tmp[,2]],index_nonzero_i[tmp[,1]])
      tmp<-tmp[tmp[,1]<tmp[,2],]
      #which(index_nonzero_i[tmp[,2]]%in%index_nonzero_i0[tmp[,1]])
      rm_index<-c()
      if(is.null(dim(tmp)[1])){
        z1<-study_info[[tmp[1]]]$Coeff^2/study_info[[tmp[1]]]$Covariance
        z2<-study_info[[tmp[2]]]$Coeff^2/study_info[[tmp[2]]]$Covariance
        if( z1 > z2){
          rm_index<-c(rm_index,tmp[2])
        }else{
          rm_index<-c(rm_index,tmp[1])
        }
      }else{
        for(i in 1:nrow(tmp)){
          z1<-study_info[[tmp[i,1]]]$Coeff^2/study_info[[tmp[i,1]]]$Covariance
          z2<-study_info[[tmp[i,2]]]$Coeff^2/study_info[[tmp[i,2]]]$Covariance
          if( z1 > z2){
            rm_index<-c(rm_index,tmp[i,2])
          }else{
            rm_index<-c(rm_index,tmp[i,1])
          }
        }
      }
      
      w_adaptive<-w_adaptive
      w_adaptive[rm_index]<-w_adaptive[rm_index]*100
      
      adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
      beta<-coef(adaptivelasso_fit,s = 'lambda.min')[-1]
      
    }else{
      break
    } 
  }
  
  beta1<-beta
  
  #### adaptive lasso with fast lasso computation 
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop#*var(UKBB_pop[,1])
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11)
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  
  ############## Beta
  pseudo_Xy<-function(
    C_half,UKBB_pop,
    var_SNP,var_GPC,
    N_SNP,theta_UKBB_GPC,study_info){
    N_Pop<-nrow(UKBB_pop)
    pseudo_X<-C_half%*%rbind(t(UKBB_pop[,-1]),t(UKBB_pop[,var_SNP]))%*%UKBB_pop[,-1]
    pseudo_y1<-crossprod(UKBB_pop[,1],UKBB_pop[,-1])
    if(0){
      pseudo_y2<-sapply(1:N_SNP, function(snp_id){
        u2_id<-UKBB_pop[,c("V",var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)
        c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
      })
    }
    pseudo_y22<-sapply(1:N_SNP, function(snp_id){
      u2_id<-UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
      c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
    })
    
    pseudo_y<-c(c(c(pseudo_y1),c(pseudo_y22))%*%C_half)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
  }
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  
  ridge_fit0<-glmnet(x= UKBB_pop[,-1],y= UKBB_pop[,1],standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01)
  #ridge_fit<-glmnet(x= UKBB_pop[,-1],y= UKBB_pop[,1],standardize=F,intercept=F,lambda = lambda_list[100],alpha = 0)
  gamma_adaptivelasso<-1/2
  w_adaptive0<-1/(abs(coef(ridge_fit0)[-1]))^gamma_adaptivelasso
  w_adaptive0[is.infinite(w_adaptive0)]<-max(w_adaptive0[!is.infinite(w_adaptive0)])*5
  #ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list[50],alpha = 0,penalty.factor = w_adaptive)
  #ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/10,alpha = 0.01,penalty.factor = w_adaptive0)
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01,penalty.factor = w_adaptive0)
  #ridge_fit<-glmnet(x= scale(ref),y= scale(pheno_EUR$T2D,scale = F),standardize=F,intercept=F,lambda = lambda_list[50]/100/mean(w_adaptive),penalty.factor = w_adaptive)
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*10
  
  
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  lambda_list<-adaptivelasso_fit$lambda
  UKBB_cor<-cor(UKBB_pop[,-1])
  beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
  
  
  while(1){
    index_nonzero_i<-which(beta!=0)
    nonzero_cor_i<-UKBB_cor[index_nonzero_i,index_nonzero_i]
    nonzero_cor_i[nonzero_cor_i == 1] = 0
    tmp<-which(abs(nonzero_cor_i)>0.75,arr.ind = T)
    if(nrow(tmp)>0){
      tmp<-cbind(index_nonzero_i[tmp[,2]],index_nonzero_i[tmp[,1]])
      tmp<-tmp[tmp[,1]<tmp[,2],]
      #which(index_nonzero_i[tmp[,2]]%in%index_nonzero_i0[tmp[,1]])
      rm_index<-c()
      if(is.null(dim(tmp)[1])){
        z1<-study_info[[tmp[1]]]$Coeff^2/study_info[[tmp[1]]]$Covariance
        z2<-study_info[[tmp[2]]]$Coeff^2/study_info[[tmp[2]]]$Covariance
        if( z1 > z2){
          rm_index<-c(rm_index,tmp[2])
        }else{
          rm_index<-c(rm_index,tmp[1])
        }
      }else{
        for(i in 1:nrow(tmp)){
          z1<-study_info[[tmp[i,1]]]$Coeff^2/study_info[[tmp[i,1]]]$Covariance
          z2<-study_info[[tmp[i,2]]]$Coeff^2/study_info[[tmp[i,2]]]$Covariance
          if( z1 > z2){
            rm_index<-c(rm_index,tmp[i,2])
          }else{
            rm_index<-c(rm_index,tmp[i,1])
          }
        }
      }
      
      w_adaptive<-w_adaptive
      w_adaptive[rm_index]<-w_adaptive[rm_index]*100
      
      adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
      beta<-coef(adaptivelasso_fit,s = 'lambda.min')[-1]
      
    }else{
      break
    } 
  }
  beta2<-beta
  
  mse1_uk<-mean((UKBB_pop[,-1]%*%beta1-UKBB_pop[,1])^2)
  mse1_ps<-mean((pseudo_X%*%beta1-pseudo_y)^2)
  
  
  mse2_uk<-mean((UKBB_pop[,-1]%*%beta2-UKBB_pop[,1])^2)
  mse2_ps<-mean((pseudo_X%*%beta2-pseudo_y)^2)
  if(type == 1){
    if(mse1_ps<mse2_ps){
      beta<-beta1
    }else{
      beta<-beta2
    }
  }else if(type == 2){
    if(mse1_uk<mse2_uk){
      beta<-beta1
    }else{
      beta<-beta2
    }
  }else if(type == 3){
    
    if(mse1_uk<mse2_uk & mse1_ps<mse2_ps){
      beta<-beta1
    }else if(mse1_uk>mse2_uk & mse1_ps>mse2_ps){
      beta<-beta2
    }else if( mse1_uk>mse2_uk & mse1_ps<mse2_ps){
      if(mse1_uk/mse2_uk > mse2_ps/mse1_ps ){
        beta<-beta2
      }else{
        beta<-beta1
      }
    }else{
      if(mse2_uk/mse1_uk > mse1_ps/mse2_ps ){
        beta<-beta1
      }else{
        beta<-beta2
      }
    }
    
    
  }
  
  index_nonzero<-which(beta!=0)
  xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
  Sigsum_half<-cbind(xtx,xtx)%*%C_half
  Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
  Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
  inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
  
  W1 = xtx/var(UKBB_pop[,1])
  W2 = xtx%*%C_22%*%xtx
  W = W1+W2
  W_nonzero = W[index_nonzero,index_nonzero]
  
  final_v<-diag(inv_Sigsum_scaled_nonzero%*%W_nonzero%*%inv_Sigsum_scaled_nonzero)
  
  aa_final<-1-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1)
  pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  pos2<-index_nonzero[which(aa_final<0.05/length(aa_final))]
  #print(pos)
  #print(beta[index_nonzero])
  #print(final_v)
  #print(index_nonzero)
  
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos2"=pos2,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive)
}



