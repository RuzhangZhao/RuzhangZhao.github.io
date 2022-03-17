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
  lambda_list<-lasso_initial$lambda
  
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list[50],alpha = 0)
  #ridge_fit<-glmnet(x= UKBB_pop[,-1],y= UKBB_pop[,1],standardize=F,intercept=F,lambda = lambda_list[100],alpha = 0)
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  #ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list[50],alpha = 0,penalty.factor = w_adaptive)
  ridge_fit<-glmnet(x= UKBB_pop[,-1],y= UKBB_pop[,1],standardize=F,intercept=F,lambda = lambda_list[100],alpha = 0,penalty.factor = w_adaptive)
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,penalty.factor = w_adaptive)
  lambda_list<-adaptivelasso_fit$lambda
  UKBB_cor<-cor(UKBB_pop[,-1])
  #nonzero_cor<-UKBB_cor[index_nonzero,index_nonzero]
  
  for( i in (which(lambda_list == adaptivelasso_fit$lambda.min)) : 2){
    beta_i<-coef(adaptivelasso_fit,s = lambda_list[i])[-1]
    index_nonzero_i<-which(abs(beta_i)>1e-3)
    index_nonzero_i<-which(beta_i!=0)
    nonzero_cor_i<-UKBB_cor[index_nonzero_i,index_nonzero_i]
    diag(nonzero_cor_i)<-0
    lambda_index<-i
    if( sum(nonzero_cor_i > 0.9) == 0 )break
  }
  lambda_index=which(lambda_list == adaptivelasso_fit$lambda.min)
  
  beta<-coef(adaptivelasso_fit,s = lambda_list[lambda_index])[-1]
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
  
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos2"=pos2,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive
    )
}
