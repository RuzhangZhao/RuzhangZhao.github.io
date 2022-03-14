
print("JSgoodhhhh")

#Nonnull_index<-c(2,130,173)
#Nonnull_index<-c(2,130,192)
#Nonnull_index<-c(139,151,211)

#Nonnull_index<-c(151,139,103)
#Nonnull_index<-sample(1:254,3)
#Nonnull_index<-c(11,115,227)#
Nonnull_index<-c(25,83,196)
library(inline,quietly = T)
library(data.table,quietly = T)
library(dplyr,quietly = T)
#library(speedglm,quietly = T)
library(Rfast,quietly = T)
#library(feather,quietly = T)
library(locfit,quietly = T)
library(stringr,quietly = T)
library(selectiveInference,quietly = T)
plus_number<-sample(1:1e4,1)
foldpath<-"/dcl01/chatterj/data/rzhao/T2D_UKBB/"

library(MultiPhen,quietly = T)
#region_name<-"T2D_FTO"
library(ggplot2,quietly = T)
library(pracma,quietly = T)
library(Matrix,quietly = T)
library(genio,quietly = T)
library(glmnet,quietly = T)

if(1){

region_name<-"T2D_FTO2"
# load whole SNP data
.<-capture.output(ref<-read.plink(paste0(foldpath,region_name)))
fam <-read.table(paste0(foldpath,region_name,".fam"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
bim <-read.table(paste0(foldpath,region_name,".bim"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

toomuchna<-rowSums(is.na(ref))
na_index<-which(toomuchna>ncol(ref)/20)
ref<-ref[-na_index,]
fam<-fam[-na_index,]
ref[is.na(ref)]<-1
duplicated_index<-which(duplicated(t(ref))==TRUE)
ref<-ref[,-duplicated_index]
bim<-bim[-duplicated_index,]

### Remove duplicated SNP

## Try to separate highly correlated SNPs. 

#cor_ref<-cor(ref)
cor_ref<-readRDS(paste0(foldpath,"cor_ref_FTO.rds"))
#library(ComplexHeatmap)
#Heatmap(cor_ref)
cor_ref_cutoff<-cor_ref
cor_ref_cutoff[which(cor_ref_cutoff<.95)]<-0



## Get the filter version of some SNPs

##### Get pairs for each SNPs
cur_pairs<-list()
for (i in 1:(dim(cor_ref_cutoff)[1]-1) ){
  cur_pairs[[i]]<-c(i)
  for (j in (i+1):(dim(cor_ref_cutoff)[1]) ){
    if(cor_ref_cutoff[i,j]!=0){
      cur_pairs[[i]]<-c(cur_pairs[[i]],j)
    }
  }
}
cur_pairs[[dim(cor_ref_cutoff)[1]]]<-c(dim(cor_ref_cutoff)[1])
##### Filter SNPs by dynamic programming
mylist<-cur_pairs
for(i in 1:length(mylist)){
  index1<-c()
  for(j in 1:length(mylist)){
    if(i %in% mylist[[j]]){
      index1<-c(index1,j) 
    }
  }
  if(length(index1)>1){
    for(k in 2:length(index1)){
      mylist[[index1[1]]]<-sort(union(mylist[[index1[1]]],mylist[[index1[k]]]))
      mylist[[index1[k]]]<-NA
    }
  }
}

##### Get final pairs without NA
cur_pairs_final<-list()
count_final<-1
for(i in 1:length(mylist)){
  if(!is.na(mylist[[i]][1])){
    cur_pairs_final[[count_final]]<-mylist[[i]]
    count_final<-count_final+1
  }
}

##############################################
########## Set Nonnull SNPs 

#Nonnull_index<-sample(1:N_SNP,3)
#Nonnull_index<-c(2,130,192)

#c(2,130,252) perfect #c(74,130,259) bad 
#c(2,130,147) #c(130,147)
#a<-cor(ref[,Nonnull_index],ref[,-Nonnull_index])
#sum(a>0.9)

##### Get unique SNP out of the list

filter_SNP_vec<-c()
for(i in 1:length(cur_pairs_final)){
  if(length(cur_pairs_final[[i]]) == 1){
    filter_SNP_vec<-c(filter_SNP_vec,cur_pairs_final[[i]])
  }else{
    cor_small<-cor_ref[cur_pairs_final[[i]],cur_pairs_final[[i]]]
    diag(cor_small)<-0
    filter_SNP_vec<-c(filter_SNP_vec,cur_pairs_final[[i]][which.max(rowSums(cor_small))])
  }
}

filter_SNP_vec<-sort(filter_SNP_vec)
#Nonnull_index_filter_less<-which(filter_SNP_vec%in%Nonnull_index_filter)

##################################
######## Whether Shrinkage########

ref<-ref[,filter_SNP_vec]
}
N_SNP<-ncol(ref)
filter_SNP_vec<-1:N_SNP
Nonnull_index_filter_less<-Nonnull_index

coef_nonnull<-log(1.25)
coef_SNP<-rep(0,ncol(ref))
coef_SNP[Nonnull_index]<-coef_nonnull
coef_SNP[Nonnull_index]<-coef_SNP[Nonnull_index]*c(1,1,1)
intercept<- -2

EAF<-colsums(ref)/nrow(ref)/2
EAF[Nonnull_index]
cur_iter<-1
  set.seed(cur_iter)
  cur_iter<-plus_number+cur_iter
  
  SNP_Pop_matrix<-ref
  SNP_Pop_matrix<-t(t(SNP_Pop_matrix)-colMeans(SNP_Pop_matrix)) 
  genetic_var<-c(coef_SNP^2%*%Rfast::colVars(SNP_Pop_matrix))
  print(paste0("genetic variance:",round(genetic_var,5) ))
  
  Phenotype_logitp<-SNP_Pop_matrix%*%coef_SNP + intercept
  Phenotype_p<-1/(1+exp(-Phenotype_logitp ))
  print(max(Phenotype_logitp))
  print(min(Phenotype_logitp))
  Phenotype<-rbinom(nrow(ref), size=1, Phenotype_p) 
  print(paste0('prevalence',sum(Phenotype)/length(Phenotype)))
  
  
  library(caret)
  index1<-createDataPartition(Phenotype,p = 0.06)[[1]]
  ref_sample<-ref[index1,filter_SNP_vec]
  ref<-ref[-index1,filter_SNP_vec]
  N_Pop<-nrow(ref_sample)
  N_SNP<-ncol(ref_sample)
  EUR<-fam[-index1,1:2]
  write.table(EUR, paste0("/dcl01/chatterj/data/rzhao/T2D_UKBB/sim_GWAS/EUR",cur_iter,".eid"),row.names =F,col.names = F)
  pheno_EUR <- data.frame(FID=fam[-index1,1], IID=fam[-index1,2],T2D=Phenotype[-index1]+1)
  readr::write_tsv(pheno_EUR, paste0("/dcl01/chatterj/data/rzhao/T2D_UKBB/sim_GWAS/pheno_EUR",cur_iter,".txt"))
  EUR_SNP<-data.frame(SNP=bim$V2[filter_SNP_vec])
  readr::write_tsv(EUR_SNP, paste0("/dcl01/chatterj/data/rzhao/T2D_UKBB/sim_GWAS/EUR_SNP",cur_iter,".txt"),col_names = F)
  
  bim_ref<-bim[filter_SNP_vec,]
  colnames(bim_ref)<-c("chr", "id", "posg", "pos", "ref", "alt")
  colnames(ref)<-bim_ref$id
  fam_ref<-fam
  fam_ref<-fam_ref[-index1,]
  colnames(fam_ref)<-c("fam", "id", "pat", "mat", "sex", "pheno")
  write_plink(paste0("/dcl01/chatterj/data/rzhao/T2D_UKBB/sim_GWAS/ukbb_pop_all",cur_iter),as.matrix(t(ref)),
    bim =  bim_ref,fam = fam_ref)
  
  arg= paste0("/dcl01/chatterj/data/tools/plink2  --extract /dcl01/chatterj/data/rzhao/T2D_UKBB/sim_GWAS/EUR_SNP",cur_iter,".txt --keep /dcl01/chatterj/data/rzhao/T2D_UKBB/sim_GWAS/EUR",cur_iter,".eid  --bfile /dcl01/chatterj/data/rzhao/T2D_UKBB/sim_GWAS/ukbb_pop_all",cur_iter," --snps-only --pheno /dcl01/chatterj/data/rzhao/T2D_UKBB/sim_GWAS/pheno_EUR",cur_iter,".txt --glm hide-covar --vif 10000 --covar-variance-standardize --out /dcl01/chatterj/data/rzhao/T2D_UKBB/sim_GWAS/GWAS_T2D_FTO_",cur_iter)
  system(arg,ignore.stdout = T)
  
  
  gwas_res<-read.table(paste0("/dcl01/chatterj/data/rzhao/T2D_UKBB/sim_GWAS/GWAS_T2D_FTO_",cur_iter,".T2D.glm.logistic"))
  colnames(gwas_res)<-c("CHROM","POS","ID","REF","ALT","TEST","OBS_CT","BETA","SE","T_STAT","P")
  
  gwas_res$BETA<-log(gwas_res$BETA)
  
  
  sim_GWAS_path<-"/dcl01/chatterj/data/rzhao/T2D_UKBB/sim_GWAS/"
  .<-capture.output(file.remove(paste0(sim_GWAS_path,"EUR",cur_iter,".eid")))
  .<-capture.output(file.remove(paste0(sim_GWAS_path,"pheno_EUR",cur_iter,".txt")))
  .<-capture.output(file.remove(paste0(sim_GWAS_path,"EUR_SNP",cur_iter,".txt")))
  .<-capture.output(file.remove(paste0(sim_GWAS_path,"ukbb_pop_all",cur_iter,".bed")))
  .<-capture.output(file.remove(paste0(sim_GWAS_path,"ukbb_pop_all",cur_iter,".fam")))
  .<-capture.output(file.remove(paste0(sim_GWAS_path,"ukbb_pop_all",cur_iter,".bim")))
  .<-capture.output(file.remove(paste0(sim_GWAS_path,"GWAS_T2D_FTO_",cur_iter,".T2D.glm.logistic")))
  .<-capture.output(file.remove(paste0(sim_GWAS_path,"GWAS_T2D_FTO_",cur_iter,".T2D.glm.logistic.id")))
  .<-capture.output(file.remove(paste0(sim_GWAS_path,"GWAS_T2D_FTO_",cur_iter,".log")))
  
  ## COJO
  
  ma_file<- data.frame(SNP=gwas_res$ID,
    A1=bim_ref$ref,A2=bim_ref$alt,
    freq=EAF[filter_SNP_vec],b=gwas_res$BETA,
    se=gwas_res$SE,p=gwas_res$P,
    N=nrow(ref))
  
  bim_sample<-bim_ref
  colnames(ref_sample)<-bim_sample$id
  
  .<-capture.output(write_plink(paste0(foldpath,"ukbb_pop_",cur_iter),as.matrix(t(ref_sample)),
    bim =  bim_sample))
  
  readr::write_tsv(ma_file,paste0(foldpath,"ukbb_pop_",cur_iter,".ma"))
  aaa<-system(paste0("gcta64  --bfile ",foldpath,"ukbb_pop_",cur_iter," --cojo-file ",foldpath,"ukbb_pop_", cur_iter,".ma --cojo-slct --cojo-p 1e-5 --out ",foldpath,"test_",cur_iter),intern = TRUE)
  
  cojo_jma<-read.table(paste0(foldpath,"test_",cur_iter,".jma.cojo"),header = T)    
  
  aft_GWAS_jma<-cojo_jma$SNP[cojo_jma$pJ<0.05/N_SNP]
  #length(aft_GWAS_jma)
  #aft_GWAS_jma
  #gwas_res$ID[Nonnull_index]
  predict_nonnull_index<-which(gwas_res$ID%in%aft_GWAS_jma)
  #filter_SNP_vec[predict_nonnull_index]
  #Nonnull_index_filter_less
  #predict_nonnull_index
  #EAF[predict_nonnull_index]
  #EAF[Nonnull_index]
  cojo_pos<-predict_nonnull_index
  
  print(paste0("COJOJMA: len:",length(cojo_pos),", true select:",sum(cojo_pos%in%Nonnull_index_filter_less)))
  

  .<-capture.output(file.remove(paste0(foldpath,"test_",cur_iter,".cma.cojo")))
  .<-capture.output(file.remove(paste0(foldpath,"test_",cur_iter,".jma.cojo")))
  .<-capture.output(file.remove(paste0(foldpath,"test_",cur_iter,".ldr.cojo")))
  .<-capture.output(file.remove(paste0(foldpath,"test_",cur_iter,".log")))
  .<-capture.output(file.remove(paste0(foldpath,"ukbb_pop_",cur_iter,".ma")))
  .<-capture.output(file.remove(paste0(foldpath,"ukbb_pop_",cur_iter,".bim")))
  .<-capture.output(file.remove(paste0(foldpath,"ukbb_pop_",cur_iter,".bed")))
  .<-capture.output(file.remove(paste0(foldpath,"ukbb_pop_",cur_iter,".fam")))
  
  
  
  ########################################################
  #######  study_info_scaled
  ########################################################
  
  study_info<-lapply(1:nrow(gwas_res), function(i){
    study.m = list(Coeff=gwas_res$BETA[i],Covariance=gwas_res$SE[i]^2,Sample_size=gwas_res$OBS_CT[i])
    study.m
  })
  study_info_scaled<-study_info
  
  for(i in 1:length(study_info)){
    cur_p<-EAF[i]
    a<-sqrt(2*cur_p*(1-cur_p))
    study_info_scaled[[i]]$Coeff<-study_info_scaled[[i]]$Coeff*a
    study_info_scaled[[i]]$Covariance<-study_info_scaled[[i]]$Covariance*a^2
  }
  
  SNP_names<-bim_sample$id
  colnames(ref_sample)<-paste0('SNP',1:length(SNP_names))
  pheno_covariates_sp<-Phenotype[index1]
  
  UKBB_pop_all<-data.frame(Phenotype=scale(pheno_covariates_sp,scale = F),scale(ref_sample,scale = F))
  # The first column is phenotype:"Y"
  colnames(UKBB_pop_all)[1]<-"Y"
  # V is the intercept term
  #UKBB_pop_all$V<-1
  #UKBB_pop_all<-UKBB_pop_all[,c(1,ncol(UKBB_pop_all),2:(ncol(UKBB_pop_all)-1))]
  UKBB_pop_all<-as.matrix(UKBB_pop_all)
  
  if(is.null(colnames(UKBB_pop_all)[1])){
    stop("The column name of UKBB matrix should be clear.(Which one is SNP/GPC/risk factor")
  }
  
  var_SNP<-paste0("SNP",1:(N_SNP))
  len_SNP<-length(var_SNP)
  var_GPC<-NULL
  len_GPC<-length(var_GPC)
  var_nonSNP<-c(var_GPC)
  len_nonSNP<-length(var_nonSNP)
  var_names_full_fit <- c(var_SNP,var_nonSNP)
  colname_UKBB<-colnames(UKBB_pop_all)
  len_U1 = len_SNP+len_nonSNP
  len_U2 = len_SNP
  len_beta = len_U1
  len_theta = len_SNP+len_GPC
  len_U = len_U1 + len_U2
  N_Pop<-nrow(UKBB_pop_all)
  
  if(0){
  xtilde<-t(UKBB_pop_all[,-1])%*%UKBB_pop_all[,-1]
  ytilde<-sapply(1:N_SNP, function(snp_id){
    u2_id<-UKBB_pop_all[,c(paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
    c(c(u2_id)%*%UKBB_pop_all[,paste0("SNP",snp_id)])
  })
  
  C_half<-diag(sapply(1:N_SNP,function(i){
    1/sqrt(study_info[[i]]$Covariance)/c(UKBB_pop_all[,c(paste0("SNP",i))]%*%UKBB_pop_all[,c(paste0("SNP",i))])
  }))
  
  xtilde1<-C_half%*%xtilde
  ytilde1<-c(C_half%*%ytilde)
  #xtilde1<-xtilde
  #ytilde1<-c(ytilde)
  
  aa<-cv.glmnet(xtilde1,ytilde1,intercept=F,standardize=F)
  
  
  lambda_initial<-aa$lambda.min*nrow(xtilde1)
  beta_initial = as.numeric(coef(aa, s = "lambda.min",exact =T))[-1]
  #which(beta_initial!=0)
  #beta_initial[which(beta_initial!=0)]
  #sigma_est<-estimateSigma(x= (xtilde1),y= (ytilde1),intercept=F,standardize = F)
  #out = fixedLassoInf(x=(pseudo_X),y= (pseudo_y),beta_initial,lambda_initial,sigma = sigma_est$sigmahat)
  #out = fixedLassoInf(x= (xtilde1),y= (ytilde1),beta_initial,lambda_initial,intercept=F,tol.beta = 1e-4,sigma = sigma_est$sigmahat)
  #lasgw_pos<-out$vars[which(out$pv<0.05/length(out$vars))]
  
  y_hat<-c(xtilde1%*%(beta_initial))
  sig<-sum((ytilde1-y_hat)^2)/(nrow(UKBB_pop_all) - sum(beta_initial!=0)-1)
  xtxinv<-diag(solve(t(xtilde1)%*%xtilde1))*sig
  
  aa_final<-1-pchisq(beta_initial^2/xtxinv,1)
  xx_final<-which(aa_final<.05/sum(beta_initial!=0))
  lasgw_pos<-xx_final
  }
  ############## Get initial value 
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
  C_half<-diag(1,len_U)
  UKBB_pop<-UKBB_pop_all
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  
  pseudo_X<-pseudo_Xy_list$pseudo_X
  pseudo_y<-pseudo_Xy_list$pseudo_y
  
  initial_sf<-N_Pop^2/nrow(pseudo_X)
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  library(selectiveInference)
  UKBB_full_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  #lambda_initial<-UKBB_full_fit$lambda.min*nrow(pseudo_X)
  beta_initial = as.numeric(coef(UKBB_full_fit, s = "lambda.min"))[-1]
  ############## Get C
  
  beta<-beta_initial
  N_Pop<-nrow(UKBB_pop)
  var_11_half<-UKBB_pop[,-1]*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11)
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop_all[,c(paste0("SNP",i))]%*%UKBB_pop_all[,c(paste0("SNP",i))]))^2
  }))
  C_22_half<-diag(sqrt(diag(C_22)))
  library(magic,quietly = T)
  C_half<-adiag(C_11_half,C_22_half)
  
  
  ############## Beta

  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop^2/nrow(pseudo_X)
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  UKBB_full_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_initial<-UKBB_full_fit$lambda.min*nrow(pseudo_X)
  beta_initial = as.numeric(coef(UKBB_full_fit, s = "lambda.min",exact = T))[-1]
  which(beta_initial!=0)
  #sigma_est<-estimateSigma(x= (pseudo_X),y= (pseudo_y),intercept=F,standardize = F)
  #sigma_est<-estimateSigma(x= (UKBB_pop[,-1]),y= (UKBB_pop[,1]),intercept=F,standardize = T)
  #out = fixedLassoInf(x= (pseudo_X),y= (pseudo_y),beta_initial,lambda_initial,intercept=F,tol.beta = 1e-4,sigma = sigma_est$sigmahat)
  #out$vars[which(out$pv< 0.05/length(out$vars))]
  
  beta<-beta_initial
  
  ######## Inference 
  study_info<-study_info_scaled
  beta<-beta_initial
  N_Pop<-nrow(UKBB_pop)
  var_11_half<-UKBB_pop[,-1]*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11)
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22_half<-diag(sqrt(diag(C_22)))
  library(magic,quietly = T)
  C_half<-adiag(C_11_half,C_22_half)
  
  
  
  xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
  M_inv_half<-cbind(xtx,xtx)%*%C_half
  M_inv<-M_inv_half%*%t(M_inv_half)*N_Pop
  M<-solve(M_inv)
  ###### When computing C related items, please always use C_half 
  hat_y<-UKBB_pop[,-1]%*%beta
  ytilde<-sapply(1:N_SNP, function(snp_id){
    u2_id<-UKBB_pop[,c(paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
    c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
  })
  D_mat<-sapply(1:N_SNP, function(snp_id){
    UKBB_pop[,c(paste0("SNP",snp_id))]%*%UKBB_pop[,paste0("SNP",snp_id)]
  })
  theta_cov<-sapply(1:N_SNP, function(snp_id){
    study_info[[snp_id]]$Covariance
  })
  a<-xtx%*%C_11%*%t(UKBB_pop[,-1])%*%(hat_y -UKBB_pop[,1])/N_Pop+
    xtx%*%C_22%*%(xtx%*%beta - ytilde/N_Pop)
  
  debias_part<-M%*%a
  
  beta_star<-beta+debias_part
  beta_star[abs(beta_star)<1e-5]<-0
  #beta<-beta_star
  hat_y2<-UKBB_pop[,-1]%*%beta_star
  est_sigma<-sum((UKBB_pop[,1]-hat_y)^2)/(N_Pop - len_SNP)
  variance_1<-M%*%xtx%*%C_11%*%xtx%*%C_11%*%xtx%*%t(M)*est_sigma
  
  variance_2<-M%*%xtx%*%C_22%*%(diag(D_mat)/N_Pop)%*%diag(theta_cov)%*%diag(D_mat)%*%C_22%*%xtx%*%t(M)
  final_v<-diag(variance_1)+diag(variance_2)
  
  which.min(final_v)
  aa_final<-1-pchisq(beta^2/final_v,1)
  xx_final<-which(aa_final<0.05/length(aa_final))
  xx_final
  print(paste0("Only GWAS: len:",length(lasgw_pos),", true select:",sum(lasgw_pos%in%Nonnull_index_filter_less)))
  
  
  UKBB_full_fit2 <-cv.glmnet(x=UKBB_pop_all[,var_names_full_fit],y=UKBB_pop_all[,"Y"])
  lambda_initial2<-UKBB_full_fit2$lambda.min*nrow(UKBB_pop_all)
  beta_initial2 = as.numeric(coef(UKBB_full_fit2, s = "lambda.min"))[-1]
  #which(beta_initial2!=0)
  out2 = fixedLassoInf(x=UKBB_pop_all[,var_names_full_fit],y=UKBB_pop_all[,"Y"],beta_initial2,lambda_initial2)
  print("individual level data prediction")
  lasso_pos<-out2$vars[which(out2$pv<0.05/length(out2$vars))]
  
  print(paste0("Only GWAS: len:",length(lasso_pos),", true select:",sum(lasso_pos%in%Nonnull_index_filter_less)))
  
  