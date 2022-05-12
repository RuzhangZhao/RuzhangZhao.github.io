if(0){
  UKBB_pop<-UKBB_pop_all
  study_info<-study_info_scaled 
  ld_cut = 0.9
  cor_cut=0.9
  filter_index=TRUE
  filter_by_EAF=FALSE
  cv_with_individual=FALSE
  C_update = "Method1"
  p_val_cut = "Method1"
  kfolds=10
  kfolds=10
  BMI<-bmi_sp
}
aa<-function(section_num,allp=FALSE){
  fold_path<-paste0("run_batch",section_num,"/result/")
  file_list<-list.files(fold_path)
  final_list<-list()
  count_final_list<-1
  for (i in file_list ){
    final_list1<-readRDS(paste0(fold_path,i)) 
    for (i in 1:length(final_list1)){
      if(!is.null(final_list1[[i]])){
        final_list[[count_final_list]]=final_list1[[i]]
        count_final_list = count_final_list+1
      }
    }
  }
  true_pos = final_list[[length(final_list)]]$true_pos
  print(paste0("true:", paste0(final_list[[length(final_list)]]$true_pos,collapse = " ") ))
  count = 0
  count_list<-c()
  for (i in 1:length(final_list)){
    
    if(!is.null(final_list[[i]])){
      count_list<-c(count_list,i)
      count = count+1
    }
  }
  print(count)
  fdr_individual<-c()
  fdr_individual2<-c()
  fdr_gwas<-c()
  fdr_gwas2<-c()
  fdr_final<-c()
  fdr_final2<-c()
  #fdr_initial<-c()
  
  tpr_individual<-c()
  tpr_individual2<-c()
  tpr_gwas<-c()
  tpr_gwas2<-c()
  tpr_final<-c()
  tpr_final2<-c()
  
  fdr_individual<-c()
  fdr_individual2<-c()
  fdr_gwas<-c()
  fdr_gwas2<-c()
  fdr_final<-c()
  fdr_final2<-c()
  #fdr_initial<-c()
  
  tpr_individual<-c()
  tpr_individual2<-c()
  tpr_gwas<-c()
  tpr_gwas2<-c()
  tpr_final<-c()
  tpr_final2<-c()
  
  ## type I error 
  for (i in 1:length(final_list)){
    if(!is.null(final_list[[i]])){
      
      fdr_gwas<-c(fdr_gwas,sum(final_list[[i]]$cojo%in% true_pos))
      fdr_final<-c(fdr_final,sum(final_list[[i]]$proposed%in% true_pos))
      fdr_final2<-c(fdr_final2,sum(final_list[[i]]$susiE%in% true_pos))
      tpr_gwas<-c(tpr_gwas,length(final_list[[i]]$cojo))
      tpr_final<-c(tpr_final,sum(!is.na(final_list[[i]]$proposed)))
      tpr_final2<-c(tpr_final2,sum(!is.na(final_list[[i]]$susiE)))
    }
  }
  x1<-(sum(tpr_gwas) - sum(fdr_gwas))/count/254
  y1<-(sum(tpr_final) - sum(fdr_final))/count/254
  z1<-(sum(tpr_final2) - sum(fdr_final2))/count/254
  
  
  ## FDR 
  
  fdr_individual<-c()
  fdr_individual2<-c()
  fdr_gwas<-c()
  fdr_gwas2<-c()
  fdr_final<-c()
  fdr_final2<-c()
  #fdr_initial<-c()
  
  tpr_individual<-c()
  tpr_individual2<-c()
  tpr_gwas<-c()
  tpr_gwas2<-c()
  tpr_final<-c()
  tpr_final2<-c()
  
  for (i in 1:length(final_list)){
    if(!is.null(final_list[[i]])){
      
      fdr_gwas<-c(fdr_gwas,sum(final_list[[i]]$cojo%in% true_pos)/length(final_list[[i]]$cojo))
      fdr_final<-c(fdr_final,sum(final_list[[i]]$proposed%in% true_pos)/max(length(final_list[[i]]$proposed),1))
      fdr_final2<-c(fdr_final2,sum(final_list[[i]]$susiE%in% true_pos)/max(1,length(final_list[[i]]$susiE)))
      
      
      tpr_gwas<-c(tpr_gwas,sum(final_list[[i]]$cojo%in% true_pos)/length(true_pos))
      tpr_final<-c(tpr_final,sum(final_list[[i]]$proposed%in% true_pos)/length(true_pos))
      tpr_final2<-c(tpr_final2,sum(final_list[[i]]$susiE%in% true_pos)/length(true_pos))
      
      
    }
  }
  
  x2<-(1-mean(fdr_gwas))
  #y2<-(1-mean(fdr_final[which(fdr_final!=0)]))
  y2<-(1-mean(fdr_final,na.rm=T))
  x3<-(mean(tpr_gwas))
  #y3<-(mean(tpr_final[which(tpr_final!=0)]))
  y3<-(mean(tpr_final))
  z2<-(1-mean(fdr_final2))
  z3<-(mean(tpr_final2))
  print(paste0("COJO:",round(x1,8),",",round(x2,8),",",round(x3,8)))
  print(paste0("Proposed_block:",round(y1,8),",",round(y2,8),",",round(y3,8)))
  print(paste0("susiE",round(z1,8),",",round(z2,8),",",round(z3,8)))
  if(allp){
    
    type_final<-c()
    typetpr_final<-c()
    fdr_final<-c()
    tpr_final<-c()
    for (i in 1:length(final_list)){
      if(!is.null(final_list[[i]])){
        type_final<-c(type_final,sum(final_list[[i]]$proposed2%in% true_pos))
        typetpr_final<-c(typetpr_final,sum(!is.na(final_list[[i]]$proposed2)))
        fdr_final<-c(fdr_final,sum(final_list[[i]]$proposed2%in% true_pos)/max(length(final_list[[i]]$proposed2),1) )
        tpr_final<-c(tpr_final,sum(final_list[[i]]$proposed2%in% true_pos)/length(true_pos))
      }
    }
    type1<-(sum(typetpr_final)-sum(type_final))/count/254
    fdr<-(1-mean(fdr_final))
    power<-(mean(tpr_final))
    print(paste0("Proposed_C:",round(type1,8),",",round(fdr,8),",",round(power,8)))
    type_final<-c()
    typetpr_final<-c()
    fdr_final<-c()
    tpr_final<-c()
    for (i in 1:length(final_list)){
      if(!is.null(final_list[[i]])){
        type_final<-c(type_final,sum(final_list[[i]]$prop_pos_nolasso%in% true_pos))
        typetpr_final<-c(typetpr_final,sum(!is.na(final_list[[i]]$prop_pos_nolasso)))
        fdr_final<-c(fdr_final,sum(final_list[[i]]$prop_pos_nolasso%in% true_pos)/max(length(final_list[[i]]$prop_pos_nolasso),1))
        tpr_final<-c(tpr_final,sum(final_list[[i]]$prop_pos_nolasso%in% true_pos)/length(true_pos))
      }
    }
    type1<-(sum(typetpr_final)-sum(type_final) )/count/254
    fdr<-(1-mean(fdr_final))
    power<-(mean(tpr_final))
    print(paste0("Proposed_block_nolasso:",round(type1,8),",",round(fdr,8),",",round(power,8)))
    type_final<-c()
    typetpr_final<-c()
    fdr_final<-c()
    tpr_final<-c()
    for (i in 1:length(final_list)){
      if(!is.null(final_list[[i]])){
        type_final<-c(type_final,sum(final_list[[i]]$prop_pos2_nolasso%in% true_pos))
        typetpr_final<-c(typetpr_final,sum(!is.na(final_list[[i]]$prop_pos2_nolasso)))
        fdr_final<-c(fdr_final,sum(final_list[[i]]$prop_pos2_nolasso%in% true_pos)/max(length(final_list[[i]]$prop_pos2_nolasso),1))
        tpr_final<-c(tpr_final,sum(final_list[[i]]$prop_pos2_nolasso%in% true_pos)/length(true_pos))
      }
    }
    type1<-(sum(typetpr_final)-sum(type_final))/count/254
    fdr<-(1-mean(fdr_final))
    power<-(mean(tpr_final))
    print(paste0("Proposed_C_nolasso:",round(type1,8),",",round(fdr,8),",",round(power,8)))
  }
}
bb2<-function(section_num,allp=FALSE){
  fold_path<-paste0("run_batch",section_num,"/result/")
  file_list<-list.files(fold_path)
  final_list<-list()
  count_final_list<-1
  for (i in file_list ){
    final_list1<-readRDS(paste0(fold_path,i)) 
    for (i in 1:length(final_list1)){
      if(!is.null(final_list1[[i]])){
        final_list[[count_final_list]]=final_list1[[i]]
        count_final_list = count_final_list+1
      }
    }
  }
  true_pos = final_list[[length(final_list)]]$true_pos
  print(paste0("true:", paste0(final_list[[length(final_list)]]$true_pos,collapse = " ") ))
  count = 0
  count_list<-c()
  for (i in 1:length(final_list)){
    
    if(!is.null(final_list[[i]])){
      count_list<-c(count_list,i)
      count = count+1
    }
  }
  print(count)
  fdr_individual<-c()
  fdr_individual2<-c()
  fdr_gwas<-c()
  fdr_gwas2<-c()
  fdr_final<-c()
  fdr_final2<-c()
  #fdr_initial<-c()
  
  tpr_individual<-c()
  tpr_individual2<-c()
  tpr_gwas<-c()
  tpr_gwas2<-c()
  tpr_final<-c()
  tpr_final2<-c()
  
  fdr_individual<-c()
  fdr_individual2<-c()
  fdr_gwas<-c()
  fdr_gwas2<-c()
  fdr_final<-c()
  fdr_final2<-c()
  #fdr_initial<-c()
  
  tpr_individual<-c()
  tpr_individual2<-c()
  tpr_gwas<-c()
  tpr_gwas2<-c()
  tpr_final<-c()
  tpr_final2<-c()
  ## type I error 
  for (i in 1:length(final_list)){
    if(!is.null(final_list[[i]])){
      
      fdr_gwas<-c(fdr_gwas,sum(final_list[[i]]$cojo%in% true_pos))
      fdr_final<-c(fdr_final,sum(final_list[[i]]$proposed%in% true_pos))
      fdr_final2<-c(fdr_final2,sum(final_list[[i]]$new%in% true_pos))
      tpr_gwas<-c(tpr_gwas,length(final_list[[i]]$cojo))
      tpr_final<-c(tpr_final,sum(!is.na(final_list[[i]]$proposed)))
      tpr_final2<-c(tpr_final2,sum(!is.na(final_list[[i]]$new)))
    }
  }
  x1<-(sum(tpr_gwas) - sum(fdr_gwas))/count/254
  y1<-(sum(tpr_final) - sum(fdr_final))/count/254
  z1<-(sum(tpr_final2) - sum(fdr_final2))/count/254
  
  
  ## FDR 
  
  fdr_individual<-c()
  fdr_individual2<-c()
  fdr_gwas<-c()
  fdr_gwas2<-c()
  fdr_final<-c()
  fdr_final2<-c()
  #fdr_initial<-c()
  
  tpr_individual<-c()
  tpr_individual2<-c()
  tpr_gwas<-c()
  tpr_gwas2<-c()
  tpr_final<-c()
  tpr_final2<-c()
  tpr_gwas11<-c()
  tpr_gwas22<-c()
  tpr_gwas33<-c()
  tpr_final11<-c()
  tpr_final22<-c()
  tpr_final33<-c()
  for (i in 1:length(final_list)){
    if(!is.null(final_list[[i]])){
      
      fdr_gwas<-c(fdr_gwas,sum(final_list[[i]]$cojo%in% true_pos)/length(final_list[[i]]$cojo))
      fdr_final<-c(fdr_final,sum(final_list[[i]]$proposed%in% true_pos)/max(length(final_list[[i]]$proposed),1))
      fdr_final2<-c(fdr_final2,sum(final_list[[i]]$new%in% true_pos)/max(1,length(final_list[[i]]$new)))
      
      tpr_gwas11<-c(tpr_gwas11,sum(final_list[[i]]$cojo%in% true_pos[1]))
      tpr_gwas22<-c(tpr_gwas22,sum(final_list[[i]]$cojo%in% true_pos[2]))
      tpr_gwas33<-c(tpr_gwas33,sum(final_list[[i]]$cojo%in% true_pos[3]))
      tpr_final11<-c(tpr_final11,sum(final_list[[i]]$proposed2%in% true_pos[1]))
      tpr_final22<-c(tpr_final22,sum(final_list[[i]]$proposed2%in% true_pos[2]))
      tpr_final33<-c(tpr_final33,sum(final_list[[i]]$proposed2%in% true_pos[3]))
      tpr_gwas<-c(tpr_gwas,sum(final_list[[i]]$cojo%in% true_pos)/length(true_pos))
      tpr_final<-c(tpr_final,sum(final_list[[i]]$proposed%in% true_pos)/length(true_pos))
      tpr_final2<-c(tpr_final2,sum(final_list[[i]]$new%in% true_pos)/length(true_pos))
      
      
    }
  }
  
  x2<-(1-mean(fdr_gwas))
  #y2<-(1-mean(fdr_final[which(fdr_final!=0)]))
  y2<-(1-mean(fdr_final,na.rm=T))
  x3<-(mean(tpr_gwas))
  print(paste0(mean(tpr_gwas11),",",mean(tpr_gwas22),",",mean(tpr_gwas33)))
  print(paste0(mean(tpr_final11),",",mean(tpr_final22),",",mean(tpr_final33)))
      print(mean(tpr_gwas11))
      #y3<-(mean(tpr_final[which(tpr_final!=0)]))
      y3<-(mean(tpr_final[which(tpr_final!=0)]))
      z2<-(1-mean(fdr_final2))
      z3<-(mean(tpr_final2))
      #print(paste0("COJO:",round(x1,8),",",round(x2,8),",",round(x3,8)))
      #print(paste0("Proposed_block:",round(y1,8),",",round(y2,8),",",round(y3,8)))
      #print(paste0(round(z1,8),",",round(z2,8),",",round(z3,8)))
     
}
bb<-function(section_num,allp=FALSE){
  fold_path<-paste0("run_batch",section_num,"/result/")
  file_list<-list.files(fold_path)
  final_list<-list()
  count_final_list<-1
  for (i in file_list ){
    final_list1<-readRDS(paste0(fold_path,i)) 
    for (i in 1:length(final_list1)){
      if(!is.null(final_list1[[i]])){
        final_list[[count_final_list]]=final_list1[[i]]
        count_final_list = count_final_list+1
      }
    }
  }
  true_pos = final_list[[length(final_list)]]$true_pos
  print(paste0("true:", paste0(final_list[[length(final_list)]]$true_pos,collapse = " ") ))
  count = 0
  count_list<-c()
  for (i in 1:length(final_list)){
    
    if(!is.null(final_list[[i]])){
      count_list<-c(count_list,i)
      count = count+1
    }
  }
  print(count)
  fdr_individual<-c()
  fdr_individual2<-c()
  fdr_gwas<-c()
  fdr_gwas2<-c()
  fdr_final<-c()
  fdr_final2<-c()
  #fdr_initial<-c()
  
  tpr_individual<-c()
  tpr_individual2<-c()
  tpr_gwas<-c()
  tpr_gwas2<-c()
  tpr_final<-c()
  tpr_final2<-c()
  
  fdr_individual<-c()
  fdr_individual2<-c()
  fdr_gwas<-c()
  fdr_gwas2<-c()
  fdr_final<-c()
  fdr_final2<-c()
  #fdr_initial<-c()
  
  tpr_individual<-c()
  tpr_individual2<-c()
  tpr_gwas<-c()
  tpr_gwas2<-c()
  tpr_final<-c()
  tpr_final2<-c()
  
  ## type I error 
  for (i in 1:length(final_list)){
    if(!is.null(final_list[[i]])){
      
      fdr_gwas<-c(fdr_gwas,sum(final_list[[i]]$cojo%in% true_pos))
      fdr_final<-c(fdr_final,sum(final_list[[i]]$proposed%in% true_pos))
      fdr_final2<-c(fdr_final2,sum(final_list[[i]]$new%in% true_pos))
      tpr_gwas<-c(tpr_gwas,length(final_list[[i]]$cojo))
      tpr_final<-c(tpr_final,sum(!is.na(final_list[[i]]$proposed)))
      tpr_final2<-c(tpr_final2,sum(!is.na(final_list[[i]]$new)))
    }
  }
  x1<-(sum(tpr_gwas) - sum(fdr_gwas))/count/254
  y1<-(sum(tpr_final) - sum(fdr_final))/count/254
  z1<-(sum(tpr_final2) - sum(fdr_final2))/count/254
  
  
  ## FDR 
  
  fdr_individual<-c()
  fdr_individual2<-c()
  fdr_gwas<-c()
  fdr_gwas2<-c()
  fdr_final<-c()
  fdr_final2<-c()
  #fdr_initial<-c()
  
  tpr_individual<-c()
  tpr_individual2<-c()
  tpr_gwas<-c()
  tpr_gwas2<-c()
  tpr_final<-c()
  tpr_final2<-c()
  
  for (i in 1:length(final_list)){
    if(!is.null(final_list[[i]])){
      
      fdr_gwas<-c(fdr_gwas,sum(final_list[[i]]$cojo%in% true_pos)/length(final_list[[i]]$cojo))
      fdr_final<-c(fdr_final,sum(final_list[[i]]$proposed%in% true_pos)/max(length(final_list[[i]]$proposed),1))
      fdr_final2<-c(fdr_final2,sum(final_list[[i]]$new%in% true_pos)/max(1,length(final_list[[i]]$new)))
      
      
      tpr_gwas<-c(tpr_gwas,sum(final_list[[i]]$cojo%in% true_pos)/length(true_pos))
      tpr_final<-c(tpr_final,sum(final_list[[i]]$proposed%in% true_pos)/length(true_pos))
      tpr_final2<-c(tpr_final2,sum(final_list[[i]]$new%in% true_pos)/length(true_pos))
      
      
    }
  }
  
  x2<-(1-mean(fdr_gwas))
  #y2<-(1-mean(fdr_final[which(fdr_final!=0)]))
  y2<-(1-mean(fdr_final,na.rm=T))
  x3<-(mean(tpr_gwas))
  #y3<-(mean(tpr_final[which(tpr_final!=0)]))
  y3<-(mean(tpr_final))
  z2<-(1-mean(fdr_final2))
  z3<-(mean(tpr_final2))
  print(paste0("COJO:",round(x1,8),",",round(x2,8),",",round(x3,8)))
  print(paste0("Proposed_block:",round(y1,8),",",round(y2,8),",",round(y3,8)))
  #print(paste0(round(z1,8),",",round(z2,8),",",round(z3,8)))
  if(allp){
    
    type_final<-c()
    typetpr_final<-c()
    fdr_final<-c()
    tpr_final<-c()
    for (i in 1:length(final_list)){
      if(!is.null(final_list[[i]])){
        type_final<-c(type_final,sum(final_list[[i]]$proposed2%in% true_pos))
        typetpr_final<-c(typetpr_final,sum(!is.na(final_list[[i]]$proposed2)))
        fdr_final<-c(fdr_final,sum(final_list[[i]]$proposed2%in% true_pos)/max(length(final_list[[i]]$proposed2),1) )
        tpr_final<-c(tpr_final,sum(final_list[[i]]$proposed2%in% true_pos)/length(true_pos))
      }
    }
    type1<-(sum(typetpr_final)-sum(type_final))/count/254
    fdr<-(1-mean(fdr_final))
    power<-(mean(tpr_final))
    print(paste0("Proposed_C:",round(type1,8),",",round(fdr,8),",",round(power,8)))
    type_final<-c()
    typetpr_final<-c()
    fdr_final<-c()
    tpr_final<-c()
    for (i in 1:length(final_list)){
      if(!is.null(final_list[[i]])){
        type_final<-c(type_final,sum(final_list[[i]]$prop_pos_nolasso%in% true_pos))
        typetpr_final<-c(typetpr_final,sum(!is.na(final_list[[i]]$prop_pos_nolasso)))
        fdr_final<-c(fdr_final,sum(final_list[[i]]$prop_pos_nolasso%in% true_pos)/max(length(final_list[[i]]$prop_pos_nolasso),1))
        tpr_final<-c(tpr_final,sum(final_list[[i]]$prop_pos_nolasso%in% true_pos)/length(true_pos))
      }
    }
    type1<-(sum(typetpr_final)-sum(type_final) )/count/254
    fdr<-(1-mean(fdr_final))
    power<-(mean(tpr_final))
    print(paste0("Proposed_block_nolasso:",round(type1,8),",",round(fdr,8),",",round(power,8)))
    type_final<-c()
    typetpr_final<-c()
    fdr_final<-c()
    tpr_final<-c()
    for (i in 1:length(final_list)){
      if(!is.null(final_list[[i]])){
        type_final<-c(type_final,sum(final_list[[i]]$prop_pos2_nolasso%in% true_pos))
        typetpr_final<-c(typetpr_final,sum(!is.na(final_list[[i]]$prop_pos2_nolasso)))
        fdr_final<-c(fdr_final,sum(final_list[[i]]$prop_pos2_nolasso%in% true_pos)/max(length(final_list[[i]]$prop_pos2_nolasso),1))
        tpr_final<-c(tpr_final,sum(final_list[[i]]$prop_pos2_nolasso%in% true_pos)/length(true_pos))
      }
    }
    type1<-(sum(typetpr_final)-sum(type_final))/count/254
    fdr<-(1-mean(fdr_final))
    power<-(mean(tpr_final))
    print(paste0("Proposed_C_nolasso:",round(type1,8),",",round(fdr,8),",",round(power,8)))
  }
}
cc<-function(section_num,allp=FALSE){
  fold_path<-paste0("run_batch",section_num,"/result/")
  file_list<-list.files(fold_path)
  final_list<-list()
  count_final_list<-1
  for (i in file_list ){
    final_list1<-readRDS(paste0(fold_path,i)) 
    for (i in 1:length(final_list1)){
      if(!is.null(final_list1[[i]])){
        final_list[[count_final_list]]=final_list1[[i]]
        count_final_list = count_final_list+1
      }
    }
  }
  true_pos = final_list[[length(final_list)]]$true_pos
  print(paste0("true:", paste0(final_list[[length(final_list)]]$true_pos,collapse = " ") ))
  count = 0
  count_list<-c()
  for (i in 1:length(final_list)){
    
    if(!is.null(final_list[[i]])){
      count_list<-c(count_list,i)
      count = count+1
    }
  }
  print(count)
  fdr_individual<-c()
  fdr_individual2<-c()
  fdr_gwas<-c()
  fdr_gwas2<-c()
  fdr_final<-c()
  fdr_final2<-c()
  #fdr_initial<-c()
  
  tpr_individual<-c()
  tpr_individual2<-c()
  tpr_gwas<-c()
  tpr_gwas2<-c()
  tpr_final<-c()
  tpr_final2<-c()
  
  fdr_individual<-c()
  fdr_individual2<-c()
  fdr_gwas<-c()
  fdr_gwas2<-c()
  fdr_final<-c()
  fdr_final2<-c()
  #fdr_initial<-c()
  
  tpr_individual<-c()
  tpr_individual2<-c()
  tpr_gwas<-c()
  tpr_gwas2<-c()
  tpr_final<-c()
  tpr_final2<-c()
  
  ## type I error 
  for (i in 1:length(final_list)){
    if(!is.null(final_list[[i]])){
      
      fdr_gwas<-c(fdr_gwas,sum(final_list[[i]]$cojo%in% true_pos))
      fdr_final<-c(fdr_final,sum(final_list[[i]]$proposed%in% true_pos))
      fdr_final2<-c(fdr_final2,sum(final_list[[i]]$individual%in% true_pos))
      tpr_gwas<-c(tpr_gwas,length(final_list[[i]]$cojo))
      tpr_final<-c(tpr_final,sum(!is.na(final_list[[i]]$proposed)))
      tpr_final2<-c(tpr_final2,sum(!is.na(final_list[[i]]$individual)))
    }
  }
  x1<-(sum(tpr_gwas) - sum(fdr_gwas))/count/254
  y1<-(sum(tpr_final) - sum(fdr_final))/count/254
  z1<-(sum(tpr_final2) - sum(fdr_final2))/count/254
  
  
  ## FDR 
  
  fdr_individual<-c()
  fdr_individual2<-c()
  fdr_gwas<-c()
  fdr_gwas2<-c()
  fdr_final<-c()
  fdr_final2<-c()
  #fdr_initial<-c()
  
  tpr_individual<-c()
  tpr_individual2<-c()
  tpr_gwas<-c()
  tpr_gwas2<-c()
  tpr_final<-c()
  tpr_final2<-c()
  
  for (i in 1:length(final_list)){
    if(!is.null(final_list[[i]])){
      
      fdr_gwas<-c(fdr_gwas,sum(final_list[[i]]$cojo%in% true_pos)/length(final_list[[i]]$cojo))
      fdr_final<-c(fdr_final,sum(final_list[[i]]$proposed%in% true_pos)/max(length(final_list[[i]]$proposed),1))
      fdr_final2<-c(fdr_final2,sum(final_list[[i]]$individual%in% true_pos)/max(1,length(final_list[[i]]$individual)))
      
      
      tpr_gwas<-c(tpr_gwas,sum(final_list[[i]]$cojo%in% true_pos)/length(true_pos))
      tpr_final<-c(tpr_final,sum(final_list[[i]]$proposed%in% true_pos)/length(true_pos))
      tpr_final2<-c(tpr_final2,sum(final_list[[i]]$individual%in% true_pos)/length(true_pos))
      
      
    }
  }
  
  x2<-(1-mean(fdr_gwas))
  #y2<-(1-mean(fdr_final[which(fdr_final!=0)]))
  y2<-(1-mean(fdr_final,na.rm=T))
  x3<-(mean(tpr_gwas))
  #y3<-(mean(tpr_final[which(tpr_final!=0)]))
  y3<-(mean(tpr_final[which(tpr_final!=0)]))
  z2<-(1-mean(fdr_final2))
  z3<-(mean(tpr_final2))
  print(paste0("COJO:",round(x1,8),",",round(x2,8),",",round(x3,8)))
  print(paste0("Proposed:",round(y1,8),",",round(y2,8),",",round(y3,8)))
  print(paste0(round(z1,8),",",round(z2,8),",",round(z3,8)))
  
}


## UKBB_pop should be scaled with SNPs, 
## but not with phenotype.
## Both SNP and phenotype should be centered. 
library(expm,quietly = T)
library(magic,quietly = T)
library(glmnet)
adaptiveGMMlasso<-function(UKBB_pop,study_info,
  ld_cut = 0.9,
  cor_cut=0.9,
  filter_index=TRUE,
  filter_by_EAF=FALSE,
  cv_with_individual=FALSE,
  C_update = "Method1",
  p_val_cut = "Method1",
  kfolds=10){
  UKBB_cor<-cor(UKBB_pop[,-1])
  diag(UKBB_cor) = 0
  if(filter_index){
    ### Trick 1 index filtering 
    index_filter<-c()
    pval_list<-c()
    for(i in 1:length(study_info)){
      #if(study_info[[i]]$P<p_val_cut){
      index_filter<-c(index_filter,i)
      z_here<-study_info[[i]]$Coeff^2/study_info[[i]]$Covariance
      pval_list<-c(pval_list,z_here)
      #}
    }
    index_filter0<-index_filter
    pval_list0<-pval_list
    index_filter2<-c()
    while (length(index_filter) > 1){
      UKBB_cor_i<-UKBB_cor[index_filter,index_filter]
      cur_i<-which.max(pval_list)
      index_filter2<-c(index_filter2,index_filter[cur_i])
      #print(index_filter[cur_i])
      rm_index<-c(cur_i,which(abs(UKBB_cor_i[cur_i,])>ld_cut))
      #print(index_filter[rm_index])
      index_filter<-index_filter[-rm_index]
      pval_list<-pval_list[-rm_index]
    }
    
    #Nonnull_index2<-which(index_filter2%in%Nonnull_index)
    UKBB_pop<-UKBB_pop[,c(1,(index_filter2+1))]
    study_info<-study_info[index_filter2]
    index_filter<-index_filter2
  }
  N_SNP<-ncol(UKBB_pop)-1
  colnames(UKBB_pop)[-1]<-paste0("SNP",1:(N_SNP))
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
  
  ############## Beta
  pseudo_Xy<-function(
    C_half,UKBB_pop,
    var_SNP,var_GPC,
    N_SNP,theta_UKBB_GPC,study_info){
    N_Pop<-nrow(UKBB_pop)
    pseudo_X<-C_half%*%rbind(t(UKBB_pop[,-1]),t(UKBB_pop[,var_SNP]))%*%UKBB_pop[,-1]
    pseudo_y1<-crossprod(UKBB_pop[,1],UKBB_pop[,-1])
    pseudo_y22<-sapply(1:N_SNP, function(snp_id){
      u2_id<-UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
      c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
    })
    
    pseudo_y<-c(c(c(pseudo_y1),c(pseudo_y22))%*%C_half)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
  }
  
  var_U_beta_theta_func<-function(
    UKBB_pop,
    beta,
    study_info,
    theta_UKBB_GPC=NULL){
    N_Pop<-nrow(UKBB_pop)
    xtbeta<-(UKBB_pop[,-1]%*%beta)
    var_11<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta))
    u2_theta_coef<-sapply(1:N_SNP, function(snp_id){
      if(is.null(var_GPC)){
        xtgamma_id<-c(UKBB_pop[,paste0("SNP",snp_id)]*study_info[[snp_id]]$Coeff)
      }else{
        xtgamma_id<-c((UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)))
      }
      xtbeta-xtgamma_id
    }) #col is SNP #row is sample 
    var_22<-crossprod(u2_theta_coef*UKBB_pop[,var_SNP],u2_theta_coef*UKBB_pop[,var_SNP])
    var_12<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),u2_theta_coef*UKBB_pop[,var_SNP])
    (1/N_Pop)*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
  }
  
  
  #### adaptive lasso with fast lasso computation 
  ## The initial is important
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop*var(UKBB_pop[,1])
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11) 
  
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01)#,penalty.factor = w_adaptive0)
  
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
  
  if(C_update == "Method1"){
    diag_theta_sd<-sapply(1:N_SNP,function(i){
      sqrt(study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))
    })
    cov_theta<-t(diag_theta_sd*cor(UKBB_pop[,var_SNP]))*diag_theta_sd/N_Pop
    C_22<-solve(cov_theta)
    C_22_half<-expm::sqrtm(C_22)
    C_half<-adiag(C_11_half,C_22_half)
  }else if(C_update == "Method2"){
    beta_initial<-coef(ridge_fit)[-1]
    
    var_1st_U_beta_theta<-var_U_beta_theta_func(UKBB_pop =UKBB_pop,
      beta = beta_initial,
      study_info = study_info)
    
    diag_theta_sd<-sapply(1:N_SNP,function(i){
      sqrt(study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))
    })
    cov_theta<-t(diag_theta_sd*cor(UKBB_pop[,var_SNP]))*diag_theta_sd/N_Pop
    
    var_2nd_grad_times_theta_hat = adiag(matrix(0,nrow = len_U1,ncol = len_U1),cov_theta)
    
    inv_C = var_1st_U_beta_theta + var_2nd_grad_times_theta_hat
    C_all<-solve(inv_C)
    C_half<-expm::sqrtm(C_all)

  }
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  lambda_list<-adaptivelasso_fit$lambda
  
  if(!cv_with_individual){
    beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
    index_nonzero<-which(beta!=0)
    if(length(index_nonzero) == 0){
      cv_with_individual = TRUE
    }else{
      xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
      Sigsum_half<-cbind(xtx,xtx)%*%C_half
      Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
      Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
      inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
      final_v<-diag(inv_Sigsum_scaled_nonzero)
      aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
      if(p_val_cut == "Method1"){
        candidate_pos<-index_nonzero[which(aa_final<0.05/N_SNP)]
      }else if(p_val_cut == "Method2"){
        candidate_pos<-index_nonzero[which(aa_final<0.05/length(index_nonzero))]
      }else{
        stop("wrong p_val_cut")
      }
      
      if(length(candidate_pos) == 0){
        cv_with_individual = TRUE
      }
    }
  }
  
  if(cv_with_individual){
    #index11<-which(adaptivelasso_fit$nzero>1)
    #index11<-unique(c(max(min(index11)-1,1),index11))
    #lambda_list<-lambda_list[index11]
    
    library(caret)
    index_fold<-createFolds(as.numeric(UKBB_pop[,1]>0),k = kfolds)
    a<-sapply(1:kfolds, function(cur_fold){
      index_test<-index_fold[[cur_fold]]
      UKBB_pop_train<-UKBB_pop[-index_test,]
      UKBB_pop_test<-UKBB_pop[index_test,]
      pseudo_Xy_list_train<-pseudo_Xy(C_half,UKBB_pop_train,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
      initial_sf_train<-nrow(UKBB_pop_train)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
      pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
      pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
      cv_beta<-lapply(lambda_list,function(cur_lam){
        cv_adaptivelasso_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive,lambda = cur_lam)
        coef(cv_adaptivelasso_fit)[-1]
      })
      cv_mse_sum<-sapply(1:length(cv_beta), function(kk){
        cur_beta<-cv_beta[[kk]]
        sum((UKBB_pop_test[,-1]%*%cur_beta - UKBB_pop_test[,1])^2)
      })
      cv_mse_sum
    })
    beta<-coef(adaptivelasso_fit,s=lambda_list[which.min(rowSums(a))])[-1]
    index_nonzero<-which(beta!=0)
    ########??????????????????
    if(length(index_nonzero) == 0){
      beta<-coef(adaptivelasso_fit,s=lambda_list[which(adaptivelasso_fit$nzero>0)[which.min(rowSums(a)[which(adaptivelasso_fit$nzero>0)])]])[-1]
      index_nonzero<-which(beta!=0)
    }
    
    xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
    Sigsum_half<-cbind(xtx,xtx)%*%C_half
    Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
    Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
    inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
    final_v<-diag(inv_Sigsum_scaled_nonzero)
    aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
    if(p_val_cut == "Method1"){
      candidate_pos<-index_nonzero[which(aa_final<0.05/N_SNP)]
    }else if(p_val_cut == "Method2"){
      candidate_pos<-index_nonzero[which(aa_final<0.05/length(index_nonzero))]
    }else{
      stop("wrong p_val_cut")
    }
  }
  
  if(length(candidate_pos) > 1){
    candidate_cor<-cor(UKBB_pop[,c(candidate_pos+1)])
    diag(candidate_cor)<-0
    if(max(abs(candidate_cor))>0.1){
      UKBB_pop_orig<-UKBB_pop[,c(candidate_pos+1)]
      for(i in 1:ncol(UKBB_pop_orig)){
        cur_snp<-UKBB_pop_orig[,i]
        UKBB_pop_orig[which(cur_snp == max(cur_snp)),i] = 2
        UKBB_pop_orig[which(cur_snp == min(cur_snp)),i] = 0
        UKBB_pop_orig[which(cur_snp > min(cur_snp) & cur_snp < max(cur_snp)),i] = 1
      }
      candidate_EAF<-colsums(UKBB_pop_orig)/nrow(UKBB_pop_orig)/2
      weak_among_candidate<-which(candidate_EAF < 0.05)
      if(filter_by_EAF){
        if(length(weak_among_candidate) > 1){
          gamma_adaptivelasso<-1/2
          w_adaptive_candidate<-1/(abs(beta[candidate_pos[weak_among_candidate] ]))^gamma_adaptivelasso
          tmp_count<-sapply(1:10, function(j){
            lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos[weak_among_candidate]+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive_candidate)
            abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
          })
          confident_pos_weak<-(candidate_pos[weak_among_candidate])[which(rowMeans(tmp_count) >= 0.8)]
          if(length(confident_pos_weak) == 0){
            z<-c()
            for(i in weak_among_candidate){
              z<-c(z,study_info[[candidate_pos[weak_among_candidate]]]$Coeff^2/study_info[[candidate_pos[weak_among_candidate]]]$Covariance)
            }
            confident_pos_weak<-weak_among_candidate[which.max(z)]
          }
        }else if(length(weak_among_candidate) == 1){
          confident_pos_weak<-candidate_pos[weak_among_candidate]
        }else{
          confident_pos_weak<-c()
        }
        true_weak<-c()
        if(length(weak_among_candidate) > 0){
          for( weakone in weak_among_candidate){
            if(max(abs(candidate_cor[,weakone]))<0.1){
              true_weak<-c(true_weak,weakone)      
            }
          }
        }
        if(length(true_weak)>0) true_weak<-candidate_pos[true_weak]
        confident_pos_weak <-union(confident_pos_weak,true_weak)
      }else{
        confident_pos_weak <-candidate_pos[weak_among_candidate]
      }
      gamma_adaptivelasso<-2
      w_adaptive_candidate<-1/(abs(beta[candidate_pos]))^gamma_adaptivelasso
      w_adaptive_candidate[which(candidate_pos%in%confident_pos_weak)]<-0
      tmp_count<-sapply(1:10, function(j){
        lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 0.9,penalty.factor = w_adaptive_candidate)
        abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
        #coef(lasso_fit_candidate,s='lambda.min')[-1]!=0
      })
      confident_pos<-candidate_pos[which(rowMeans(tmp_count) >= 0.8)]
      confident_pos<-union(confident_pos,confident_pos_weak)
    }else{
      print("NOT USED")
      confident_pos<-candidate_pos
    }
  }else{
    confident_pos<-candidate_pos
  }
  
  if(filter_index){
    pos<-index_filter[confident_pos]
    pos_bf<-index_filter[candidate_pos]
  }else{
    pos<-confident_pos
    pos_bf<-candidate_pos
  }
  #print(paste0("candidate",paste0(index_filter[candidate_pos],collapse = " ")))
  #print(paste0("confid",paste0(pos,collapse = " ")))
  #print(pos)
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos_bf"=pos_bf,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive,
    "index_filter"=index_filter)
}



adaptiveGMMlassoBMI<-function(UKBB_pop,BMI,study_info,
  ld_cut = 0.9,
  cor_cut=0.9,
  filter_index=TRUE,
  filter_by_EAF=FALSE,
  cv_with_individual=FALSE,
  C_update = "Method1",
  p_val_cut = "Method1",
  kfolds=10){
  UKBB_cor<-cor(UKBB_pop[,-1])
  diag(UKBB_cor) = 0
  if(filter_index){
    ### Trick 1 index filtering 
    index_filter<-c()
    pval_list<-c()
    for(i in 1:length(study_info)){
      #if(study_info[[i]]$P<p_val_cut){
      index_filter<-c(index_filter,i)
      z_here<-study_info[[i]]$Coeff^2/study_info[[i]]$Covariance
      pval_list<-c(pval_list,z_here)
      #}
    }
    index_filter0<-index_filter
    pval_list0<-pval_list
    index_filter2<-c()
    while (length(index_filter) > 1){
      UKBB_cor_i<-UKBB_cor[index_filter,index_filter]
      cur_i<-which.max(pval_list)
      index_filter2<-c(index_filter2,index_filter[cur_i])
      #print(index_filter[cur_i])
      rm_index<-c(cur_i,which(abs(UKBB_cor_i[cur_i,])>ld_cut))
      #print(index_filter[rm_index])
      index_filter<-index_filter[-rm_index]
      pval_list<-pval_list[-rm_index]
    }
    
    #Nonnull_index2<-which(index_filter2%in%Nonnull_index)
    UKBB_pop<-UKBB_pop[,c(1,(index_filter2+1))]
    study_info<-study_info[index_filter2]
    index_filter<-index_filter2
  }
  
  N_SNP<-ncol(UKBB_pop)-1
  colnames(UKBB_pop)[-1]<-paste0("SNP",1:(N_SNP))
  colnames(UKBB_pop)[1]<-"Y"
  var_SNP<-paste0("SNP",1:(N_SNP))
  len_SNP<-length(var_SNP)
  var_GPC<-NULL
  len_GPC<-length(var_GPC)
  var_nonSNP<-c(var_GPC,"BMI")
  len_nonSNP<-length(var_nonSNP)
  var_names_full_fit <- c(var_SNP,var_nonSNP)
  colname_UKBB<-colnames(UKBB_pop)
  len_U1 = len_SNP+len_nonSNP
  len_U2 = len_SNP
  len_beta = len_U1
  len_theta = len_SNP+len_GPC
  len_U = len_U1 + len_U2
  N_Pop<-nrow(UKBB_pop)
  UKBB_pop<-cbind(UKBB_pop,BMI)
  colnames(UKBB_pop)[ncol(UKBB_pop)]<-"BMI"
  ############## Beta
  pseudo_Xy<-function(
    C_half,UKBB_pop,
    var_SNP,var_GPC,
    N_SNP,theta_UKBB_GPC,study_info){
    N_Pop<-nrow(UKBB_pop)
    pseudo_X<-C_half%*%rbind(t(UKBB_pop[,-1]),t(UKBB_pop[,var_SNP]))%*%UKBB_pop[,-1]
    pseudo_y1<-crossprod(UKBB_pop[,1],UKBB_pop[,-1])
    pseudo_y22<-sapply(1:N_SNP, function(snp_id){
      u2_id<-UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
      c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
    })
    
    pseudo_y<-c(c(c(pseudo_y1),c(pseudo_y22))%*%C_half)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
  }
  var_U_beta_theta_func<-function(
    UKBB_pop,
    beta,
    study_info,
    theta_UKBB_GPC=NULL){
    N_Pop<-nrow(UKBB_pop)
    xtbeta<-(UKBB_pop[,-1]%*%beta)
    var_11<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta))
    u2_theta_coef<-sapply(1:N_SNP, function(snp_id){
      if(is.null(var_GPC)){
        xtgamma_id<-c(UKBB_pop[,paste0("SNP",snp_id)]*study_info[[snp_id]]$Coeff)
      }else{
        xtgamma_id<-c((UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)))
      }
      xtbeta-xtgamma_id
    }) #col is SNP #row is sample 
    var_22<-crossprod(u2_theta_coef*UKBB_pop[,var_SNP],u2_theta_coef*UKBB_pop[,var_SNP])
    var_12<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),u2_theta_coef*UKBB_pop[,var_SNP])
    (1/N_Pop)*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
  }
  
  #### adaptive lasso with fast lasso computation 
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta0)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop*var(UKBB_pop[,1]) 
  ## Seems to be correlation matrix.
  ## Have taken the scaling 
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11)
  var_U2<-diag(sapply(1:N_SNP,function(i){
    (study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2/(N_Pop)
  }))
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  non_penalize_BMI<-c(rep(1,N_SNP),0)
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01,penalty.factor = non_penalize_BMI)
  
  gamma_adaptivelasso<-1/2
  w_adaptive<-c(1/(abs(coef(ridge_fit)[-1])[1:N_SNP])^gamma_adaptivelasso,0)
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
  beta_initial<-coef(ridge_fit)[-1]
  
  
  var_1st_U_beta_theta<-var_U_beta_theta_func(UKBB_pop = UKBB_pop,
    beta = beta_initial,
    study_info = study_info)
  diag_theta_sd<-sapply(1:N_SNP,function(i){
    sqrt(study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))
  })
  cov_theta<-t(diag_theta_sd*cor(UKBB_pop[,var_SNP]))*diag_theta_sd/N_Pop
  var_2nd_grad_times_theta_hat = adiag(matrix(0,nrow = len_U1,ncol = len_U1),cov_theta)
  inv_C = var_1st_U_beta_theta + var_2nd_grad_times_theta_hat
  C_all<-solve(inv_C)
  C_half<-expm::sqrtm(C_all)
  
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  lambda_list<-adaptivelasso_fit$lambda
  
  
  
  
  if(!cv_with_individual){
    beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
    index_nonzero<-c(which(beta[1:N_SNP]!=0),length(beta))
    
    if(length(index_nonzero) == 1){
      cv_with_individual = TRUE
    }else{
      xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
      xtx2<-t(UKBB_pop[,-1])%*%UKBB_pop[,var_SNP]/N_Pop
      Sigsum_half<-cbind(xtx,xtx2)%*%C_half
      Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
      Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
      inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
      
      index_nonzero<-index_nonzero[-length(index_nonzero)]
      final_v<-diag(inv_Sigsum_scaled_nonzero)
      final_v<-final_v[-length(final_v)]
      aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
      
      if(p_val_cut == "Method1"){
        candidate_pos<-index_nonzero[which(aa_final<0.05/N_SNP)]
      }else if(p_val_cut == "Method2"){
        candidate_pos<-index_nonzero[which(aa_final<0.05/length(index_nonzero))]
      }else{
        stop("wrong p_val_cut")
      }
      
      if(length(candidate_pos) == 0){
        cv_with_individual = TRUE
      }
    }
  }
  
  
  
  
  
  if(cv_with_individual){
    
    library(caret)
    index_fold<-createFolds(as.numeric(UKBB_pop[,1]>0),k = kfolds)
    a<-sapply(1:kfolds, function(cur_fold){
      index_test<-index_fold[[cur_fold]]
      UKBB_pop_train<-UKBB_pop[-index_test,]
      UKBB_pop_test<-UKBB_pop[index_test,]
      
      pseudo_Xy_list_train<-pseudo_Xy(C_half,UKBB_pop_train,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
      initial_sf_train<-nrow(UKBB_pop_train)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
      pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
      pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
      cv_beta<-lapply(lambda_list,function(cur_lam){
        cv_adaptivelasso_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive,lambda = cur_lam)
        coef(cv_adaptivelasso_fit)[-1]
      })
      cv_mse_sum<-sapply(1:length(cv_beta), function(kk){
        cur_beta<-cv_beta[[kk]]
        sum((UKBB_pop_test[,-1]%*%cur_beta - UKBB_pop_test[,1])^2)
      })
      cv_mse_sum
    })
    beta<-coef(adaptivelasso_fit,s=lambda_list[which.min(rowSums(a))])[-1]
    index_nonzero<-c(which(beta[1:N_SNP]!=0),length(beta))
    if(length(index_nonzero) == 1){
      beta<-coef(adaptivelasso_fit,s=lambda_list[which(adaptivelasso_fit$nzero>1)[which.min(rowSums(a)[which(adaptivelasso_fit$nzero>1)])]])[-1]
      index_nonzero<-c(which(beta[1:N_SNP]!=0),length(beta))
    }
    
    xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
    xtx2<-t(UKBB_pop[,-1])%*%UKBB_pop[,var_SNP]/N_Pop
    Sigsum_half<-cbind(xtx,xtx2)%*%C_half
    Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
    Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
    inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
    
    index_nonzero<-index_nonzero[-length(index_nonzero)]
    final_v<-diag(inv_Sigsum_scaled_nonzero)
    final_v<-final_v[-length(final_v)]
    aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
    
    
    if(p_val_cut == "Method1"){
      candidate_pos<-index_nonzero[which(aa_final<0.05/N_SNP)]
    }else if(p_val_cut == "Method2"){
      candidate_pos<-index_nonzero[which(aa_final<0.05/length(index_nonzero))]
    }else{
      stop("wrong p_val_cut")
    }
  }
  
  
  
  if(length(candidate_pos) > 1){
    candidate_cor<-cor(UKBB_pop[,c(candidate_pos+1)])
    diag(candidate_cor)<-0
    
    
    if(max(abs(candidate_cor))>0.1){
      
      UKBB_pop_orig<-UKBB_pop[,c(candidate_pos+1)]
      
      for(i in 1:ncol(UKBB_pop_orig)){
        cur_snp<-UKBB_pop_orig[,i]
        UKBB_pop_orig[which(cur_snp == max(cur_snp)),i] = 2
        UKBB_pop_orig[which(cur_snp == min(cur_snp)),i] = 0
        UKBB_pop_orig[which(cur_snp > min(cur_snp) & cur_snp < max(cur_snp)),i] = 1
      }
      candidate_EAF<-colsums(UKBB_pop_orig)/nrow(UKBB_pop_orig)/2
      
      weak_among_candidate<-which(candidate_EAF < 0.05)
      if(filter_by_EAF){
        if(length(weak_among_candidate) > 1){
          gamma_adaptivelasso<-1/2
          w_adaptive_candidate<-c(1/(abs(beta[candidate_pos[weak_among_candidate] ]))^gamma_adaptivelasso,0)
          
          tmp_count<-sapply(1:10, function(j){
            lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,c((candidate_pos[weak_among_candidate]+1),ncol(UKBB_pop)) ],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive_candidate)
            tmptmp<-abs(coef(lasso_fit_candidate,s='lambda.min')[-1])
            tmptmp<-tmptmp[-length(tmptmp)]
            tmptmp>1e-4
          })
          
          confident_pos_weak<-(candidate_pos[weak_among_candidate])[which(rowMeans(tmp_count) >= 0.8)]
          if(length(confident_pos_weak) == 0){
            z<-c()
            for(i in weak_among_candidate){
              z<-c(z,study_info[[candidate_pos[weak_among_candidate]]]$Coeff^2/study_info[[candidate_pos[weak_among_candidate]]]$Covariance)
            }
            confident_pos_weak<-weak_among_candidate[which.max(z)]
          }
        }else if(length(weak_among_candidate) == 1){
          confident_pos_weak<-candidate_pos[weak_among_candidate]
        }else{
          confident_pos_weak<-c()
        }
        
        true_weak<-c()
        if(length(weak_among_candidate) > 0){
          for( weakone in weak_among_candidate){
            if(max(abs(candidate_cor[,weakone]))<0.1){
              true_weak<-c(true_weak,weakone)      
            }
          }
        }
        if(length(true_weak)>0) true_weak<-candidate_pos[true_weak]
        
        confident_pos_weak <-union(confident_pos_weak,true_weak)
      }else{
        confident_pos_weak <-candidate_pos[weak_among_candidate]
      }
      
      gamma_adaptivelasso<-2
      w_adaptive_candidate<-c(1/(abs(beta[candidate_pos]))^gamma_adaptivelasso,0)
      w_adaptive_candidate[which(candidate_pos%in%confident_pos_weak)]<-0
      tmp_count<-sapply(1:10, function(j){
        lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,c( (candidate_pos+1),ncol(UKBB_pop)) ] ,y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 0.9,penalty.factor = w_adaptive_candidate)
        tmptmp<-abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
        tmptmp<-tmptmp[-length(tmptmp)]
        tmptmp>1e-4
        #coef(lasso_fit_candidate,s='lambda.min')[-1]!=0
      })
      confident_pos<-candidate_pos[which(rowMeans(tmp_count) >= 0.8)]
      
      confident_pos<-union(confident_pos,confident_pos_weak)
      
    }else{
      #print("NOT USED")
      confident_pos<-candidate_pos
    }
  }else{
    confident_pos<-candidate_pos
  }
  
  if(filter_index){
    pos2<-index_filter[confident_pos]
    pos<-index_filter[confident_pos]
  }else{
    pos<-confident_pos
    pos2<-confident_pos
  }
  #print(paste0("candidate",paste0(index_filter[candidate_pos],collapse = " ")))
  #print(paste0("confid",paste0(pos,collapse = " ")))
  #print(pos)
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos2"=pos2,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive,
    "index_filter"=index_filter)
}





adaptiveGMMlassoCBMI<-function(UKBB_pop,BMI,study_info,ld_cut = 0.9,cor_cut=0.9,filter_index=TRUE,filter_by_EAF=FALSE){
  UKBB_cor<-cor(UKBB_pop[,-1])
  diag(UKBB_cor) = 0
  if(filter_index){
    ### Trick 1 index filtering 
    index_filter<-c()
    pval_list<-c()
    for(i in 1:length(study_info)){
      #if(study_info[[i]]$P<p_val_cut){
      index_filter<-c(index_filter,i)
      z_here<-study_info[[i]]$Coeff^2/study_info[[i]]$Covariance
      pval_list<-c(pval_list,z_here)
      #}
    }
    index_filter0<-index_filter
    pval_list0<-pval_list
    index_filter2<-c()
    while (length(index_filter) > 1){
      UKBB_cor_i<-UKBB_cor[index_filter,index_filter]
      cur_i<-which.max(pval_list)
      index_filter2<-c(index_filter2,index_filter[cur_i])
      #print(index_filter[cur_i])
      rm_index<-c(cur_i,which(abs(UKBB_cor_i[cur_i,])>ld_cut))
      #print(index_filter[rm_index])
      index_filter<-index_filter[-rm_index]
      pval_list<-pval_list[-rm_index]
    }
    
    #Nonnull_index2<-which(index_filter2%in%Nonnull_index)
    UKBB_pop<-UKBB_pop[,c(1,(index_filter2+1))]
    study_info<-study_info[index_filter2]
    index_filter<-index_filter2
  }
  
  N_SNP<-ncol(UKBB_pop)-1
  colnames(UKBB_pop)[-1]<-paste0("SNP",1:(N_SNP))
  colnames(UKBB_pop)[1]<-"Y"
  var_SNP<-paste0("SNP",1:(N_SNP))
  len_SNP<-length(var_SNP)
  var_GPC<-NULL
  len_GPC<-length(var_GPC)
  var_nonSNP<-c(var_GPC,"BMI")
  len_nonSNP<-length(var_nonSNP)
  var_names_full_fit <- c(var_SNP,var_nonSNP)
  colname_UKBB<-colnames(UKBB_pop)
  len_U1 = len_SNP+len_nonSNP
  len_U2 = len_SNP
  len_beta = len_U1
  len_theta = len_SNP+len_GPC
  len_U = len_U1 + len_U2
  N_Pop<-nrow(UKBB_pop)
  UKBB_pop<-cbind(UKBB_pop,BMI)
  colnames(UKBB_pop)[ncol(UKBB_pop)]<-"BMI"
  ############## Beta
  pseudo_Xy<-function(
    C_half,UKBB_pop,
    var_SNP,var_GPC,
    N_SNP,theta_UKBB_GPC,study_info){
    N_Pop<-nrow(UKBB_pop)
    pseudo_X<-C_half%*%rbind(t(UKBB_pop[,-1]),t(UKBB_pop[,var_SNP]))%*%UKBB_pop[,-1]
    pseudo_y1<-crossprod(UKBB_pop[,1],UKBB_pop[,-1])
    pseudo_y22<-sapply(1:N_SNP, function(snp_id){
      u2_id<-UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
      c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
    })
    
    pseudo_y<-c(c(c(pseudo_y1),c(pseudo_y22))%*%C_half)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
  }
  var_U_beta_theta_func<-function(
    UKBB_pop,
    beta,
    study_info,
    theta_UKBB_GPC=NULL){
    N_Pop<-nrow(UKBB_pop)
    xtbeta<-(UKBB_pop[,-1]%*%beta)
    var_11<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta))
    u2_theta_coef<-sapply(1:N_SNP, function(snp_id){
      if(is.null(var_GPC)){
        xtgamma_id<-c(UKBB_pop[,paste0("SNP",snp_id)]*study_info[[snp_id]]$Coeff)
      }else{
        xtgamma_id<-c((UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)))
      }
      xtbeta-xtgamma_id
    }) #col is SNP #row is sample 
    var_22<-crossprod(u2_theta_coef*UKBB_pop[,var_SNP],u2_theta_coef*UKBB_pop[,var_SNP])
    var_12<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),u2_theta_coef*UKBB_pop[,var_SNP])
    (1/N_Pop)*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
  }
  
  #### adaptive lasso with fast lasso computation 
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta0)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop*var(UKBB_pop[,1]) 
  ## Seems to be correlation matrix.
  ## Have taken the scaling 
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11)
  var_U2<-diag(sapply(1:N_SNP,function(i){
    (study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2/(N_Pop)
  }))
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  non_penalize_BMI<-c(rep(1,N_SNP),0)
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01,penalty.factor = non_penalize_BMI)
  
  gamma_adaptivelasso<-1/2
  w_adaptive<-c(1/(abs(coef(ridge_fit)[-1])[1:N_SNP])^gamma_adaptivelasso,0)
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
  beta_initial<-coef(ridge_fit)[-1]
  
  
  var_1st_U_beta_theta<-var_U_beta_theta_func(UKBB_pop = UKBB_pop,
    beta = beta_initial,
    study_info = study_info)
  #inv_C_22<-diag(sapply(1:N_SNP,function(i){
  #  (study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2/(N_Pop)
  #}))
  
  diag_theta_sd<-sapply(1:N_SNP,function(i){
    sqrt(study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))
  })
  cov_theta<-t(diag_theta_sd*cor(UKBB_pop[,var_SNP]))*diag_theta_sd/N_Pop
  
  var_2nd_grad_times_theta_hat = adiag(matrix(0,nrow = len_U1,ncol = len_U1),cov_theta)
  
  inv_C = var_1st_U_beta_theta + var_2nd_grad_times_theta_hat
  C_all<-solve(inv_C)
  C_half<-expm::sqrtm(C_all)
  
  
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  lambda_list<-adaptivelasso_fit$lambda
  UKBB_cor<-cor(UKBB_pop[,var_SNP])
  beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
  
  
  while(1){
    index_nonzero_i<-which(beta[1:N_SNP]!=0)
    nonzero_cor_i<-UKBB_cor[index_nonzero_i,index_nonzero_i]
    nonzero_cor_i[nonzero_cor_i == 1] = 0
    tmp<-which(abs(nonzero_cor_i)>cor_cut,arr.ind = T)
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
  
  index_nonzero<-c(which(beta[1:N_SNP]!=0),length(beta))
  xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
  xtx2<-t(UKBB_pop[,-1])%*%UKBB_pop[,var_SNP]/N_Pop
  Sigsum_half<-cbind(xtx,xtx2)%*%C_half
  Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
  Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
  inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
  
  #W1 = xtx/var(UKBB_pop[,1])
  #W2 = xtx%*%C_22%*%xtx
  #W = W1+W2
  #W_nonzero = W[index_nonzero,index_nonzero]
  #print(beta[index_nonzero])
  #final_v<-diag(inv_Sigsum_scaled_nonzero%*%W_nonzero%*%inv_Sigsum_scaled_nonzero)
  index_nonzero<-index_nonzero[-length(index_nonzero)]
  final_v<-diag(inv_Sigsum_scaled_nonzero)
  final_v<-final_v[-length(final_v)]
  aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
  
  
  #candidate_pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  candidate_pos<-index_nonzero[which(aa_final<0.05/length(index_nonzero))]
  if(length(candidate_pos) > 1){
    candidate_cor<-cor(UKBB_pop[,c(candidate_pos+1)])
    diag(candidate_cor)<-0
    
    
    if(max(abs(candidate_cor))>0.1){
      
      UKBB_pop_orig<-UKBB_pop[,c(candidate_pos+1)]
      
      for(i in 1:ncol(UKBB_pop_orig)){
        cur_snp<-UKBB_pop_orig[,i]
        UKBB_pop_orig[which(cur_snp == max(cur_snp)),i] = 2
        UKBB_pop_orig[which(cur_snp == min(cur_snp)),i] = 0
        UKBB_pop_orig[which(cur_snp > min(cur_snp) & cur_snp < max(cur_snp)),i] = 1
      }
      candidate_EAF<-colsums(UKBB_pop_orig)/nrow(UKBB_pop_orig)/2
      
      weak_among_candidate<-which(candidate_EAF < 0.05)
      if(filter_by_EAF){
        if(length(weak_among_candidate) > 1){
          gamma_adaptivelasso<-1/2
          w_adaptive_candidate<-c(1/(abs(beta[candidate_pos[weak_among_candidate] ]))^gamma_adaptivelasso,0)
          
          tmp_count<-sapply(1:10, function(j){
            lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,c((candidate_pos[weak_among_candidate]+1),ncol(UKBB_pop)) ],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive_candidate)
            tmptmp<-abs(coef(lasso_fit_candidate,s='lambda.min')[-1])
            tmptmp<-tmptmp[-length(tmptmp)]
            tmptmp>1e-4
          })
          
          confident_pos_weak<-(candidate_pos[weak_among_candidate])[which(rowMeans(tmp_count) >= 0.8)]
          if(length(confident_pos_weak) == 0){
            z<-c()
            for(i in weak_among_candidate){
              z<-c(z,study_info[[candidate_pos[weak_among_candidate]]]$Coeff^2/study_info[[candidate_pos[weak_among_candidate]]]$Covariance)
            }
            confident_pos_weak<-weak_among_candidate[which.max(z)]
          }
        }else if(length(weak_among_candidate) == 1){
          confident_pos_weak<-candidate_pos[weak_among_candidate]
        }else{
          confident_pos_weak<-c()
        }
        
        true_weak<-c()
        if(length(weak_among_candidate) > 0){
          for( weakone in weak_among_candidate){
            if(max(abs(candidate_cor[,weakone]))<0.1){
              true_weak<-c(true_weak,weakone)      
            }
          }
        }
        if(length(true_weak)>0) true_weak<-candidate_pos[true_weak]
        
        confident_pos_weak <-union(confident_pos_weak,true_weak)
      }else{
        confident_pos_weak <-candidate_pos[weak_among_candidate]
      }
      
      gamma_adaptivelasso<-2
      w_adaptive_candidate<-c(1/(abs(beta[candidate_pos]))^gamma_adaptivelasso,0)
      w_adaptive_candidate[which(candidate_pos%in%confident_pos_weak)]<-0
      tmp_count<-sapply(1:10, function(j){
        lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,c( (candidate_pos+1),ncol(UKBB_pop)) ] ,y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 0.9,penalty.factor = w_adaptive_candidate)
        tmptmp<-abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
        tmptmp<-tmptmp[-length(tmptmp)]
        tmptmp>1e-4
        #coef(lasso_fit_candidate,s='lambda.min')[-1]!=0
      })
      confident_pos<-candidate_pos[which(rowMeans(tmp_count) >= 0.8)]
      
      confident_pos<-union(confident_pos,confident_pos_weak)
      
    }else{
      print("NOT USED")
      confident_pos<-candidate_pos
    }
  }else{
    confident_pos<-candidate_pos
  }
  
  if(filter_index){
    pos2<-index_filter[confident_pos]
    pos<-index_filter[confident_pos]
  }else{
    pos<-confident_pos
    pos2<-confident_pos
  }
  print(paste0("candidate",paste0(index_filter[candidate_pos],collapse = " ")))
  print(paste0("confid",paste0(pos,collapse = " ")))
  #print(pos)
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos2"=pos2,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive)
}



adaptiveGMMlassoC<-function(UKBB_pop,study_info,
  ld_cut = 0.9,
  cor_cut=0.9,
  filter_index=TRUE,
  filter_by_EAF=FALSE,
  cv_with_individual=FALSE,
  kfolds=10){
  UKBB_cor<-cor(UKBB_pop[,-1])
  diag(UKBB_cor) = 0
  if(filter_index){
    ### Trick 1 index filtering 
    index_filter<-c()
    pval_list<-c()
    for(i in 1:length(study_info)){
      #if(study_info[[i]]$P<p_val_cut){
      index_filter<-c(index_filter,i)
      z_here<-study_info[[i]]$Coeff^2/study_info[[i]]$Covariance
      pval_list<-c(pval_list,z_here)
      #}
    }
    index_filter0<-index_filter
    pval_list0<-pval_list
    index_filter2<-c()
    while (length(index_filter) > 1){
      UKBB_cor_i<-UKBB_cor[index_filter,index_filter]
      cur_i<-which.max(pval_list)
      index_filter2<-c(index_filter2,index_filter[cur_i])
      #print(index_filter[cur_i])
      rm_index<-c(cur_i,which(abs(UKBB_cor_i[cur_i,])>ld_cut))
      #print(index_filter[rm_index])
      index_filter<-index_filter[-rm_index]
      pval_list<-pval_list[-rm_index]
    }
    
    #Nonnull_index2<-which(index_filter2%in%Nonnull_index)
    UKBB_pop<-UKBB_pop[,c(1,(index_filter2+1))]
    study_info<-study_info[index_filter2]
    index_filter<-index_filter2
  }
  
  N_SNP<-ncol(UKBB_pop)-1
  colnames(UKBB_pop)[-1]<-paste0("SNP",1:(N_SNP))
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
  
  ############## Beta
  pseudo_Xy<-function(
    C_half,UKBB_pop,
    var_SNP,var_GPC,
    N_SNP,theta_UKBB_GPC,study_info){
    N_Pop<-nrow(UKBB_pop)
    pseudo_X<-C_half%*%rbind(t(UKBB_pop[,-1]),t(UKBB_pop[,var_SNP]))%*%UKBB_pop[,-1]
    pseudo_y1<-crossprod(UKBB_pop[,1],UKBB_pop[,-1])
    pseudo_y22<-sapply(1:N_SNP, function(snp_id){
      u2_id<-UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
      c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
    })
    
    pseudo_y<-c(c(c(pseudo_y1),c(pseudo_y22))%*%C_half)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
  }
  
  
  var_U_beta_theta_func<-function(
    UKBB_pop,
    beta,
    study_info,
    theta_UKBB_GPC=NULL){
    N_Pop<-nrow(UKBB_pop)
    xtbeta<-(UKBB_pop[,-1]%*%beta)
    var_11<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta))
    u2_theta_coef<-sapply(1:N_SNP, function(snp_id){
      if(is.null(var_GPC)){
        xtgamma_id<-c(UKBB_pop[,paste0("SNP",snp_id)]*study_info[[snp_id]]$Coeff)
      }else{
        xtgamma_id<-c((UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)))
      }
      xtbeta-xtgamma_id
    }) #col is SNP #row is sample 
    var_22<-crossprod(u2_theta_coef*UKBB_pop[,var_SNP],u2_theta_coef*UKBB_pop[,var_SNP])
    var_12<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),u2_theta_coef*UKBB_pop[,var_SNP])
    (1/N_Pop)*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
  }
  
  #### adaptive lasso with fast lasso computation 
  #var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  #var_U1<-crossprod(var_11_half,var_11_half)/N_Pop
  #C_11<-solve(var_U1)
  #C_11_half<-expm::sqrtm(C_11)
  
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop*var(UKBB_pop[,1])
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11) 
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01)#,penalty.factor = w_adaptive0)
  
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
  beta_initial<-coef(ridge_fit)[-1]
  
  var_1st_U_beta_theta<-var_U_beta_theta_func(UKBB_pop =UKBB_pop,
    beta = beta_initial,
    study_info = study_info)
  
  diag_theta_sd<-sapply(1:N_SNP,function(i){
    sqrt(study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))
  })
  cov_theta<-t(diag_theta_sd*cor(UKBB_pop[,var_SNP]))*diag_theta_sd/N_Pop
  
  var_2nd_grad_times_theta_hat = adiag(matrix(0,nrow = len_U1,ncol = len_U1),cov_theta)
  
  inv_C = var_1st_U_beta_theta + var_2nd_grad_times_theta_hat
  C_all<-solve(inv_C)
  C_half<-expm::sqrtm(C_all)
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  
  lambda_list<-adaptivelasso_fit$lambda
  
  if(!cv_with_individual){
    beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
    index_nonzero<-which(beta!=0)
    if(length(index_nonzero) == 0){
      cv_with_individual = TRUE
    }else{
      xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
      Sigsum_half<-cbind(xtx,xtx)%*%C_half
      Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
      Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
      inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
      final_v<-diag(inv_Sigsum_scaled_nonzero)
      aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
      candidate_pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
      if(length(candidate_pos) == 0){
        cv_with_individual = TRUE
      }
    }
  }
  
  
  if(cv_with_individual){
    #index11<-which(adaptivelasso_fit$nzero>1)
    #index11<-unique(c(max(min(index11)-1,1),index11))
    #lambda_list<-lambda_list[index11]
    
    library(caret)
    index_fold<-createFolds(as.numeric(UKBB_pop[,1]>0),k = kfolds)
    a<-sapply(1:kfolds, function(cur_fold){
      index_test<-index_fold[[cur_fold]]
      UKBB_pop_train<-UKBB_pop[-index_test,]
      UKBB_pop_test<-UKBB_pop[index_test,]
      pseudo_Xy_list_train<-pseudo_Xy(C_half,UKBB_pop_train,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
      initial_sf_train<-nrow(UKBB_pop_train)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
      pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
      pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
      cv_beta<-lapply(lambda_list,function(cur_lam){
        cv_adaptivelasso_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive,lambda = cur_lam)
        coef(cv_adaptivelasso_fit)[-1]
      })
      cv_mse_sum<-sapply(1:length(cv_beta), function(kk){
        cur_beta<-cv_beta[[kk]]
        mean((UKBB_pop_test[,-1]%*%cur_beta - UKBB_pop_test[,1])^2)
      })
      cv_mse_sum
    })

    beta<-coef(adaptivelasso_fit,s=lambda_list[which.min(rowSums(a))])[-1]
    index_nonzero<-which(beta!=0)
    ########??????????????????
    if(length(index_nonzero) == 0){
      beta<-coef(adaptivelasso_fit,s=lambda_list[which(adaptivelasso_fit$nzero>0)[which.min(rowSums(a)[which(adaptivelasso_fit$nzero>0)])]])[-1]
      index_nonzero<-which(beta!=0)
    }
    xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
    Sigsum_half<-cbind(xtx,xtx)%*%C_half
    Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
    Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
    inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
    final_v<-diag(inv_Sigsum_scaled_nonzero)
    aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
    candidate_pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  }
  
  if(length(candidate_pos) > 1){
    candidate_cor<-cor(UKBB_pop[,c(candidate_pos+1)])
    diag(candidate_cor)<-0
    
    
    if(max(abs(candidate_cor))>0.1){
      
      
      
      UKBB_pop_orig<-UKBB_pop[,c(candidate_pos+1)]
      
      for(i in 1:ncol(UKBB_pop_orig)){
        cur_snp<-UKBB_pop_orig[,i]
        UKBB_pop_orig[which(cur_snp == max(cur_snp)),i] = 2
        UKBB_pop_orig[which(cur_snp == min(cur_snp)),i] = 0
        UKBB_pop_orig[which(cur_snp > min(cur_snp) & cur_snp < max(cur_snp)),i] = 1
      }
      candidate_EAF<-colsums(UKBB_pop_orig)/nrow(UKBB_pop_orig)/2
      
      weak_among_candidate<-which(candidate_EAF < 0.05)
      if(filter_by_EAF){
        if(length(weak_among_candidate) > 1){
          gamma_adaptivelasso<-1/2
          w_adaptive_candidate<-1/(abs(beta[candidate_pos[weak_among_candidate] ]))^gamma_adaptivelasso
          
          
          tmp_count<-sapply(1:10, function(j){
            lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos[weak_among_candidate]+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive_candidate)
            abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
          })
          
          confident_pos_weak<-(candidate_pos[weak_among_candidate])[which(rowMeans(tmp_count) >= 0.8)]
          if(length(confident_pos_weak) == 0){
            z<-c()
            for(i in weak_among_candidate){
              z<-c(z,study_info[[candidate_pos[weak_among_candidate]]]$Coeff^2/study_info[[candidate_pos[weak_among_candidate]]]$Covariance)
            }
            confident_pos_weak<-weak_among_candidate[which.max(z)]
          }
        }else if(length(weak_among_candidate) == 1){
          confident_pos_weak<-candidate_pos[weak_among_candidate]
        }else{
          confident_pos_weak<-c()
        }
        
        true_weak<-c()
        if(length(weak_among_candidate) > 0){
          for( weakone in weak_among_candidate){
            if(max(abs(candidate_cor[,weakone]))<0.1){
              true_weak<-c(true_weak,weakone)      
            }
          }
        }
        if(length(true_weak)>0) true_weak<-candidate_pos[true_weak]
        
        confident_pos_weak <-union(confident_pos_weak,true_weak)
      }else{
        confident_pos_weak <-candidate_pos[weak_among_candidate]
      }
      
      gamma_adaptivelasso<-2
      w_adaptive_candidate<-1/(abs(beta[candidate_pos]))^gamma_adaptivelasso
      w_adaptive_candidate[which(candidate_pos%in%confident_pos_weak)]<-0
      tmp_count<-sapply(1:10, function(j){
        lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 0.9,penalty.factor = w_adaptive_candidate)
        abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
        #coef(lasso_fit_candidate,s='lambda.min')[-1]!=0
      })
      confident_pos<-candidate_pos[which(rowMeans(tmp_count) >= 0.8)]
      
      confident_pos<-union(confident_pos,confident_pos_weak)
      
    }else{
      print("NOT USED")
      confident_pos<-candidate_pos
    }
  }else{
    confident_pos<-candidate_pos
  }
  
  if(filter_index){
    pos<-index_filter[confident_pos]
    pos_bf<-index_filter[candidate_pos]
  }else{
    pos<-confident_pos
    pos_bf<-candidate_pos
  }
  #print(paste0("candidate",paste0(index_filter[candidate_pos],collapse = " ")))
  #print(paste0("confid",paste0(pos,collapse = " ")))
  #print(pos)
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos_bf"=pos_bf,
    "aa_final"=aa_final,
    "candidate"=index_filter[candidate_pos],
    "final_v"=final_v,
    "w_adaptive"=w_adaptive)
}


adaptiveGMMlasso_v2_ind<-function(UKBB_pop,study_info,ld_cut = 0.8,cor_cut=0.8,filter_index=TRUE,filter_by_EAF=FALSE,cv_with_individual=TRUE){
  UKBB_cor<-cor(UKBB_pop[,-1])
  diag(UKBB_cor) = 0
  if(filter_index){
    ### Trick 1 index filtering 
    index_filter<-c()
    pval_list<-c()
    for(i in 1:length(study_info)){
      #if(study_info[[i]]$P<p_val_cut){
      index_filter<-c(index_filter,i)
      z_here<-study_info[[i]]$Coeff^2/study_info[[i]]$Covariance
      pval_list<-c(pval_list,z_here)
      #}
    }
    index_filter0<-index_filter
    pval_list0<-pval_list
    index_filter2<-c()
    while (length(index_filter) > 1){
      UKBB_cor_i<-UKBB_cor[index_filter,index_filter]
      cur_i<-which.max(pval_list)
      index_filter2<-c(index_filter2,index_filter[cur_i])
      #print(index_filter[cur_i])
      rm_index<-c(cur_i,which(abs(UKBB_cor_i[cur_i,])>ld_cut))
      #print(index_filter[rm_index])
      index_filter<-index_filter[-rm_index]
      pval_list<-pval_list[-rm_index]
    }
    
    #Nonnull_index2<-which(index_filter2%in%Nonnull_index)
    UKBB_pop<-UKBB_pop[,c(1,(index_filter2+1))]
    study_info<-study_info[index_filter2]
    index_filter<-index_filter2
  }
  
  N_SNP<-ncol(UKBB_pop)-1
  colnames(UKBB_pop)[-1]<-paste0("SNP",1:(N_SNP))
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
  
  #### adaptive lasso with fast lasso computation 
  ## The initial is important
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop*var(UKBB_pop[,1])
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11) 
  
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01)#,penalty.factor = w_adaptive0)
  
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
  diag_theta_sd<-sapply(1:N_SNP,function(i){
    sqrt(study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))
  })
  cov_theta<-t(diag_theta_sd*cor(UKBB_pop[,var_SNP]))*diag_theta_sd/N_Pop
  C_22<-solve(cov_theta)
  C_22_half<-expm::sqrtm(C_22)
  C_half<-adiag(C_11_half,C_22_half)
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  lambda_list<-adaptivelasso_fit$lambda
  #UKBB_cor<-cor(UKBB_pop[,-1])
  
  if(cv_with_individual){
    index11<-which(adaptivelasso_fit$nzero>1)
    index11<-unique(c(max(min(index11)-1,1),index11))
    lambda_list<-lambda_list[index11]
    
    kfolds<-10
    library(caret)
    index_fold<-createFolds(as.numeric(UKBB_pop[,1]>0),k = kfolds)
    
    a<-sapply(1:kfolds, function(cur_fold){
      index_test<-index_fold[[cur_fold]]
      UKBB_pop_train<-UKBB_pop[-index_test,]
      UKBB_pop_test<-UKBB_pop[index_test,]
      pseudo_Xy_list_train<-pseudo_Xy(C_half,UKBB_pop_train,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
      initial_sf_train<-nrow(UKBB_pop_train)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
      pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
      pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
      cv_beta<-lapply(lambda_list,function(cur_lam){
        cv_adaptivelasso_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive,lambda = cur_lam)
        coef(cv_adaptivelasso_fit)[-1]
      })
      cv_mse_sum<-sapply(1:length(cv_beta), function(kk){
        cur_beta<-cv_beta[[kk]]
        sum((UKBB_pop_test[,-1]%*%cur_beta - UKBB_pop_test[,1])^2)
      })
      cv_mse_sum
    })
    
    beta<-coef(adaptivelasso_fit,s=lambda_list[which.min(rowSums(a))])[-1]
  }else{
    beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
  }
  
  
  index_nonzero<-which(beta!=0)
  xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
  Sigsum_half<-cbind(xtx,xtx)%*%C_half
  Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
  Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
  inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
  
  #W1 = xtx/var(UKBB_pop[,1])
  #W2 = xtx%*%C_22%*%xtx
  #W = W1+W2
  #W_nonzero = W[index_nonzero,index_nonzero]
  #print(beta[index_nonzero])
  #final_v<-diag(inv_Sigsum_scaled_nonzero%*%W_nonzero%*%inv_Sigsum_scaled_nonzero)
  final_v<-diag(inv_Sigsum_scaled_nonzero)
  aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
  
  
  candidate_pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  #candidate_pos<-index_nonzero[which(aa_final<0.05/length(index_nonzero))]
  if(length(candidate_pos) > 1){
    candidate_cor<-cor(UKBB_pop[,c(candidate_pos+1)])
    diag(candidate_cor)<-0
    
    
    if(max(abs(candidate_cor))>0.1){
      
      
      
      UKBB_pop_orig<-UKBB_pop[,c(candidate_pos+1)]
      
      for(i in 1:ncol(UKBB_pop_orig)){
        cur_snp<-UKBB_pop_orig[,i]
        UKBB_pop_orig[which(cur_snp == max(cur_snp)),i] = 2
        UKBB_pop_orig[which(cur_snp == min(cur_snp)),i] = 0
        UKBB_pop_orig[which(cur_snp > min(cur_snp) & cur_snp < max(cur_snp)),i] = 1
      }
      candidate_EAF<-colsums(UKBB_pop_orig)/nrow(UKBB_pop_orig)/2
      
      weak_among_candidate<-which(candidate_EAF < 0.05)
      if(filter_by_EAF){
        if(length(weak_among_candidate) > 1){
          gamma_adaptivelasso<-1/2
          w_adaptive_candidate<-1/(abs(beta[candidate_pos[weak_among_candidate] ]))^gamma_adaptivelasso
          
          
          tmp_count<-sapply(1:10, function(j){
            lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos[weak_among_candidate]+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive_candidate)
            abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
          })
          
          confident_pos_weak<-(candidate_pos[weak_among_candidate])[which(rowMeans(tmp_count) >= 0.8)]
          if(length(confident_pos_weak) == 0){
            z<-c()
            for(i in weak_among_candidate){
              z<-c(z,study_info[[candidate_pos[weak_among_candidate]]]$Coeff^2/study_info[[candidate_pos[weak_among_candidate]]]$Covariance)
            }
            confident_pos_weak<-weak_among_candidate[which.max(z)]
          }
        }else if(length(weak_among_candidate) == 1){
          confident_pos_weak<-candidate_pos[weak_among_candidate]
        }else{
          confident_pos_weak<-c()
        }
        
        true_weak<-c()
        if(length(weak_among_candidate) > 0){
          for( weakone in weak_among_candidate){
            if(max(abs(candidate_cor[,weakone]))<0.1){
              true_weak<-c(true_weak,weakone)      
            }
          }
        }
        if(length(true_weak)>0) true_weak<-candidate_pos[true_weak]
        
        confident_pos_weak <-union(confident_pos_weak,true_weak)
      }else{
        confident_pos_weak <-candidate_pos[weak_among_candidate]
      }
      
      gamma_adaptivelasso<-2
      w_adaptive_candidate<-1/(abs(beta[candidate_pos]))^gamma_adaptivelasso
      w_adaptive_candidate[which(candidate_pos%in%confident_pos_weak)]<-0
      tmp_count<-sapply(1:10, function(j){
        lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 0.9,penalty.factor = w_adaptive_candidate)
        abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
        #coef(lasso_fit_candidate,s='lambda.min')[-1]!=0
      })
      confident_pos<-candidate_pos[which(rowMeans(tmp_count) >= 0.8)]
      
      confident_pos<-union(confident_pos,confident_pos_weak)
      
    }else{
      print("NOT USED")
      confident_pos<-candidate_pos
    }
  }else{
    confident_pos<-candidate_pos
  }
  
  if(filter_index){
    pos<-index_filter[confident_pos]
    pos_bf<-index_filter[candidate_pos]
  }else{
    pos<-confident_pos
    pos_bf<-candidate_pos
  }
  #print(paste0("candidate",paste0(index_filter[candidate_pos],collapse = " ")))
  #print(paste0("confid",paste0(pos,collapse = " ")))
  #print(pos)
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos_bf"=pos_bf,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive)
}

adaptiveGMMlassoC_v2_ind<-function(UKBB_pop,study_info,ld_cut = 0.8,cor_cut=0.8,filter_index=TRUE,filter_by_EAF=FALSE,cv_with_individual=TRUE){
  UKBB_cor<-cor(UKBB_pop[,-1])
  diag(UKBB_cor) = 0
  if(filter_index){
    ### Trick 1 index filtering 
    index_filter<-c()
    pval_list<-c()
    for(i in 1:length(study_info)){
      #if(study_info[[i]]$P<p_val_cut){
      index_filter<-c(index_filter,i)
      z_here<-study_info[[i]]$Coeff^2/study_info[[i]]$Covariance
      pval_list<-c(pval_list,z_here)
      #}
    }
    index_filter0<-index_filter
    pval_list0<-pval_list
    index_filter2<-c()
    while (length(index_filter) > 1){
      UKBB_cor_i<-UKBB_cor[index_filter,index_filter]
      cur_i<-which.max(pval_list)
      index_filter2<-c(index_filter2,index_filter[cur_i])
      #print(index_filter[cur_i])
      rm_index<-c(cur_i,which(abs(UKBB_cor_i[cur_i,])>ld_cut))
      #print(index_filter[rm_index])
      index_filter<-index_filter[-rm_index]
      pval_list<-pval_list[-rm_index]
    }
    
    #Nonnull_index2<-which(index_filter2%in%Nonnull_index)
    UKBB_pop<-UKBB_pop[,c(1,(index_filter2+1))]
    study_info<-study_info[index_filter2]
    index_filter<-index_filter2
  }
  
  N_SNP<-ncol(UKBB_pop)-1
  colnames(UKBB_pop)[-1]<-paste0("SNP",1:(N_SNP))
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
  
  
  var_U_beta_theta_func<-function(
    UKBB_pop,
    beta,
    study_info,
    theta_UKBB_GPC=NULL){
    N_Pop<-nrow(UKBB_pop)
    xtbeta<-(UKBB_pop[,-1]%*%beta)
    var_11<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta))
    u2_theta_coef<-sapply(1:N_SNP, function(snp_id){
      if(is.null(var_GPC)){
        xtgamma_id<-c(UKBB_pop[,paste0("SNP",snp_id)]*study_info[[snp_id]]$Coeff)
      }else{
        xtgamma_id<-c((UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)))
      }
      xtbeta-xtgamma_id
    }) #col is SNP #row is sample 
    var_22<-crossprod(u2_theta_coef*UKBB_pop[,var_SNP],u2_theta_coef*UKBB_pop[,var_SNP])
    var_12<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),u2_theta_coef*UKBB_pop[,var_SNP])
    (1/N_Pop)*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
  }
  
  #### adaptive lasso with fast lasso computation 
  #var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  #var_U1<-crossprod(var_11_half,var_11_half)/N_Pop
  #C_11<-solve(var_U1)
  #C_11_half<-expm::sqrtm(C_11)
  
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop*var(UKBB_pop[,1])
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11) 
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01)#,penalty.factor = w_adaptive0)
  
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
  beta_initial<-coef(ridge_fit)[-1]
  
  var_1st_U_beta_theta<-var_U_beta_theta_func(UKBB_pop =UKBB_pop,
    beta = beta_initial,
    study_info = study_info)
  #inv_C_22<-diag(sapply(1:N_SNP,function(i){
  #  (study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2/(N_Pop)
  #}))
  
  diag_theta_sd<-sapply(1:N_SNP,function(i){
    sqrt(study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))
  })
  cov_theta<-t(diag_theta_sd*cor(UKBB_pop[,var_SNP]))*diag_theta_sd/N_Pop
  
  var_2nd_grad_times_theta_hat = adiag(matrix(0,nrow = len_U1,ncol = len_U1),cov_theta)
  
  inv_C = var_1st_U_beta_theta + var_2nd_grad_times_theta_hat
  C_all<-solve(inv_C)
  C_half<-expm::sqrtm(C_all)
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  
  lambda_list<-adaptivelasso_fit$lambda
  
  if(cv_with_individual){
    index11<-which(adaptivelasso_fit$nzero>1)
    index11<-unique(c(max(min(index11)-1,1),index11))
    lambda_list<-lambda_list[index11]
    
    kfolds<-10
    library(caret)
    index_fold<-createFolds(as.numeric(UKBB_pop[,1]>0),k = kfolds)
    
    a<-sapply(1:kfolds, function(cur_fold){
      index_test<-index_fold[[cur_fold]]
      UKBB_pop_train<-UKBB_pop[-index_test,]
      UKBB_pop_test<-UKBB_pop[index_test,]
      pseudo_Xy_list_train<-pseudo_Xy(C_half,UKBB_pop_train,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
      initial_sf_train<-nrow(UKBB_pop_train)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
      pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
      pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
      cv_beta<-lapply(lambda_list,function(cur_lam){
        cv_adaptivelasso_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive,lambda = cur_lam)
        coef(cv_adaptivelasso_fit)[-1]
      })
      cv_mse_sum<-sapply(1:length(cv_beta), function(kk){
        cur_beta<-cv_beta[[kk]]
        sum((UKBB_pop_test[,-1]%*%cur_beta - UKBB_pop_test[,1])^2)
      })
      cv_mse_sum
    })
    
    beta<-coef(adaptivelasso_fit,s=lambda_list[which.min(rowSums(a))])[-1]
  }else{
    beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
  }
  
  #UKBB_cor<-cor(UKBB_pop[,-1])
  #beta<-coef(adaptivelasso_fit,s =lambda_list[which(adaptivelasso_fit$nzero>3)[which.min(adaptivelasso_fit$cvm[which(adaptivelasso_fit$nzero>3)])]])[-1]
  
  index_nonzero<-which(beta!=0)
  xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
  Sigsum_half<-cbind(xtx,xtx)%*%C_half
  Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
  Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
  inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
  
  #W1 = xtx/var(UKBB_pop[,1])
  #W2 = xtx%*%C_22%*%xtx
  #W = W1+W2
  #W_nonzero = W[index_nonzero,index_nonzero]
  #print(beta[index_nonzero])
  #final_v<-diag(inv_Sigsum_scaled_nonzero%*%W_nonzero%*%inv_Sigsum_scaled_nonzero)
  final_v<-diag(inv_Sigsum_scaled_nonzero)
  aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
  
  
  candidate_pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  #candidate_pos<-index_nonzero[which(aa_final<0.05/length(index_nonzero))]
  if(length(candidate_pos) > 1){
    candidate_cor<-cor(UKBB_pop[,c(candidate_pos+1)])
    diag(candidate_cor)<-0
    
    
    if(max(abs(candidate_cor))>0.1){
      
      
      
      UKBB_pop_orig<-UKBB_pop[,c(candidate_pos+1)]
      
      for(i in 1:ncol(UKBB_pop_orig)){
        cur_snp<-UKBB_pop_orig[,i]
        UKBB_pop_orig[which(cur_snp == max(cur_snp)),i] = 2
        UKBB_pop_orig[which(cur_snp == min(cur_snp)),i] = 0
        UKBB_pop_orig[which(cur_snp > min(cur_snp) & cur_snp < max(cur_snp)),i] = 1
      }
      candidate_EAF<-colsums(UKBB_pop_orig)/nrow(UKBB_pop_orig)/2
      
      weak_among_candidate<-which(candidate_EAF < 0.05)
      if(filter_by_EAF){
        if(length(weak_among_candidate) > 1){
          gamma_adaptivelasso<-1/2
          w_adaptive_candidate<-1/(abs(beta[candidate_pos[weak_among_candidate] ]))^gamma_adaptivelasso
          
          
          tmp_count<-sapply(1:10, function(j){
            lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos[weak_among_candidate]+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive_candidate)
            abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
          })
          
          confident_pos_weak<-(candidate_pos[weak_among_candidate])[which(rowMeans(tmp_count) >= 0.8)]
          if(length(confident_pos_weak) == 0){
            z<-c()
            for(i in weak_among_candidate){
              z<-c(z,study_info[[candidate_pos[weak_among_candidate]]]$Coeff^2/study_info[[candidate_pos[weak_among_candidate]]]$Covariance)
            }
            confident_pos_weak<-weak_among_candidate[which.max(z)]
          }
        }else if(length(weak_among_candidate) == 1){
          confident_pos_weak<-candidate_pos[weak_among_candidate]
        }else{
          confident_pos_weak<-c()
        }
        
        true_weak<-c()
        if(length(weak_among_candidate) > 0){
          for( weakone in weak_among_candidate){
            if(max(abs(candidate_cor[,weakone]))<0.1){
              true_weak<-c(true_weak,weakone)      
            }
          }
        }
        if(length(true_weak)>0) true_weak<-candidate_pos[true_weak]
        
        confident_pos_weak <-union(confident_pos_weak,true_weak)
      }else{
        confident_pos_weak <-candidate_pos[weak_among_candidate]
      }
      
      gamma_adaptivelasso<-2
      w_adaptive_candidate<-1/(abs(beta[candidate_pos]))^gamma_adaptivelasso
      w_adaptive_candidate[which(candidate_pos%in%confident_pos_weak)]<-0
      tmp_count<-sapply(1:10, function(j){
        lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 0.9,penalty.factor = w_adaptive_candidate)
        abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
        #coef(lasso_fit_candidate,s='lambda.min')[-1]!=0
      })
      confident_pos<-candidate_pos[which(rowMeans(tmp_count) >= 0.8)]
      
      confident_pos<-union(confident_pos,confident_pos_weak)
      
    }else{
      print("NOT USED")
      confident_pos<-candidate_pos
    }
  }else{
    confident_pos<-candidate_pos
  }
  
  if(filter_index){
    pos<-index_filter[confident_pos]
    pos_bf<-index_filter[candidate_pos]
  }else{
    pos<-confident_pos
    pos_bf<-candidate_pos
  }
  #print(paste0("candidate",paste0(index_filter[candidate_pos],collapse = " ")))
  #print(paste0("confid",paste0(pos,collapse = " ")))
  #print(pos)
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos_bf"=pos_bf,
    "aa_final"=aa_final,
    "candidate"=index_filter[candidate_pos],
    "final_v"=final_v,
    "w_adaptive"=w_adaptive)
}

adaptiveGMMlasso_v1<-function(UKBB_pop,study_info,ld_cut = 0.8,cor_cut=0.8,filter_index=TRUE,filter_by_EAF=FALSE){
  UKBB_cor<-cor(UKBB_pop[,-1])
  diag(UKBB_cor) = 0
  if(filter_index){
    ### Trick 1 index filtering 
    index_filter<-c()
    pval_list<-c()
    for(i in 1:length(study_info)){
      #if(study_info[[i]]$P<p_val_cut){
      index_filter<-c(index_filter,i)
      z_here<-study_info[[i]]$Coeff^2/study_info[[i]]$Covariance
      pval_list<-c(pval_list,z_here)
      #}
    }
    index_filter0<-index_filter
    pval_list0<-pval_list
    index_filter2<-c()
    while (length(index_filter) > 1){
      UKBB_cor_i<-UKBB_cor[index_filter,index_filter]
      cur_i<-which.max(pval_list)
      index_filter2<-c(index_filter2,index_filter[cur_i])
      #print(index_filter[cur_i])
      rm_index<-c(cur_i,which(UKBB_cor_i[cur_i,]>ld_cut))
      #print(index_filter[rm_index])
      index_filter<-index_filter[-rm_index]
      pval_list<-pval_list[-rm_index]
    }
    
    #Nonnull_index2<-which(index_filter2%in%Nonnull_index)
    UKBB_pop<-UKBB_pop[,c(1,(index_filter2+1))]
    study_info<-study_info[index_filter2]
    index_filter<-index_filter2
  }
  
  N_SNP<-ncol(UKBB_pop)-1
  colnames(UKBB_pop)[-1]<-paste0("SNP",1:(N_SNP))
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
  
  #### adaptive lasso with fast lasso computation 
  ## The initial is important
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop*var(UKBB_pop[,1])
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11) 
  
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01)#,penalty.factor = w_adaptive0)
  
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
  diag_theta_sd<-sapply(1:N_SNP,function(i){
    sqrt(study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))
  })
  cov_theta<-t(diag_theta_sd*cor(UKBB_pop[,var_SNP]))*diag_theta_sd/N_Pop
  C_22<-solve(cov_theta)
  C_22_half<-expm::sqrtm(C_22)
  C_half<-adiag(C_11_half,C_22_half)
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  lambda_list<-adaptivelasso_fit$lambda
  UKBB_cor<-cor(UKBB_pop[,-1])
  beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
  #beta<-coef(adaptivelasso_fit,s =lambda_list[30])[-1]
  beta<-coef(adaptivelasso_fit,s =lambda_list[which(adaptivelasso_fit$nzero>1)[which.min(adaptivelasso_fit$cvm[which(adaptivelasso_fit$nzero>1)])]])[-1]
  
  
  while(1){
    index_nonzero_i<-which(beta!=0)
    nonzero_cor_i<-UKBB_cor[index_nonzero_i,index_nonzero_i]
    nonzero_cor_i[nonzero_cor_i == 1] = 0
    tmp<-which(abs(nonzero_cor_i)>cor_cut,arr.ind = T)
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
  
  #W1 = xtx/var(UKBB_pop[,1])
  #W2 = xtx%*%C_22%*%xtx
  #W = W1+W2
  #W_nonzero = W[index_nonzero,index_nonzero]
  #print(beta[index_nonzero])
  #final_v<-diag(inv_Sigsum_scaled_nonzero%*%W_nonzero%*%inv_Sigsum_scaled_nonzero)
  final_v<-diag(inv_Sigsum_scaled_nonzero)
  aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
  
  
  candidate_pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  #candidate_pos<-index_nonzero[which(aa_final<0.05/length(index_nonzero))]
  if(length(candidate_pos) > 1){
    candidate_cor<-cor(UKBB_pop[,c(candidate_pos+1)])
    diag(candidate_cor)<-0
    
    
    if(max(abs(candidate_cor))>0.1){
      
      
      
      UKBB_pop_orig<-UKBB_pop[,c(candidate_pos+1)]
      
      for(i in 1:ncol(UKBB_pop_orig)){
        cur_snp<-UKBB_pop_orig[,i]
        UKBB_pop_orig[which(cur_snp == max(cur_snp)),i] = 2
        UKBB_pop_orig[which(cur_snp == min(cur_snp)),i] = 0
        UKBB_pop_orig[which(cur_snp > min(cur_snp) & cur_snp < max(cur_snp)),i] = 1
      }
      candidate_EAF<-colsums(UKBB_pop_orig)/nrow(UKBB_pop_orig)/2
      
      weak_among_candidate<-which(candidate_EAF < 0.05)
      if(filter_by_EAF){
        if(length(weak_among_candidate) > 1){
          gamma_adaptivelasso<-1/2
          w_adaptive_candidate<-1/(abs(beta[candidate_pos[weak_among_candidate] ]))^gamma_adaptivelasso
          
          
          tmp_count<-sapply(1:10, function(j){
            lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos[weak_among_candidate]+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive_candidate)
            abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
          })
          
          confident_pos_weak<-(candidate_pos[weak_among_candidate])[which(rowMeans(tmp_count) >= 0.8)]
          if(length(confident_pos_weak) == 0){
            z<-c()
            for(i in weak_among_candidate){
              z<-c(z,study_info[[candidate_pos[weak_among_candidate]]]$Coeff^2/study_info[[candidate_pos[weak_among_candidate]]]$Covariance)
            }
            confident_pos_weak<-weak_among_candidate[which.max(z)]
          }
        }else if(length(weak_among_candidate) == 1){
          confident_pos_weak<-candidate_pos[weak_among_candidate]
        }else{
          confident_pos_weak<-c()
        }
        
        true_weak<-c()
        if(length(weak_among_candidate) > 0){
          for( weakone in weak_among_candidate){
            if(max(abs(candidate_cor[,weakone]))<0.1){
              true_weak<-c(true_weak,weakone)      
            }
          }
        }
        if(length(true_weak)>0) true_weak<-candidate_pos[true_weak]
        
        confident_pos_weak <-union(confident_pos_weak,true_weak)
      }else{
        confident_pos_weak <-candidate_pos[weak_among_candidate]
      }
      
      gamma_adaptivelasso<-2
      w_adaptive_candidate<-1/(abs(beta[candidate_pos]))^gamma_adaptivelasso
      w_adaptive_candidate[which(candidate_pos%in%confident_pos_weak)]<-0
      tmp_count<-sapply(1:10, function(j){
        lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 0.9,penalty.factor = w_adaptive_candidate)
        abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
        #coef(lasso_fit_candidate,s='lambda.min')[-1]!=0
      })
      confident_pos<-candidate_pos[which(rowMeans(tmp_count) >= 0.8)]
      
      confident_pos<-union(confident_pos,confident_pos_weak)
      
    }else{
      print("NOT USED")
      confident_pos<-candidate_pos
    }
  }else{
    confident_pos<-candidate_pos
  }
  
  if(filter_index){
    pos<-index_filter[confident_pos]
    pos_bf<-index_filter[candidate_pos]
  }else{
    pos<-confident_pos
    pos_bf<-candidate_pos
  }
  #print(paste0("candidate",paste0(index_filter[candidate_pos],collapse = " ")))
  #print(paste0("confid",paste0(pos,collapse = " ")))
  #print(pos)
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos_bf"=pos_bf,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive)
}

adaptiveGMMlassoC_v1<-function(UKBB_pop,study_info,ld_cut = 0.8,cor_cut=0.8,filter_index=TRUE,filter_by_EAF=FALSE){
  UKBB_cor<-cor(UKBB_pop[,-1])
  diag(UKBB_cor) = 0
  if(filter_index){
    ### Trick 1 index filtering 
    index_filter<-c()
    pval_list<-c()
    for(i in 1:length(study_info)){
      #if(study_info[[i]]$P<p_val_cut){
      index_filter<-c(index_filter,i)
      z_here<-study_info[[i]]$Coeff^2/study_info[[i]]$Covariance
      pval_list<-c(pval_list,z_here)
      #}
    }
    index_filter0<-index_filter
    pval_list0<-pval_list
    index_filter2<-c()
    while (length(index_filter) > 1){
      UKBB_cor_i<-UKBB_cor[index_filter,index_filter]
      cur_i<-which.max(pval_list)
      index_filter2<-c(index_filter2,index_filter[cur_i])
      #print(index_filter[cur_i])
      rm_index<-c(cur_i,which(UKBB_cor_i[cur_i,]>ld_cut))
      #print(index_filter[rm_index])
      index_filter<-index_filter[-rm_index]
      pval_list<-pval_list[-rm_index]
    }
    
    #Nonnull_index2<-which(index_filter2%in%Nonnull_index)
    UKBB_pop<-UKBB_pop[,c(1,(index_filter2+1))]
    study_info<-study_info[index_filter2]
    index_filter<-index_filter2
  }
  
  N_SNP<-ncol(UKBB_pop)-1
  colnames(UKBB_pop)[-1]<-paste0("SNP",1:(N_SNP))
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
  
  
  var_U_beta_theta_func<-function(
    UKBB_pop,
    beta,
    study_info,
    theta_UKBB_GPC=NULL){
    N_Pop<-nrow(UKBB_pop)
    xtbeta<-(UKBB_pop[,-1]%*%beta)
    var_11<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta))
    u2_theta_coef<-sapply(1:N_SNP, function(snp_id){
      if(is.null(var_GPC)){
        xtgamma_id<-c(UKBB_pop[,paste0("SNP",snp_id)]*study_info[[snp_id]]$Coeff)
      }else{
        xtgamma_id<-c((UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]%*%c(theta_UKBB_GPC,study_info[[snp_id]]$Coeff)))
      }
      xtbeta-xtgamma_id
    }) #col is SNP #row is sample 
    var_22<-crossprod(u2_theta_coef*UKBB_pop[,var_SNP],u2_theta_coef*UKBB_pop[,var_SNP])
    var_12<-crossprod(UKBB_pop[,-1]*c(UKBB_pop[,1]-xtbeta),u2_theta_coef*UKBB_pop[,var_SNP])
    (1/N_Pop)*rbind(cbind(var_11,var_12),cbind(t(var_12),var_22))
  }
  
  #### adaptive lasso with fast lasso computation 
  #var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  #var_U1<-crossprod(var_11_half,var_11_half)/N_Pop
  #C_11<-solve(var_U1)
  #C_11_half<-expm::sqrtm(C_11)

  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop*var(UKBB_pop[,1])
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11) 
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01)#,penalty.factor = w_adaptive0)
  
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
  beta_initial<-coef(ridge_fit)[-1]
  
  var_1st_U_beta_theta<-var_U_beta_theta_func(UKBB_pop =UKBB_pop,
    beta = beta_initial,
    study_info = study_info)
  #inv_C_22<-diag(sapply(1:N_SNP,function(i){
  #  (study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2/(N_Pop)
  #}))
  
  diag_theta_sd<-sapply(1:N_SNP,function(i){
    sqrt(study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))
  })
  cov_theta<-t(diag_theta_sd*cor(UKBB_pop[,var_SNP]))*diag_theta_sd/N_Pop
  
  var_2nd_grad_times_theta_hat = adiag(matrix(0,nrow = len_U1,ncol = len_U1),cov_theta)
  
  inv_C = var_1st_U_beta_theta + var_2nd_grad_times_theta_hat
  C_all<-solve(inv_C)
  C_half<-expm::sqrtm(C_all)
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  
  adaptivelasso_fit<-cv.glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive)
  lambda_list<-adaptivelasso_fit$lambda
  UKBB_cor<-cor(UKBB_pop[,-1])
  beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
  #beta<-coef(adaptivelasso_fit,s =lambda_list[which(adaptivelasso_fit$nzero>3)[which.min(adaptivelasso_fit$cvm[which(adaptivelasso_fit$nzero>3)])]])[-1]
  
  
  while(1){
    index_nonzero_i<-which(beta!=0)
    nonzero_cor_i<-UKBB_cor[index_nonzero_i,index_nonzero_i]
    nonzero_cor_i[nonzero_cor_i == 1] = 0
    tmp<-which(abs(nonzero_cor_i)>cor_cut,arr.ind = T)
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
  
  #W1 = xtx/var(UKBB_pop[,1])
  #W2 = xtx%*%C_22%*%xtx
  #W = W1+W2
  #W_nonzero = W[index_nonzero,index_nonzero]
  #print(beta[index_nonzero])
  #final_v<-diag(inv_Sigsum_scaled_nonzero%*%W_nonzero%*%inv_Sigsum_scaled_nonzero)
  final_v<-diag(inv_Sigsum_scaled_nonzero)
  aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
  
  
  candidate_pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  #candidate_pos<-index_nonzero[which(aa_final<0.05/length(index_nonzero))]
  if(length(candidate_pos) > 1){
    candidate_cor<-cor(UKBB_pop[,c(candidate_pos+1)])
    diag(candidate_cor)<-0
    
    
    if(max(abs(candidate_cor))>0.1){
      
      
      
      UKBB_pop_orig<-UKBB_pop[,c(candidate_pos+1)]
      
      for(i in 1:ncol(UKBB_pop_orig)){
        cur_snp<-UKBB_pop_orig[,i]
        UKBB_pop_orig[which(cur_snp == max(cur_snp)),i] = 2
        UKBB_pop_orig[which(cur_snp == min(cur_snp)),i] = 0
        UKBB_pop_orig[which(cur_snp > min(cur_snp) & cur_snp < max(cur_snp)),i] = 1
      }
      candidate_EAF<-colsums(UKBB_pop_orig)/nrow(UKBB_pop_orig)/2
      
      weak_among_candidate<-which(candidate_EAF < 0.05)
      if(filter_by_EAF){
        if(length(weak_among_candidate) > 1){
          gamma_adaptivelasso<-1/2
          w_adaptive_candidate<-1/(abs(beta[candidate_pos[weak_among_candidate] ]))^gamma_adaptivelasso
          
          
          tmp_count<-sapply(1:10, function(j){
            lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos[weak_among_candidate]+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive_candidate)
            abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
          })
          
          confident_pos_weak<-(candidate_pos[weak_among_candidate])[which(rowMeans(tmp_count) >= 0.8)]
          if(length(confident_pos_weak) == 0){
            z<-c()
            for(i in weak_among_candidate){
              z<-c(z,study_info[[candidate_pos[weak_among_candidate]]]$Coeff^2/study_info[[candidate_pos[weak_among_candidate]]]$Covariance)
            }
            confident_pos_weak<-weak_among_candidate[which.max(z)]
          }
        }else if(length(weak_among_candidate) == 1){
          confident_pos_weak<-candidate_pos[weak_among_candidate]
        }else{
          confident_pos_weak<-c()
        }
        
        true_weak<-c()
        if(length(weak_among_candidate) > 0){
          for( weakone in weak_among_candidate){
            if(max(abs(candidate_cor[,weakone]))<0.1){
              true_weak<-c(true_weak,weakone)      
            }
          }
        }
        if(length(true_weak)>0) true_weak<-candidate_pos[true_weak]
        
        confident_pos_weak <-union(confident_pos_weak,true_weak)
      }else{
        confident_pos_weak <-candidate_pos[weak_among_candidate]
      }
      
      gamma_adaptivelasso<-2
      w_adaptive_candidate<-1/(abs(beta[candidate_pos]))^gamma_adaptivelasso
      w_adaptive_candidate[which(candidate_pos%in%confident_pos_weak)]<-0
      tmp_count<-sapply(1:10, function(j){
        lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 0.9,penalty.factor = w_adaptive_candidate)
        abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
        #coef(lasso_fit_candidate,s='lambda.min')[-1]!=0
      })
      confident_pos<-candidate_pos[which(rowMeans(tmp_count) >= 0.8)]
      
      confident_pos<-union(confident_pos,confident_pos_weak)
      
    }else{
      print("NOT USED")
      confident_pos<-candidate_pos
    }
  }else{
    confident_pos<-candidate_pos
  }
  
  if(filter_index){
    pos<-index_filter[confident_pos]
    pos_bf<-index_filter[candidate_pos]
  }else{
    pos<-confident_pos
    pos_bf<-candidate_pos
  }
  #print(paste0("candidate",paste0(index_filter[candidate_pos],collapse = " ")))
  #print(paste0("confid",paste0(pos,collapse = " ")))
  #print(pos)
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos_bf"=pos_bf,
    "aa_final"=aa_final,
    "candidate"=index_filter[candidate_pos],
    "final_v"=final_v,
    "w_adaptive"=w_adaptive)
}

adaptiveGMMlasso35BMI<-function(UKBB_pop,BMI,study_info,ld_cut = 0.8,cor_cut=0.75,filter_index=TRUE,filter_by_EAF=FALSE){
  UKBB_cor<-cor(UKBB_pop[,-1])
  diag(UKBB_cor) = 0
  if(filter_index){
    ### Trick 1 index filtering 
    index_filter<-c()
    pval_list<-c()
    for(i in 1:length(study_info)){
      #if(study_info[[i]]$P<p_val_cut){
      index_filter<-c(index_filter,i)
      z_here<-study_info[[i]]$Coeff^2/study_info[[i]]$Covariance
      pval_list<-c(pval_list,z_here)
      #}
    }
    index_filter0<-index_filter
    pval_list0<-pval_list
    index_filter2<-c()
    while (length(index_filter) > 1){
      UKBB_cor_i<-UKBB_cor[index_filter,index_filter]
      cur_i<-which.max(pval_list)
      index_filter2<-c(index_filter2,index_filter[cur_i])
      #print(index_filter[cur_i])
      rm_index<-c(cur_i,which(UKBB_cor_i[cur_i,]>ld_cut))
      #print(index_filter[rm_index])
      index_filter<-index_filter[-rm_index]
      pval_list<-pval_list[-rm_index]
    }
    
    #Nonnull_index2<-which(index_filter2%in%Nonnull_index)
    UKBB_pop<-UKBB_pop[,c(1,(index_filter2+1))]
    study_info<-study_info[index_filter2]
    index_filter<-index_filter2
  }
  
  N_SNP<-ncol(UKBB_pop)-1
  colnames(UKBB_pop)[-1]<-paste0("SNP",1:(N_SNP))
  colnames(UKBB_pop)[1]<-"Y"
  var_SNP<-paste0("SNP",1:(N_SNP))
  len_SNP<-length(var_SNP)
  var_GPC<-NULL
  len_GPC<-length(var_GPC)
  var_nonSNP<-c(var_GPC,"BMI")
  len_nonSNP<-length(var_nonSNP)
  var_names_full_fit <- c(var_SNP,var_nonSNP)
  colname_UKBB<-colnames(UKBB_pop)
  len_U1 = len_SNP+len_nonSNP
  len_U2 = len_SNP
  len_beta = len_U1
  len_theta = len_SNP+len_GPC
  len_U = len_U1 + len_U2
  N_Pop<-nrow(UKBB_pop)
  UKBB_pop<-cbind(UKBB_pop,BMI)
  colnames(UKBB_pop)[ncol(UKBB_pop)]<-"BMI"
  ############## Beta
  pseudo_Xy<-function(
    C_half,UKBB_pop,
    var_SNP,var_GPC,
    N_SNP,theta_UKBB_GPC,study_info){
    N_Pop<-nrow(UKBB_pop)
    pseudo_X<-C_half%*%rbind(t(UKBB_pop[,-1]),t(UKBB_pop[,var_SNP]))%*%UKBB_pop[,-1]
    pseudo_y1<-crossprod(UKBB_pop[,1],UKBB_pop[,-1])
    pseudo_y22<-sapply(1:N_SNP, function(snp_id){
      u2_id<-UKBB_pop[,c(var_GPC,paste0("SNP",snp_id))]*(study_info[[snp_id]]$Coeff)
      c(c(u2_id)%*%UKBB_pop[,paste0("SNP",snp_id)])
    })
    
    pseudo_y<-c(c(c(pseudo_y1),c(pseudo_y22))%*%C_half)
    newList<-list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
    newList
  }
  
  #### adaptive lasso with fast lasso computation 
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta0)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop*var(UKBB_pop[,1]) 
  ## Seems to be correlation matrix.
  ## Have taken the scaling 
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11)
  var_U2<-diag(sapply(1:N_SNP,function(i){
    (study_info[[i]]$Covariance)*(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2/(N_Pop)
  }))
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  non_penalize_BMI<-c(rep(1,N_SNP),0)
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01,penalty.factor = non_penalize_BMI)
  
  gamma_adaptivelasso<-1/2
  w_adaptive<-c(1/(abs(coef(ridge_fit)[-1])[1:N_SNP])^gamma_adaptivelasso,0)
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
  
  
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
  UKBB_cor<-cor(UKBB_pop[,var_SNP])
  beta<-coef(adaptivelasso_fit,s ='lambda.min')[-1]
  
  
  while(1){
    index_nonzero_i<-which(beta[1:N_SNP]!=0)
    nonzero_cor_i<-UKBB_cor[index_nonzero_i,index_nonzero_i]
    nonzero_cor_i[nonzero_cor_i == 1] = 0
    tmp<-which(abs(nonzero_cor_i)>cor_cut,arr.ind = T)
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
  
  index_nonzero<-c(which(beta[1:N_SNP]!=0),length(beta))
  xtx<-t(UKBB_pop[,-1])%*%UKBB_pop[,-1]/N_Pop
  xtx2<-t(UKBB_pop[,-1])%*%UKBB_pop[,var_SNP]/N_Pop
  Sigsum_half<-cbind(xtx,xtx2)%*%C_half
  Sigsum_scaled<-Sigsum_half%*%t(Sigsum_half)
  Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero]
  inv_Sigsum_scaled_nonzero<-solve(Sigsum_scaled_nonzero)
  
  #W1 = xtx/var(UKBB_pop[,1])
  #W2 = xtx%*%C_22%*%xtx
  #W = W1+W2
  #W_nonzero = W[index_nonzero,index_nonzero]
  #print(beta[index_nonzero])
  #final_v<-diag(inv_Sigsum_scaled_nonzero%*%W_nonzero%*%inv_Sigsum_scaled_nonzero)
  index_nonzero<-index_nonzero[-length(index_nonzero)]
  final_v<-diag(inv_Sigsum_scaled_nonzero)
  final_v<-final_v[-length(final_v)]
  aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
  
  
  #candidate_pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  candidate_pos<-index_nonzero[which(aa_final<0.05/length(index_nonzero))]
  if(length(candidate_pos) > 1){
  candidate_cor<-cor(UKBB_pop[,c(candidate_pos+1)])
  diag(candidate_cor)<-0
  
  
  if(max(abs(candidate_cor))>0.1){
    
    UKBB_pop_orig<-UKBB_pop[,c(candidate_pos+1)]
    
    for(i in 1:ncol(UKBB_pop_orig)){
      cur_snp<-UKBB_pop_orig[,i]
      UKBB_pop_orig[which(cur_snp == max(cur_snp)),i] = 2
      UKBB_pop_orig[which(cur_snp == min(cur_snp)),i] = 0
      UKBB_pop_orig[which(cur_snp > min(cur_snp) & cur_snp < max(cur_snp)),i] = 1
    }
    candidate_EAF<-colsums(UKBB_pop_orig)/nrow(UKBB_pop_orig)/2
    
    weak_among_candidate<-which(candidate_EAF < 0.05)
    if(filter_by_EAF){
      if(length(weak_among_candidate) > 1){
        gamma_adaptivelasso<-1/2
        w_adaptive_candidate<-c(1/(abs(beta[candidate_pos[weak_among_candidate] ]))^gamma_adaptivelasso,0)
        
        tmp_count<-sapply(1:10, function(j){
          lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,c((candidate_pos[weak_among_candidate]+1),ncol(UKBB_pop)) ],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive_candidate)
          tmptmp<-abs(coef(lasso_fit_candidate,s='lambda.min')[-1])
          tmptmp<-tmptmp[-length(tmptmp)]
          tmptmp>1e-4
        })
        
        confident_pos_weak<-(candidate_pos[weak_among_candidate])[which(rowMeans(tmp_count) >= 0.8)]
        if(length(confident_pos_weak) == 0){
          z<-c()
          for(i in weak_among_candidate){
            z<-c(z,study_info[[candidate_pos[weak_among_candidate]]]$Coeff^2/study_info[[candidate_pos[weak_among_candidate]]]$Covariance)
          }
          confident_pos_weak<-weak_among_candidate[which.max(z)]
        }
      }else if(length(weak_among_candidate) == 1){
        confident_pos_weak<-candidate_pos[weak_among_candidate]
      }else{
        confident_pos_weak<-c()
      }
      
      true_weak<-c()
      if(length(weak_among_candidate) > 0){
        for( weakone in weak_among_candidate){
          if(max(abs(candidate_cor[,weakone]))<0.1){
            true_weak<-c(true_weak,weakone)      
          }
        }
      }
      if(length(true_weak)>0) true_weak<-candidate_pos[true_weak]
      
      confident_pos_weak <-union(confident_pos_weak,true_weak)
    }else{
      confident_pos_weak <-candidate_pos[weak_among_candidate]
    }
    
    gamma_adaptivelasso<-2
    w_adaptive_candidate<-c(1/(abs(beta[candidate_pos]))^gamma_adaptivelasso,0)
    w_adaptive_candidate[which(candidate_pos%in%confident_pos_weak)]<-0
    tmp_count<-sapply(1:10, function(j){
      lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,c( (candidate_pos+1),ncol(UKBB_pop)) ] ,y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 0.9,penalty.factor = w_adaptive_candidate)
      tmptmp<-abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
      tmptmp<-tmptmp[-length(tmptmp)]
      tmptmp>1e-4
      #coef(lasso_fit_candidate,s='lambda.min')[-1]!=0
    })
    confident_pos<-candidate_pos[which(rowMeans(tmp_count) >= 0.8)]
    
    confident_pos<-union(confident_pos,confident_pos_weak)
    
  }else{
    print("NOT USED")
    confident_pos<-candidate_pos
  }
  }else{
    confident_pos<-candidate_pos
  }
  
  if(filter_index){
    pos2<-index_filter[confident_pos]
    pos<-index_filter[confident_pos]
  }else{
    pos<-confident_pos
    pos2<-confident_pos
  }
  print(paste0("candidate",paste0(index_filter[candidate_pos],collapse = " ")))
  print(paste0("confid",paste0(pos,collapse = " ")))
  #print(pos)
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos2"=pos2,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive)
}

adaptiveGMMlasso35adhoc<-function(UKBB_pop,study_info,ld_cut = 0.8,cor_cut=0.75,filter_index=TRUE,filter_by_EAF=FALSE){
  UKBB_cor<-cor(UKBB_pop[,-1])
  diag(UKBB_cor) = 0
  if(filter_index){
    ### Trick 1 index filtering 
    index_filter<-c()
    pval_list<-c()
    for(i in 1:length(study_info)){
      #if(study_info[[i]]$P<p_val_cut){
      index_filter<-c(index_filter,i)
      z_here<-study_info[[i]]$Coeff^2/study_info[[i]]$Covariance
      pval_list<-c(pval_list,z_here)
      #}
    }
    index_filter0<-index_filter
    pval_list0<-pval_list
    index_filter2<-c()
    while (length(index_filter) > 1){
      UKBB_cor_i<-UKBB_cor[index_filter,index_filter]
      cur_i<-which.max(pval_list)
      index_filter2<-c(index_filter2,index_filter[cur_i])
      #print(index_filter[cur_i])
      rm_index<-c(cur_i,which(UKBB_cor_i[cur_i,]>ld_cut))
      #print(index_filter[rm_index])
      index_filter<-index_filter[-rm_index]
      pval_list<-pval_list[-rm_index]
    }
    
    #Nonnull_index2<-which(index_filter2%in%Nonnull_index)
    UKBB_pop<-UKBB_pop[,c(1,(index_filter2+1))]
    study_info<-study_info[index_filter2]
    index_filter<-index_filter2
  }
  
  N_SNP<-ncol(UKBB_pop)-1
  colnames(UKBB_pop)[-1]<-paste0("SNP",1:(N_SNP))
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
  
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01)#,penalty.factor = w_adaptive0)
  
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
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
    tmp<-which(abs(nonzero_cor_i)>cor_cut,arr.ind = T)
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
  
  #W1 = xtx/var(UKBB_pop[,1])
  #W2 = xtx%*%C_22%*%xtx
  #W = W1+W2
  #W_nonzero = W[index_nonzero,index_nonzero]
  #print(beta[index_nonzero])
  #final_v<-diag(inv_Sigsum_scaled_nonzero%*%W_nonzero%*%inv_Sigsum_scaled_nonzero)
  final_v<-diag(inv_Sigsum_scaled_nonzero)
  aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
  
  
  candidate_pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  #candidate_pos<-index_nonzero[which(aa_final<0.05/length(index_nonzero))]
  if(length(candidate_pos) > 1){
    candidate_cor<-cor(UKBB_pop[,c(candidate_pos+1)])
    diag(candidate_cor)<-0
    
    
    if(max(abs(candidate_cor))>0.1){
      
      
      
      UKBB_pop_orig<-UKBB_pop[,c(candidate_pos+1)]
      
      for(i in 1:ncol(UKBB_pop_orig)){
        cur_snp<-UKBB_pop_orig[,i]
        UKBB_pop_orig[which(cur_snp == max(cur_snp)),i] = 2
        UKBB_pop_orig[which(cur_snp == min(cur_snp)),i] = 0
        UKBB_pop_orig[which(cur_snp > min(cur_snp) & cur_snp < max(cur_snp)),i] = 1
      }
      candidate_EAF<-colsums(UKBB_pop_orig)/nrow(UKBB_pop_orig)/2
      
      weak_among_candidate<-which(candidate_EAF < 0.05)
      if(filter_by_EAF){
        if(length(weak_among_candidate) > 1){
          gamma_adaptivelasso<-1/2
          w_adaptive_candidate<-1/(abs(beta[candidate_pos[weak_among_candidate] ]))^gamma_adaptivelasso
          
          
          tmp_count<-sapply(1:10, function(j){
            lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos[weak_among_candidate]+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive_candidate)
            abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
          })
          
          confident_pos_weak<-(candidate_pos[weak_among_candidate])[which(rowMeans(tmp_count) >= 0.8)]
          if(length(confident_pos_weak) == 0){
            z<-c()
            for(i in weak_among_candidate){
              z<-c(z,study_info[[candidate_pos[weak_among_candidate]]]$Coeff^2/study_info[[candidate_pos[weak_among_candidate]]]$Covariance)
            }
            confident_pos_weak<-weak_among_candidate[which.max(z)]
          }
        }else if(length(weak_among_candidate) == 1){
          confident_pos_weak<-candidate_pos[weak_among_candidate]
        }else{
          confident_pos_weak<-c()
        }
        
        true_weak<-c()
        if(length(weak_among_candidate) > 0){
          for( weakone in weak_among_candidate){
            if(max(abs(candidate_cor[,weakone]))<0.1){
              true_weak<-c(true_weak,weakone)      
            }
          }
        }
        if(length(true_weak)>0) true_weak<-candidate_pos[true_weak]
        
        confident_pos_weak <-union(confident_pos_weak,true_weak)
      }else{
        confident_pos_weak <-candidate_pos[weak_among_candidate]
      }
      
      gamma_adaptivelasso<-2
      w_adaptive_candidate<-1/(abs(beta[candidate_pos]))^gamma_adaptivelasso
      w_adaptive_candidate[which(candidate_pos%in%confident_pos_weak)]<-0
      tmp_count<-sapply(1:10, function(j){
        lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 0.9,penalty.factor = w_adaptive_candidate)
        abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
        #coef(lasso_fit_candidate,s='lambda.min')[-1]!=0
      })
      confident_pos<-candidate_pos[which(rowMeans(tmp_count) >= 0.8)]
      
      confident_pos<-union(confident_pos,confident_pos_weak)
      
    }else{
      print("NOT USED")
      confident_pos<-candidate_pos
    }
  }else{
    confident_pos<-candidate_pos
  }
  
  if(filter_index){
    pos2<-index_filter[confident_pos]
    pos<-index_filter[confident_pos]
  }else{
    pos<-confident_pos
    pos2<-confident_pos
  }
  print(paste0("candidate",paste0(index_filter[candidate_pos],collapse = " ")))
  print(paste0("confid",paste0(pos,collapse = " ")))
  #print(pos)
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos2"=pos2,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive)
}
## C_22 wrong
adaptiveGMMlasso35<-function(UKBB_pop,study_info,ld_cut = 0.8,cor_cut=0.75,filter_index=TRUE,filter_by_EAF=FALSE){
  UKBB_cor<-cor(UKBB_pop[,-1])
  diag(UKBB_cor) = 0
  if(filter_index){
    ### Trick 1 index filtering 
    index_filter<-c()
    pval_list<-c()
    for(i in 1:length(study_info)){
      #if(study_info[[i]]$P<p_val_cut){
      index_filter<-c(index_filter,i)
      z_here<-study_info[[i]]$Coeff^2/study_info[[i]]$Covariance
      pval_list<-c(pval_list,z_here)
      #}
    }
    index_filter0<-index_filter
    pval_list0<-pval_list
    index_filter2<-c()
    while (length(index_filter) > 1){
      UKBB_cor_i<-UKBB_cor[index_filter,index_filter]
      cur_i<-which.max(pval_list)
      index_filter2<-c(index_filter2,index_filter[cur_i])
      #print(index_filter[cur_i])
      rm_index<-c(cur_i,which(UKBB_cor_i[cur_i,]>ld_cut))
      #print(index_filter[rm_index])
      index_filter<-index_filter[-rm_index]
      pval_list<-pval_list[-rm_index]
    }
    
    #Nonnull_index2<-which(index_filter2%in%Nonnull_index)
    UKBB_pop<-UKBB_pop[,c(1,(index_filter2+1))]
    study_info<-study_info[index_filter2]
    index_filter<-index_filter2
  }
  
  N_SNP<-ncol(UKBB_pop)-1
  colnames(UKBB_pop)[-1]<-paste0("SNP",1:(N_SNP))
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
  
  #### adaptive lasso with fast lasso computation 
  var_11_half<-UKBB_pop[,-1]#*c(UKBB_pop[,1]-UKBB_pop[,-1]%*%beta)
  var_U1<-crossprod(var_11_half,var_11_half)/N_Pop*var(UKBB_pop[,1])
  C_11<-solve(var_U1)
  C_11_half<-expm::sqrtm(C_11)
  C_22<-diag(sapply(1:N_SNP,function(i){
    (N_Pop)/(study_info[[i]]$Covariance)/(c(UKBB_pop[,c(paste0("SNP",i))]%*%UKBB_pop[,c(paste0("SNP",i))]))^2
  }))
  C_22<-C_22*study_info[[1]]$Sample_size/N_Pop
  C_22_half<-diag(sqrt(diag(C_22)))
  C_half<-adiag(C_11_half,C_22_half)
  
  
  pseudo_Xy_list<-pseudo_Xy(C_half,UKBB_pop,var_SNP,var_GPC,N_SNP,theta_UKBB_GPC,study_info)
  initial_sf<-N_Pop/sqrt(nrow(pseudo_Xy_list$pseudo_X))
  pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
  pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
  
  lasso_initial<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F)
  lambda_list0<-lasso_initial$lambda
  ridge_fit<-glmnet(x= (pseudo_X),y= (pseudo_y),standardize=F,intercept=F,lambda = lambda_list0[50]/100,alpha = 0.01)#,penalty.factor = w_adaptive0)
  
  gamma_adaptivelasso<-1/2
  w_adaptive<-1/(abs(coef(ridge_fit)[-1]))^gamma_adaptivelasso
  w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
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
    tmp<-which(abs(nonzero_cor_i)>cor_cut,arr.ind = T)
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
  
  #W1 = xtx/var(UKBB_pop[,1])
  #W2 = xtx%*%C_22%*%xtx
  #W = W1+W2
  #W_nonzero = W[index_nonzero,index_nonzero]
  #print(beta[index_nonzero])
  #final_v<-diag(inv_Sigsum_scaled_nonzero%*%W_nonzero%*%inv_Sigsum_scaled_nonzero)
  final_v<-diag(inv_Sigsum_scaled_nonzero)
  aa_final<-pchisq(N_Pop*beta[index_nonzero]^2/final_v,1,lower.tail = F)
  
  
  candidate_pos<-index_nonzero[which(aa_final<0.05/ncol(UKBB_pop))]
  #candidate_pos<-index_nonzero[which(aa_final<0.05/length(index_nonzero))]
  if(length(candidate_pos) > 1){
    candidate_cor<-cor(UKBB_pop[,c(candidate_pos+1)])
    diag(candidate_cor)<-0
    
    
    if(max(abs(candidate_cor))>0.1){
      
      
      
      UKBB_pop_orig<-UKBB_pop[,c(candidate_pos+1)]
      
      for(i in 1:ncol(UKBB_pop_orig)){
        cur_snp<-UKBB_pop_orig[,i]
        UKBB_pop_orig[which(cur_snp == max(cur_snp)),i] = 2
        UKBB_pop_orig[which(cur_snp == min(cur_snp)),i] = 0
        UKBB_pop_orig[which(cur_snp > min(cur_snp) & cur_snp < max(cur_snp)),i] = 1
      }
      candidate_EAF<-colsums(UKBB_pop_orig)/nrow(UKBB_pop_orig)/2
      
      weak_among_candidate<-which(candidate_EAF < 0.05)
      if(filter_by_EAF){
        if(length(weak_among_candidate) > 1){
          gamma_adaptivelasso<-1/2
          w_adaptive_candidate<-1/(abs(beta[candidate_pos[weak_among_candidate] ]))^gamma_adaptivelasso
          
          
          tmp_count<-sapply(1:10, function(j){
            lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos[weak_among_candidate]+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 1,penalty.factor = w_adaptive_candidate)
            abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
          })
          
          confident_pos_weak<-(candidate_pos[weak_among_candidate])[which(rowMeans(tmp_count) >= 0.8)]
          if(length(confident_pos_weak) == 0){
            z<-c()
            for(i in weak_among_candidate){
              z<-c(z,study_info[[candidate_pos[weak_among_candidate]]]$Coeff^2/study_info[[candidate_pos[weak_among_candidate]]]$Covariance)
            }
            confident_pos_weak<-weak_among_candidate[which.max(z)]
          }
        }else if(length(weak_among_candidate) == 1){
          confident_pos_weak<-candidate_pos[weak_among_candidate]
        }else{
          confident_pos_weak<-c()
        }
        
        true_weak<-c()
        if(length(weak_among_candidate) > 0){
          for( weakone in weak_among_candidate){
            if(max(abs(candidate_cor[,weakone]))<0.1){
              true_weak<-c(true_weak,weakone)      
            }
          }
        }
        if(length(true_weak)>0) true_weak<-candidate_pos[true_weak]
        
        confident_pos_weak <-union(confident_pos_weak,true_weak)
      }else{
        confident_pos_weak <-candidate_pos[weak_among_candidate]
      }
      
      gamma_adaptivelasso<-2
      w_adaptive_candidate<-1/(abs(beta[candidate_pos]))^gamma_adaptivelasso
      w_adaptive_candidate[which(candidate_pos%in%confident_pos_weak)]<-0
      tmp_count<-sapply(1:10, function(j){
        lasso_fit_candidate<-cv.glmnet(x= UKBB_pop[,(candidate_pos+1)],y= UKBB_pop[,1],standardize=F,intercept=F,alpha = 0.9,penalty.factor = w_adaptive_candidate)
        abs(coef(lasso_fit_candidate,s='lambda.min')[-1])>1e-4
        #coef(lasso_fit_candidate,s='lambda.min')[-1]!=0
      })
      confident_pos<-candidate_pos[which(rowMeans(tmp_count) >= 0.8)]
      
      confident_pos<-union(confident_pos,confident_pos_weak)
      
    }else{
      print("NOT USED")
      confident_pos<-candidate_pos
    }
  }else{
    confident_pos<-candidate_pos
  }
  
  if(filter_index){
    pos2<-index_filter[confident_pos]
    pos<-index_filter[confident_pos]
  }else{
    pos<-confident_pos
    pos2<-confident_pos
  }
  print(paste0("candidate",paste0(index_filter[candidate_pos],collapse = " ")))
  print(paste0("confid",paste0(pos,collapse = " ")))
  #print(pos)
  newList<-list("beta"=beta,
    "pos"=pos,
    "pos2"=pos2,
    "aa_final"=aa_final,
    "final_v"=final_v,
    "w_adaptive"=w_adaptive)
}





real_data_func<-function(region_name="T2D_FTO_v3"){
  
  library(inline)
  library(data.table)
  library(dplyr)
  library(speedglm)
  library(Rfast)
  library(feather)
  library(locfit)
  library(bigmemory)
  library(stringr)
  plus_number<-sample(1:1e7,1)
  library(MultiPhen)
  #region_name<-"T2D_FTO"
  
  library(pracma)
  library(Matrix)
  library(genio)
  library(glmnet)
  library(inline)
  library(data.table)
  library(dplyr)
  library(speedglm)
  library(Rfast)
  library(feather)
  library(locfit)
  library(bigmemory)
  library(stringr)
  library(MultiPhen)
  
  #foldpath<-"/dcl01/chatterj/data/rzhao/T2D_UKBB/"
  #region_name<-"T2D_FTO2"
  foldpath<-"/users/rzhao1/T2D_UKBB/"
  
  # load whole SNP data
  .<-capture.output(ref<-read.plink(paste0(foldpath,region_name)))
  fam <-read.table(paste0(foldpath,region_name,".fam"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  bim <-read.table(paste0(foldpath,region_name,".bim"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  colnames(bim)<-c("chr", "id", "posg", "pos", "ref", "alt")
  colnames(fam)<-c("fam", "id", "pat", "mat", "sex", "pheno")
  toomuchna<-rowSums(is.na(ref))
  na_index<-which(toomuchna>ncol(ref)/20)
  ref<-ref[-na_index,]
  fam<-fam[-na_index,]
  duplicated_index<-which(duplicated(t(ref))==TRUE)
  ref<-ref[,-duplicated_index]
  bim<-bim[-duplicated_index,]
  dim(ref)
  na_index_snp<-which(colSums(is.na(ref))/nrow(ref)>0.01)
  ref<-ref[,-na_index_snp]
  bim<-bim[-na_index_snp,]
  ref[is.na(ref)]<-0
  dim(ref)
  foldpath1<-"~/T2D_UKBB_data/"
  summary_statistics<-readRDS(paste0(foldpath1,"summary_stat_Chr16.rds"))
  keep_index<-which(paste0("16:",bim$pos)%in%summary_statistics$SNP[which(summary_statistics$Chr=="16")])
  length(keep_index)
  
  ref<-ref[,keep_index]
  bim<-bim[keep_index,]
  which(bim$id == "rs1121980")
  
  
  foldpath2<-"/dcl01/chatterj/data/rzhao/T2D_UKBB/"
  pheno_covariates<-read.csv(paste0(foldpath2,"t2d_UKBB.csv"))
  index_sample<-which(fam$fam%in%pheno_covariates$eid)
  fam<-fam[index_sample,]
  ref<-ref[index_sample,]
  dim(ref)
  
  
  pheno_covariates<-pheno_covariates[which(pheno_covariates$eid %in% fam$fam),]
  rownames(pheno_covariates)<-pheno_covariates$eid
  pheno_covariates<-pheno_covariates[as.character(fam$fam),]
  ref<-ref[as.character(fam$fam),]
  sum(rownames(ref) != as.character(fam$fam))
  sum(rownames(pheno_covariates)!=as.character(fam$fam))
  pheno_covariates$eid<-as.numeric(pheno_covariates$eid)
  
  N_SNP<-ncol(ref)
  
  ########################
  
  colnames(ref)<-bim$id
  foldpath<-"/dcs04/nilanjan/data/rzhao/T2D_UKBB/"
  cur_iter<-0
  .<-capture.output(write_plink(paste0(foldpath,"ukbb_pop_",cur_iter),as.matrix(t(ref)),
    bim =  bim,fam = fam))
  gwas_res0<-summary_statistics[which(summary_statistics$Chr=="16")[which(summary_statistics$SNP[which(summary_statistics$Chr=="16")]%in%paste0("16:",bim$pos))],]
  sum(paste0("16:",bim$pos)%in%gwas_res0$SNP)
  ma_file<- data.frame(SNP=bim$id,
    A1=bim$ref,A2=bim$alt,
    freq=gwas_res0$EAF,b=gwas_res0$Beta,
    se=gwas_res0$SE,p=gwas_res0$Pvalue,
    N=gwas_res0$Neff[1])
  
  readr::write_tsv(ma_file,paste0(foldpath,"ukbb_pop_",cur_iter,".ma"))
  
  aaa<-system(paste0("gcta64  --bfile ",foldpath,"ukbb_pop_",cur_iter," --cojo-file ",foldpath,"ukbb_pop_", cur_iter,".ma --cojo-slct --cojo-p 1e-5 --out ",foldpath,"test_",cur_iter),intern = TRUE)
  
  cojo_jma<-read.table(paste0(foldpath,"test_",cur_iter,".jma.cojo"),header = T)    
  
  aft_GWAS_jma<-cojo_jma$SNP[cojo_jma$pJ<0.05/N_SNP]
  predict_nonnull_index<-which(bim$id%in%aft_GWAS_jma)
  print(paste0("COJOJMA: length: ",length(aft_GWAS_jma)))
  file.remove(paste0("test_",cur_iter,".cma.cojo"))
  file.remove(paste0("test_",cur_iter,".jma.cojo"))
  file.remove(paste0("test_",cur_iter,".ldr.cojo"))
  file.remove(paste0("test_",cur_iter,".log"))
  file.remove(paste0("ukbb_pop_",cur_iter,".ma"))
  file.remove(paste0("ukbb_pop_",cur_iter,".bim"))
  file.remove(paste0("ukbb_pop_",cur_iter,".bed"))
  file.remove(paste0("ukbb_pop_",cur_iter,".fam"))
  
  ###### Run GWAS on this data 
  
  if(1){
    ####### My own summary statistics 
    Phenotype=pheno_covariates$Y_t2d
    cur_iter<-0
    EUR<-fam[,1:2]
    write.table(EUR, paste0("/dcs04/nilanjan/data/rzhao/T2D_UKBB/sim_GWAS/EUR",cur_iter,".eid"),row.names =F,col.names = F)
    pheno_EUR <- data.frame(FID=fam[,1], IID=fam[,2],T2D=Phenotype+1)
    readr::write_tsv(pheno_EUR, paste0("/dcs04/nilanjan/data/rzhao/T2D_UKBB/sim_GWAS/pheno_EUR",cur_iter,".txt"))
    EUR_SNP<-data.frame(SNP=bim$id)
    readr::write_tsv(EUR_SNP, paste0("/dcs04/nilanjan/data/rzhao/T2D_UKBB/sim_GWAS/EUR_SNP",cur_iter,".txt"),col_names = F)
    
    
    colnames(ref)<-bim$id
    
    colnames(fam)<-c("fam", "id", "pat", "mat", "sex", "pheno")
    library(genio)
    write_plink(paste0("/dcs04/nilanjan/data/rzhao/T2D_UKBB/sim_GWAS/ukbb_pop_all",cur_iter),as.matrix(t(ref)),
      bim =  bim,fam = fam)
    
    arg= paste0("/dcl01/chatterj/data/tools/plink2  --extract /dcs04/nilanjan/data/rzhao/T2D_UKBB/sim_GWAS/EUR_SNP",cur_iter,".txt --keep /dcs04/nilanjan/data/rzhao/T2D_UKBB/sim_GWAS/EUR",cur_iter,".eid  --bfile /dcs04/nilanjan/data/rzhao/T2D_UKBB/sim_GWAS/ukbb_pop_all",cur_iter," --snps-only --pheno /dcs04/nilanjan/data/rzhao/T2D_UKBB/sim_GWAS/pheno_EUR",cur_iter,".txt --glm hide-covar --vif 10000 --covar-variance-standardize --out /dcs04/nilanjan/data/rzhao/T2D_UKBB/sim_GWAS/GWAS_T2D_FTO_",cur_iter)
    system(arg,ignore.stdout = T)
    
    
    gwas_res<-read.table(paste0("/dcs04/nilanjan/data/rzhao/T2D_UKBB/sim_GWAS/GWAS_T2D_FTO_",cur_iter,".T2D.glm.logistic"))
    colnames(gwas_res)<-c("CHROM","POS","ID","REF","ALT","TEST","OBS_CT","BETA","SE","T_STAT","P")
    
    gwas_res$BETA<-log(gwas_res$BETA)
    
    
    sim_GWAS_path<-"/dcs04/nilanjan/data/rzhao/T2D_UKBB/sim_GWAS/"
    .<-capture.output(file.remove(paste0(sim_GWAS_path,"EUR",cur_iter,".eid")))
    .<-capture.output(file.remove(paste0(sim_GWAS_path,"pheno_EUR",cur_iter,".txt")))
    .<-capture.output(file.remove(paste0(sim_GWAS_path,"EUR_SNP",cur_iter,".txt")))
    .<-capture.output(file.remove(paste0(sim_GWAS_path,"ukbb_pop_all",cur_iter,".bed")))
    .<-capture.output(file.remove(paste0(sim_GWAS_path,"ukbb_pop_all",cur_iter,".fam")))
    .<-capture.output(file.remove(paste0(sim_GWAS_path,"ukbb_pop_all",cur_iter,".bim")))
    .<-capture.output(file.remove(paste0(sim_GWAS_path,"GWAS_T2D_FTO_",cur_iter,".T2D.glm.logistic")))
    .<-capture.output(file.remove(paste0(sim_GWAS_path,"GWAS_T2D_FTO_",cur_iter,".T2D.glm.logistic.id")))
    .<-capture.output(file.remove(paste0(sim_GWAS_path,"GWAS_T2D_FTO_",cur_iter,".log")))
    
    
    gwas_res_P_meta<-gwas_res$P
    gwas_res_BETA_meta<-gwas_res$BETA
    gwas_res_SE_meta<-gwas_res$SE
    
    for( i in 1:N_SNP){
      w1<-1/gwas_res$SE[i]^2
      w2<-1/gwas_res0$SE[i]^2
      
      gwas_res_BETA_meta[i]<-(w1*gwas_res$BETA[i] + w2*gwas_res0$Beta[i])/(w1+w2)
      gwas_res_SE_meta[i]<-sqrt(1/(w1+w2))
      gwas_res_P_meta[i]<-pchisq(gwas_res_BETA_meta[i]^2/gwas_res_SE_meta[i]^2,1,lower.tail = F)
    }
    foldpath<-"/dcs04/nilanjan/data/rzhao/T2D_UKBB/"
    cur_iter<-0
    .<-capture.output(write_plink(paste0(foldpath,"ukbb_pop_",cur_iter),as.matrix(t(ref)),
      bim =  bim,fam = fam))
    ## COJO
    ma_file<- data.frame(SNP=bim$id,
      A1=bim$ref,A2=bim$alt,
      freq=gwas_res0$EAF,b=gwas_res_BETA_meta,
      se=gwas_res_SE_meta,p=gwas_res_P_meta,
      N=gwas_res0$Neff[1]+gwas_res$OBS_CT[1])
    
    readr::write_tsv(ma_file,paste0(foldpath,"ukbb_pop_",cur_iter,".ma"))
    
    aaa<-system(paste0("gcta64  --bfile ",foldpath,"ukbb_pop_",cur_iter," --cojo-file ",foldpath,"ukbb_pop_", cur_iter,".ma --cojo-slct --cojo-p 1e-5 --out ",foldpath,"test_",cur_iter),intern = TRUE)
    
    cojo_jma<-read.table(paste0(foldpath,"test_",cur_iter,".jma.cojo"),header = T)    
    
    aft_GWAS_jma<-cojo_jma$SNP[cojo_jma$pJ<0.05/N_SNP]
    #length(aft_GWAS_jma)
    #aft_GWAS_jma
    #gwas_res$ID[Nonnull_index]
    predict_nonnull_index2<-which(gwas_res$ID%in%aft_GWAS_jma)
    
    #gwas_res$P[predict_nonnull_index]
    #gwas_res0$Pvalue[predict_nonnull_index]
  }
  
  
  
  ### Run our method, 
  if(0){
    study_info<-lapply(1:nrow(gwas_res), function(i){
      study.m = list(Coeff=gwas_res$BETA[i],Covariance=gwas_res$SE[i]^2,Sample_size=gwas_res$OBS_CT[i],P=gwas_res$P[i])
      study.m
    })
  }
  study_info<-lapply(1:nrow(gwas_res0), function(i){
    study.m = list(Coeff=gwas_res0$Beta[i],Covariance=gwas_res0$SE[i]^2,Sample_size=gwas_res0$Neff[i],P=gwas_res0$Pvalue[i])
    study.m
  })
  study_info_scaled<-study_info
  EAF<-colsums(ref)/nrow(ref)/2
  for(i in 1:length(study_info)){
    cur_p<-gwas_res0$EAF[i]
    #cur_p<-EAF[i]
    a<-sqrt(2*cur_p*(1-cur_p))
    study_info_scaled[[i]]$Coeff<-study_info_scaled[[i]]$Coeff*a
    study_info_scaled[[i]]$Covariance<-study_info_scaled[[i]]$Covariance*a^2
  }
  
  source("~/multiSNP/adaptiveGMMLasso2.R")
  dim(ref)
  dim(pheno_covariates)
  index1<-which(pheno_covariates$Y_t2d == 1)
  index2<-sample(which(pheno_covariates$Y_t2d == 0),10*length(index1))
  
  ref_sp<-ref[c(index1,index2),]
  colnames(ref_sp)<-paste0('SNP',1:ncol(ref_sp))
  pheno_covariates_sp<-pheno_covariates[c(index1,index2),]
  
  #UKBB_pop_all<-data.frame(Phenotype=scale(pheno_covariates$Y_t2d,scale = F),scale(ref))
  UKBB_pop_all<-data.frame(Phenotype=scale(pheno_covariates_sp$Y_t2d,scale = F)/var(pheno_covariates_sp$Y_t2d),scale(ref_sp))
  colnames(UKBB_pop_all)[1]<-"Y"
  UKBB_pop_all<-as.matrix(UKBB_pop_all)
  #UKBB_full_fit2 <-cv.glmnet(x=cbind(UKBB_pop_all[,-1]),y=factor(pheno_covariates_sp$Y_t2d),family="binomial")
  #beta_initial2 = as.numeric(coef(UKBB_full_fit2, s = "lambda.min"))
  #bmi<-scale(pheno_covariates_sp$bmi)
  #aC<-adaptiveGMMlassoCBMI(as.matrix(UKBB_pop_all),bmi,study_info_scaled,cor_cut=0.8)
  #prop_pos<-aC$pos
  a<-adaptiveGMMlassoC(as.matrix(UKBB_pop_all),study_info_scaled,ld_cut = 0.9,cv_with_individual=F)
  prop_pos<-a$pos
  prop_pos
  
  newList<-list("cojo"=predict_nonnull_index,
    "cojometa"=predict_nonnull_index2,
    "proposed"=prop_pos)
  
  return(newList)
  #sum(bim$id%in%c(c("rs9930333","rs1558902","rs1121980","rs9936385","rs11642841")))
}
