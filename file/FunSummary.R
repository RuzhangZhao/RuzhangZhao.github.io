library(Seurat)
library(reticulate)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(ggrepel)
#library(sparsepca)
library(Rfast)
library(pdist)
library(Spectrum)
library(mclust)
library(stats)
library(utils)


#' savis
#'
#' savis: single-cell RNAseq adaptive visualiztaion
#'
#' @details This function argument to the function
#'
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA
#' @importFrom Rfast nth
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
#' @examples
#' a<-1
#'
savis<-function(
  expr_matrix,
  is_count_matrix=TRUE,
  npcs = 10,
  nfeatures = 2000,
  cluster_method = "louvain",
  resolution = 0.5,
  resolution_sub = NULL,
  process_min_size = NULL,
  process_min_count = NULL,
  verbose = TRUE,
  adaptive = FALSE,
  show_cluster = FALSE,
  return_cluster = FALSE,
  verbose_more = FALSE
){
  if(is.null(resolution_sub)){
    resolution_sub<-resolution
  }
  if (nrow(expr_matrix) < nfeatures){
    stop("nfeatures should be smaller than
      the number of features in expression
      matrix")
  }
  if (ncol(expr_matrix) < npcs){
    stop("npcs(number of PC) should be smaller
      than the number of samples in expression
      matrix")
  }
  if(is_count_matrix){
    if(verbose){
      cat('\n')
      print("Normalizing Expression Matrix...")
      pb <- txtProgressBar(min = 0, max = 20, style = 3, file = stderr())
    }
    
    expr_matrix<-NormalizeData(
      expr_matrix,
      verbose = verbose_more)
  }
  expr_matrix_process<-expr_matrix
  if(verbose){
    cat('\n')
    print("Finding Variable Features...")
    setTxtProgressBar(pb = pb, value = 1)
  }
  expr_matrix_hvg <- FindVariableFeatures(
    expr_matrix_process,
    verbose = verbose_more)$vst.variance.standardized
  hvg<- nth(x = expr_matrix_hvg,
    k = nfeatures,
    num.of.nths = 2,
    descending = T,
    index.return = T)[,1]
  #hvg<- which(expr_hvg %in% sort(expr_hvg,decreasing = T)[1:2000])
  expr_matrix_process<-expr_matrix_process[hvg,]
  if(verbose){
    cat('\n')
    print("Scaling Expression Matrix...")
    setTxtProgressBar(pb = pb, value = 2)
  }
  expr_matrix_process <- ScaleData(
    expr_matrix_process,
    verbose = verbose_more)
  
  if(verbose){
    cat('\n')
    print("Calculating Global PCA...")
    setTxtProgressBar(pb = pb, value = 3)
  }
  suppressWarnings(expr_matrix_pca <- RunPCA(
    object = expr_matrix_process,
    features = rownames(expr_matrix_process),
    npcs = npcs,
    verbose = verbose_more)@cell.embeddings)
  rm(expr_matrix_process)
  expr_matrix_pca<-data.frame(expr_matrix_pca)
  if(verbose){
    cat('\n')
    print("Doing Clustering...")
    setTxtProgressBar(pb = pb, value = 5)
  }
  res_cluster<-DoCluster(
    pc_embedding = expr_matrix_pca,
    method = cluster_method,
    resolution = resolution)
  size_cluster<-c()
  for ( i in unique(res_cluster$cluster)){
    size_cluster<-c(size_cluster,
      sum(res_cluster$cluster == i))
  }
  if(verbose){
    cat('\n')
    print(paste0("Clustering Results:",
      length(unique(res_cluster$cluster)),
      " clusters."))
    setTxtProgressBar(pb = pb, value = 6)
  }
  if(verbose){
    if(show_cluster){
      cat('\n')
      print("Size of Cluster:")
      print(size_cluster)
    }
  }
  cluster_label<-res_cluster$cluster
  
  if(verbose){
    cat('\n')
    print("Calculating Local PCA...")
    setTxtProgressBar(pb = pb, value = 8)
  }
  sub_PC_list<-SubPCEmbedding(
    expr_matrix = expr_matrix,
    cluster_label = cluster_label,
    npcs = npcs,
    nfeatures = nfeatures)
  #rm(expr_matrix)
  combined_embedding<-CombinePC(
    PC_embedding = expr_matrix_pca,
    cluster_label = cluster_label,
    sub_PC_list = sub_PC_list)
  if(verbose){
    cat('\n')
    print("Scaling Global PCs and Local PCs...")
    setTxtProgressBar(pb = pb, value = 9)
  }
  combined_res<-ScaleFactor(
    combined_embedding = combined_embedding,
    cluster_label = cluster_label,
    npcs = npcs,
    center_method = "mean")
  combined_scale_factor<-combined_res$scale_factor
  if(adaptive){
    if (!is.null(process_min_count)){
      process_min_size<-sort(size_cluster,decreasing = T)[process_min_count]
    }else{
      if (is.null(process_min_size)){
        process_min_size<-mean(size_cluster)
      }
    }
    if(verbose){
      cat('\n')
      print(paste0("Process_min_size: ",process_min_size))
      setTxtProgressBar(pb = pb, value = 9.5)
    }
    combined_embedding_sub<-combined_embedding[,((npcs+1):(2*npcs))]
    # sorted unique cluster label
    label_index<-sort(as.numeric(
      unique(as.character(cluster_label))))
    # number of cluster label
    N_label<-length(label_index)
    
    if(verbose){
      cat('\n')
      print("Exploring if clusters can be separated further...")
      setTxtProgressBar(pb = pb, value = 10)
    }
    
    cluster_label_i<-lapply(1:N_label, function(i){
      index_i<-which(cluster_label == label_index[i])
      res_cluster_i<-DoCluster(
        pc_embedding = combined_embedding_sub[index_i,],
        method = cluster_method,
        resolution = resolution_sub)
      res_cluster_i$cluster
    })
    
    continue_adaptive<-FALSE
    size_cluster_sub<-c()
    for ( i in 1:length(cluster_label_i)){
      if(length(unique(cluster_label_i[[i]]))>1){
        continue_adaptive<-TRUE
      }
      size_cluster_sub<-c(size_cluster_sub,
        length(unique(cluster_label_i[[i]])))
    }
    
    if(continue_adaptive){
      if(verbose){
        cat('\n')
        print("Can be separated further...")
        setTxtProgressBar(pb = pb, value = 11)
        if(show_cluster){
          cat('\n')
          print("Separation of Each Cluster:")
          print(size_cluster_sub)
        }
      }
      combined_embedding<-AdaptiveCombine(expr_matrix = expr_matrix,
        combined_embedding = combined_embedding,
        npcs = npcs,
        nfeatures = nfeatures,
        process_min_size = process_min_size,
        cluster_label = cluster_label,
        cluster_label_i = cluster_label_i,
        resolution = resolution_sub,
        combined_scale_factor = combined_scale_factor)
      cluster_label_sub<-rep(-1,length(cluster_label))
      for ( i in 1:N_label){
        if (length(unique(cluster_label_i[[i]]))>1){
          index_i<-which(cluster_label == label_index[i])
          cluster_label_sub[index_i]<-cluster_label_i[[i]]
        }else{
          index_i<-which(cluster_label == label_index[i])
          cluster_label_sub[index_i]<- -1
        }
        
      }
      combined_embedding<-data.frame(
        "cluster_label"=cluster_label,
        "cluster_label_sub"=cluster_label_sub,
        combined_embedding)
      if(verbose){
        cat('\n')
        print("Last Step:Running Adaptive UMAP...")
        setTxtProgressBar(pb = pb, value = 12)
      }
      umap_embedding<-RunAdaUMAP(X = combined_embedding,
        metric = 'euclidean2')
      if(verbose){
        cat('\n')
        print("Finished...")
        setTxtProgressBar(pb = pb, value = 20)
      }
      if(return_cluster){
        newList<-list("umap_embedding"=umap_embedding,
          "cluster_label"=cluster_label,
          "cluster_label_sub"=cluster_label_sub)
        return(newList)
      }else{
        return(umap_embedding )
      }
    }else{
      if(verbose){
        cat('\n')
        print("Cannot be separated further...")
        setTxtProgressBar(pb = pb, value = 11)
      }
      combined_embedding<-combined_res$combined_embedding
      combined_embedding<-data.frame(
        "cluster_label"=cluster_label,
        combined_embedding)
      if(verbose){
        cat('\n')
        print("Last Step:Running Adaptive UMAP...")
        setTxtProgressBar(pb = pb, value = 12)
      }
      umap_embedding<-RunAdaUMAP(X = combined_embedding)
      if(verbose){
        cat('\n')
        print("Finished...")
        setTxtProgressBar(pb = pb, value = 20)
      }
      if(return_cluster){
        newList<-list("umap_embedding"=umap_embedding,
          "cluster_label"=cluster_label)
        return(newList)
      }else{
        return(umap_embedding)
      }
    }
    
    
  }else{
    
    combined_embedding<-combined_res$combined_embedding
    combined_embedding<-data.frame(
      "cluster_label"=cluster_label,
      combined_embedding)
    if(verbose){
      cat('\n')
      print("Last Step:Running Adaptive UMAP...")
      setTxtProgressBar(pb = pb, value = 10)
    }
    umap_embedding<-RunAdaUMAP(X = combined_embedding)
    if(verbose){
      cat('\n')
      print("Finished...")
      setTxtProgressBar(pb = pb, value = 20)
    }
    if(return_cluster){
      newList<-list("umap_embedding"=umap_embedding,
        "cluster_label"=cluster_label)
      return(newList)
    }else{
      return(umap_embedding)
    }
  }
  
  
  
}






#' RunAdaUMAP
#'
#' Use Adaptive Distance Metric to Run UMAP
#'
#' @details This function argument to the function
#'
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom reticulate source_python py_module_available py_set_seed import
#'
#' @export
#'
#' @examples
#' a<-1
#'
RunAdaUMAP<-function(
  X,
  n.neighbors = 30L,
  n.components = 2L,
  metric = 'euclidean',
  n.epochs = NULL,
  learning.rate = 1.0,
  min.dist = 0.3,
  spread = 1.0,
  set.op.mix.ratio = 1.0,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5,
  a = NULL,
  b = NULL,
  uwot.sgd = FALSE,
  seed.use = 42L,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  verbose = FALSE){
  if (!py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn) or through reticulate package (e.g. reticulate::py_module_available(module = 'umap')")
  }
  adaptive_distance_url<-"http://ruzhangzhao.github.io/file/adaptive_distance.py"
  source_python(adaptive_distance_url)
  if (!is.null(x = seed.use)) {
    py_set_seed(seed = seed.use)
  }
  if (typeof(x = n.epochs) == "double") {
    n.epochs <- as.integer(x = n.epochs)
  }
  
  if (metric == "euclidean"){
    adaptive_metric<-adaptive_euclidean_grad
  }
  if (metric == "euclidean2"){
    adaptive_metric<-adaptive_euclidean2_grad
  }
  umap_import <- import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(x = n.neighbors),
    n_components = as.integer(x = n.components),
    metric = adaptive_metric,
    n_epochs = n.epochs,
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a,
    b = b,
    metric_kwds = metric.kwds,
    angular_rp_forest = angular.rp.forest,
    random_state=seed.use,
    verbose = verbose
  )
  umap_embedding<-umap$fit_transform(as.matrix(x = X))
  umap_embedding<-data.frame(umap_embedding)
  colnames(umap_embedding)<-paste0("UMAP_",1:n.components)
  umap_embedding
}




#' ScaleFactor
#'
#' Combined PC embedding with scale factor for subPC
#'
#' @details Adapts source code from: https://github.com/rstudio/reticulate/blob/master/R/source.R
#' @details Now source_python_ supports loading from http
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom Rfast Dist
#' @importFrom pdist pdist
#'
#' @export
#'
#' @examples
#' a<-1
#'
#'
ScaleFactor<-function(
  combined_embedding,
  cluster_label,
  npcs=NULL,
  center_method = "mean",
  return_factor=TRUE){
  
  # sorted unique cluster label
  label_index<-sort(as.numeric(
    unique(as.character(cluster_label))))
  # number of cluster label
  N_label<-length(label_index)
  if (is.null(npcs)){
    npcs<-ncol(combined_embedding)/2
    if (npcs %% 1 != 0){
      stop("combined_embedding combines global PC
        and local PC, column number of combined
        embedding should be even.")
    }
  }else{
    if (2*npcs != ncol(combined_embedding)){
      stop("column number of combined embedding
        is not equal to 2*npcs")
    }
  }
  
  if (center_method == "geometry"){
    
    cluster_center<-t(sapply(1:N_label, function(i){
      index_i<-which(cluster_label == label_index[i])
      tmp<-sapply(1:(2*npcs), function(j){
        (max(combined_embedding[index_i,j]) +
            min(combined_embedding[index_i,j]))/2
      })
      names(tmp)<-c(paste0("PC_",1:npcs),
        paste0("subPC_",1:npcs))
      tmp
    }))
    
  }else if (center_method == "mean"){
    # cluster center by mean
    cluster_center<-t(sapply(1:N_label, function(i){
      index_i<-which(cluster_label == label_index[i])
      colMeans(combined_embedding[index_i,1:(2*npcs)])
    }))
    colnames(cluster_center)<-c(paste0("PC_",1:npcs),
      paste0("subPC_",1:npcs))
    
  }else{
    stop("center_method is chosen from geometry or mean")
  }
  
  # distance matrix for cluster center
  cluster_center_dist<-Dist(cluster_center[,1:npcs])
  diag(cluster_center_dist)<-NA
  
  # pdist is used here, which is fast
  cluster_semi_d<-sapply(1:N_label, function(i){
    index_i<-which(cluster_label == label_index[i])
    max((pdist(cluster_center[i,(npcs+1):(2*npcs)],
      combined_embedding[index_i,(npcs+1):(2*npcs)]))@dist)
  })
  
  scale_factor<-sapply(1:N_label, function(i){
    min(cluster_center_dist[i,],na.rm = T)/3/cluster_semi_d[i]
  })
  
  for ( i in 1:N_label){
    index_i<-which(cluster_label == label_index[i])
    combined_embedding[index_i,(1+npcs):(2*npcs)]<-
      combined_embedding[index_i,(1+npcs):(2*npcs)]*scale_factor[i]
  }
  if(return_factor){
    newList<-list("combined_embedding"=combined_embedding,
      "scale_factor"=scale_factor)
    return(newList)
  }
  combined_embedding
}

#' DoCluster
#'
#' Do cluster for PC embedding
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom Spectrum Spectrum
#'
#' @export
#'
#' @examples
#' a<-1
#'
#'
DoCluster<-function(
  pc_embedding,
  method = "louvain",
  resolution = 0.5,
  subcluster = FALSE,
  fixk = NULL,
  seed.use = 42,
  verbose = FALSE,
  verbose_more = FALSE){
  if (verbose_more){
    verbose <- TRUE
  }
  set.seed(seed.use)
  if (method == "louvain"){
    if (verbose){
      print("Finding Neighbors...")
    }
    snn_<- FindNeighbors(object = pc_embedding,
      verbose = verbose_more)$snn
    if (verbose){
      print("Finding Clusters...")
    }
    cluster_ <- FindClusters(snn_,
      resolution = resolution,
      verbose = verbose_more)[[1]]
    cluster_ <- as.numeric(as.character(cluster_))
    if (subcluster){
      
      N_cluster<-length(unique(cluster_))
      subcluster_<-rep(NA,length(cluster_))
      for ( i in unique(cluster_)){
        if (verbose){
          print(paste0("Processing cluster ",i))
        }
        index_i<-which(cluster_ == i)
        if (verbose_more){
          print("Finding Neighbors...")
        }
        snn_<- FindNeighbors(object = pc_embedding[index_i,],
          verbose = verbose_more)$snn
        if (verbose_more){
          print("Finding Clusters...")
        }
        subcluster_label <- FindClusters(snn_,
          resolution = resolution,
          verbose = verbose_more)[[1]]
        subcluster_[index_i]<-paste0(i,subcluster_label)
      }
      subcluster_ <- as.numeric(as.character(subcluster_))
    }
  }
  if (method == "spectral"){
    num_example<-nrow(pc_embedding)
    N_center<-min(1000,num_example%/%5)
    set.seed(42)
    suppressMessages(
      spectral_eg<-Spectrum(t(
        as.matrix(pc_embedding)),
        method = 1,showres = verbose_more,
        FASP = T,FASPk = N_center))
    cluster_<-spectral_eg$allsample_assignments
    cluster_<-as.numeric(as.character(cluster_))
    if (subcluster){
      N_cluster<-length(unique(cluster_))
      subcluster_<-rep(NA,length(cluster_))
      for ( i in unique(cluster_)){
        print(paste0("Processing cluster ",i))
        index_i<-which(cluster_ == i)
        if (length(index_i) > min(100,N_center)){
          suppressMessages(
            spectral_eg<-Spectrum(t(
              as.matrix(pc_embedding[index_i,])),
              method = 1,FASP = T,
              FASPk = N_center,
              showres = F))
        }else{
          suppressMessages(
            spectral_eg<-Spectrum(t(
              as.matrix(pc_embedding[index_i,])),
              method = 1,FASP = T,
              FASPk = length(index_i)%/%3,
              showres = F))
        }
        subcluster_[index_i]<-paste0(i,
          spectral_eg$allsample_assignments)
      }
      subcluster_<-as.numeric(as.character(subcluster_))
    }
  }
  if (subcluster){
    newList<-list("cluster" = cluster_,
      "subcluster" = subcluster_)
    return(newList)
  }else{
    newList<-list("cluster" = cluster_)
    return(newList)
  }
}

#' SubPCEmbedding
#'
#' Based on clustering, calculate the subPC embedding for each cluster
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom Seurat FindVariableFeatures ScaleData RunPCA
#' @importFrom Rfast nth
#'
#' @export
#'
#' @examples
#' a<-1
#'
#'
SubPCEmbedding<-function(
  expr_matrix,
  cluster_label,
  npcs=10,
  nfeatures =2000,
  return_hvg=FALSE,
  verbose = FALSE
  #  process_min_size = 0
){
  # get unique cluster label
  label_index<-sort(as.numeric(
    unique(as.character(cluster_label))))
  # number of cluster
  N_label<-length(label_index)
  
  # get PC for each cluster
  expr_cluster_pc<-list()
  hvg_list<-list()
  hvg_name_list<-list()
  for ( i in c(1:N_label)){
    index_i<-which(cluster_label == label_index[i])
    #if(length(index_i) > process_min_size){
    if (verbose){
      print(paste0("Process cluster ",
        label_index[i]," Yes",
        "(size=",
        length(index_i),")"))
    }
    # expr_matrix is initial expression matrix (gene*cell)
    expr_tmp<-expr_matrix[,index_i]
    #expr_tmp<- NormalizeData(expr_tmp,
    #  verbose = F)
    expr_tmp_hvg <- FindVariableFeatures(
      expr_tmp,verbose = F)$vst.variance.standardized
    #tmp_hvg<- which(expr_tmp_hvg %in%
    #    sort(expr_tmp_hvg,decreasing = T)[1:nfeatures])
    tmp_hvg<- nth(x = expr_tmp_hvg,
      k = nfeatures,
      num.of.nths = 2,
      descending = T,
      index.return = T)[,1]
    hvg_list[[i]]<-tmp_hvg
    hvg_name_list[[i]]<-rownames(expr_matrix)[tmp_hvg]
    expr_tmp<-expr_tmp[tmp_hvg,]
    expr_tmp <-ScaleData(expr_tmp,verbose = F)
    rm(expr_tmp_hvg,tmp_hvg,index_i)
    
    suppressWarnings(
      expr_tmp <- RunPCA(
        object = expr_tmp,
        features = rownames(expr_matrix),
        npcs = npcs,
        verbose = F)@cell.embeddings)
    
    # Get subgroup PC with scale factor(default to be 1)
    expr_cluster_pc[[i]]<-expr_tmp[,c(1:npcs)]
    
    
    #}
    #else{
    #  if (verbose){
    #    print(paste0("Process cluster ",
    #      label_index[i]," No",
    #      "(size=",
    #      length(index_i),")"))
    #  }
    
    # if the cluster size less than process_min_size,
    # just fill in NA
    #  expr_cluster_pc[[i]]<-NA
    #expr_cluster_pc[[i]]<-matrix(0,
    # nrow = length(index_i),
    # ncol = npcs)
    #}
    #####
  }
  if (return_hvg){
    newList<-list("sub_PC_list" = expr_cluster_pc,
      "hvg_index"=hvg_list,
      "hvg_name"=hvg_name_list)
    return(newList)
  }else{
    return(expr_cluster_pc)
  }
}

#' CombinePC
#'
#' Combine PC and subPCs
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#'
#' @export
#'
#' @examples
#' a<-1
#'
#'
CombinePC<-function(PC_embedding,
  cluster_label,
  sub_PC_list,
  scale_factor=NULL
){
  label_index<-sort(as.numeric(unique(
    as.character(cluster_label))))
  N_label<-length(label_index)
  if (N_label != length(sub_PC_list)){
    stop("cluster number is not equal to sub PC list length")
  }
  npcs<-ncol(sub_PC_list[[1]])
  PC_embedding_combined<-PC_embedding
  
  sub_PC_list_supp<-matrix(0,nrow = nrow(PC_embedding),
    ncol = npcs)
  sub_PC_list_supp<-data.frame(sub_PC_list_supp)
  rownames(sub_PC_list_supp)<-rownames(PC_embedding)
  
  
  for (i in c(1:N_label)){
    index_i<-which(cluster_label == label_index[i])
    if (is.null(scale_factor[1])){
      if (!is.na(sub_PC_list[[i]][1])){
        sub_PC_list_supp[index_i,]<-sub_PC_list[[i]]
      }else{
        sub_PC_list_supp[index_i,]<-0
        #PC_embedding[index_i,]
      }
    }else{
      if (!is.na(sub_PC_list[[i]][1])){
        sub_PC_list_supp[index_i,]<-sub_PC_list[[i]]*scale_factor[i]
      }else{
        sub_PC_list_supp[index_i,]<-0
        #PC_embedding[index_i,]
      }
    }
    
  }
  
  PC_embedding_combined<-cbind(PC_embedding_combined,
    sub_PC_list_supp)
  
  return(PC_embedding_combined)
}


#' AdaptiveCombine
#'
#' Adaptively Combine PC and subPCs
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#'
#' @export
#'
#' @examples
#' a<-1
#'
#'
#'
AdaptiveCombine<-function(expr_matrix,
  combined_embedding,
  cluster_label,
  cluster_label_i,
  combined_scale_factor,
  npcs=10,
  nfeatures=2000,
  process_min_size=0,
  resolution=0.5
){
  expr_matrix_pca<-combined_embedding[,(1:(npcs))]
  combined_embedding_sub<-combined_embedding[,((npcs+1):(2*npcs))]
  # sorted unique cluster label
  label_index<-sort(as.numeric(
    unique(as.character(cluster_label))))
  # number of cluster label
  N_label<-length(label_index)
  
  sub_PC_list_i<-lapply(1:N_label, function(i){
    
    index_i<-which(cluster_label == label_index[i])
    if(length(unique(cluster_label_i[[i]]))>1 &
        length(index_i) >= process_min_size){
      sub_PC_list_i<-SubPCEmbedding(
        expr_matrix = expr_matrix[,index_i],
        cluster_label = cluster_label_i[[i]],
        npcs = npcs,
        nfeatures = nfeatures)
      
      combined_embedding_i<-CombinePC(
        PC_embedding = combined_embedding_sub[index_i,],
        cluster_label = cluster_label_i[[i]],
        sub_PC_list = sub_PC_list_i)
      
      combined_embedding_i<-ScaleFactor(
        combined_embedding = combined_embedding_i,
        cluster_label = cluster_label_i[[i]],
        npcs = npcs,
        center_method = "mean")$combined_embedding
      combined_embedding_i<-as.matrix(combined_embedding_i)
      colnames(combined_embedding_i)<-c(paste0("subPC",1:npcs),
        paste0("subsubPC",1:npcs))
    }else{
      combined_embedding_i<-matrix(0,
        nrow = length(index_i),
        ncol = npcs)
      combined_embedding_i<-cbind(combined_embedding_sub[index_i,],
        combined_embedding_i)
      combined_embedding_i<-as.matrix(combined_embedding_i)
      colnames(combined_embedding_i)<-c(paste0("subPC",1:npcs),
        paste0("subsubPC",1:npcs))
    }
    combined_embedding_i
  })
  
  combined_embedding<-CombinePC(PC_embedding = expr_matrix_pca,
    cluster_label = cluster_label,
    sub_PC_list = sub_PC_list_i,
    scale_factor = combined_scale_factor)
  
  
  return(combined_embedding)
}



#' adaDimPlot
#'
#' Adaptively Plot the UMAP Embedding
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom dplyr group_by summarise
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats median
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' a<-1
#'
#'
#'
adaDimPlot<-function(
  umap_embedding,
  label,
  pt.size=0,
  show.legend=TRUE
){
  umap_embedding<-data.frame(umap_embedding)
  colnames(umap_embedding)<-paste0("UMAP_",1:ncol(umap_embedding))
  umap_embedding$label<-factor(label)
  
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal,
    qual_col_pals$maxcolors,rownames(qual_col_pals)))
  
  
  gg<-ggplot(umap_embedding)+
    geom_point(aes(x = UMAP_1,
      y = UMAP_2,
      color = label),
      size = pt.size)+
    theme(legend.title = element_blank())+
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.key=element_blank())+
    guides(color = guide_legend(override.aes =
        list(size=3)))+
    scale_colour_manual(values =
        col_vector[c(1:length(unique(label)))])
  if (!show.legend){
    gg<-gg+theme(legend.position="none")
  }
  
  centers<-group_by(umap_embedding,label)
  centers<-summarise(centers,x = median(x = UMAP_1),
    y = median(x = UMAP_2),.groups = 'drop')
  
  gg <- gg +
    geom_text_repel(data = centers,
      mapping = aes(x = x, y = y,
        label = label), size = 4)
  gg
}

#' ARIEvaluate
#'
#' Adaptively Plot the UMAP Embedding
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom Spectrum Spectrum
#' @importFrom stats kmeans
#' @importFrom mclust adjustedRandIndex
#'
#' @export
#'
#' @examples
#' a<-1
#'
#'
#'
#'
#'
ARIEvaluate<-function(
  test_data,
  label,
  method = "louvain",
  resolution = 0.5,
  fixk = NULL
){
  if (method == "louvain"){
    snn_<- FindNeighbors(object = test_data,
      verbose = F)$snn
    cluster_label <- FindClusters(snn_,
      resolution = resolution,
      verbose = F)[[1]]
  }
  if (method == "spectral"){
    
    set.seed(42)
    N_center<-nrow(test_data)
    if (N_center > 5000){
      N_center <- 1000
    }else{
      N_center<-ceiling(nrow(test_data)/5)
    }
    if (is.null(fixk)){
      suppressMessages(
        spectrual_eg<-Spectrum(t(as.matrix(
          test_data)),
          method = 1,showres = F,
          FASP = T,FASPk = N_center))
    }else{
      suppressMessages(
        spectrual_eg<-Spectrum(t(as.matrix(
          test_data)),
          method = 3,showres = F,
          FASP = T,FASPk = N_center,fixk = fixk))
    }
    
    cluster_label<-spectrual_eg$allsample_assignments
  }
  if (method == "kmeans"){
    res<-kmeans(x = as.matrix(
      test_data),centers = fixk,iter.max = 100)
    cluster_label<-res$cluster
  }
  
  
  c(length(unique(cluster_label)),
    adjustedRandIndex(cluster_label,label))
}



