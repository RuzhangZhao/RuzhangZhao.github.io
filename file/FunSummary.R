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
library(uwot)
library(glue)
library(MASS)
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
  distance_metric = "euclidean",
  cluster_method = "louvain",
  resolution = 0.5,
  resolution_sub = 0,
  adaptive = FALSE,
  max_stratification = 3,
  scale_factor_separation =3,
  process_min_size = NULL,
  process_min_count = NULL,
  center_method = "mean",
  verbose = TRUE,
  show_cluster = FALSE,
  return_cluster = FALSE,
  verbose_more = FALSE,
  run_adaUMAP = TRUE,
  adjust_UMAP = TRUE,
  adjust_method = "umap",
  adjust_rotate = TRUE,
  check_differential = TRUE,
  seed.use = 42L
){
  if(max_stratification == 1){
    stop("Please directly use umap: savis 
      supports adaptive visualization for 
      max_stratification larger than 1.")
  }
  
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
  if(is.null(rownames(expr_matrix)[1])){
    rownames(expr_matrix)<-c(1:nrow(expr_matrix))
  }
  if(is.null(colnames(expr_matrix)[1])){
    colnames(expr_matrix)<-c(1:ncol(expr_matrix))
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
  }else{
    pb <- txtProgressBar(min = 0, max = 20, style = 3, file = stderr())
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
  cluster_label<-DoCluster(
    pc_embedding = expr_matrix_pca,
    method = cluster_method,
    resolution = resolution)$cluster
  # sorted unique cluster label
  label_index<-sort(as.numeric(
    unique(as.character(cluster_label))))
  # number of cluster label
  N_label<-length(label_index)
  
  size_cluster<-c()
  for ( i in 1:N_label){
    size_cluster<-c(size_cluster,
      sum(cluster_label == label_index[i]))
  }
  if(verbose){
    cat('\n')
    print(paste0("Clustering Results:",
      length(unique(cluster_label)),
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
  
  if(verbose){
    cat('\n')
    print("Calculating Local PCA...")
    setTxtProgressBar(pb = pb, value = 8)
  }
  if(max_stratification == 2 | adaptive == FALSE){
    
    combined_embedding<-FormCombinedEmbedding(
      expr_matrix=expr_matrix,
      expr_matrix_pca=expr_matrix_pca,
      cluster_label=cluster_label,
      npcs=npcs,
      nfeatures =nfeatures,
      center_method = center_method,
      scale_factor_separation=scale_factor_separation)
    adaptive<-FALSE
  }
  
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
    
    
    if(verbose){
      cat('\n')
      print("Exploring if clusters can be separated further...")
      setTxtProgressBar(pb = pb, value = 10)
    }
    
    umap_res<-FormAdaptiveCombineList(
      expr_matrix=expr_matrix,
      expr_matrix_pca=expr_matrix_pca,
      max_stratification=max_stratification,
      stratification_count=1,
      scale_factor_separation = scale_factor_separation,
      resolution=resolution_sub,
      cluster_method=cluster_method,
      npcs=npcs,
      nfeatures=nfeatures,
      process_min_size=process_min_size,
      do_cluster = FALSE,
      cluster_label = cluster_label,
      check_differential = check_differential,
      verbose = verbose_more)
    cluster_label<-umap_res$cluster_label
    if(is.null(dim(cluster_label)[1])){
      combined_embedding<-data.frame(
        "Layer2"=cluster_label,
        umap_res$combined_embedding)
      if(ncol(umap_res$combined_embedding)!=(2*npcs)){
        stop("label and combined embedding size do not match")
      }
      if(run_adaUMAP){
        if(verbose){
          cat('\n')
          print("Running Adaptive UMAP...")
          setTxtProgressBar(pb = pb, value = 12)
        }
        umap_embedding<-RunAdaUMAP(
          X = combined_embedding,
          metric = distance_metric,
          metric_count = 1,
          py_envir = parent.frame(),
          seed.use = seed.use)
      }
      
      
      if(run_adaUMAP){
        if(adjust_UMAP){
          if(verbose){
            cat('\n')
            print("Adjusting UMAP...")
            setTxtProgressBar(pb = pb, value = 18)
          }
          expr_matrix_umap<-umap(
            X = expr_matrix_pca,
            a = 1.8956, 
            b = 0.8006, 
            metric = distance_metric
          )
          umap_embedding<-adjustUMAP(
            pca_embedding = expr_matrix_pca,
            umap_embedding = umap_embedding,
            global_umap_embedding = expr_matrix_umap,
            adjust_method = adjust_method,
            rotate = adjust_rotate,
            seed.use = seed.use)
        }
        if(return_cluster){
          newList<-list("umap_embedding"=umap_embedding,
            "cluster_label"=cluster_label)
          
        }else{
          newList<-umap_embedding
        }
      }else{
        if(return_cluster){
          newList<-list("combined_embedding"=combined_embedding,
            "cluster_label"=cluster_label)
          
        }else{
          newList<-combined_embedding
        }
      }
    }else{
      if (ncol(cluster_label) == 2) {
        colnames(cluster_label)<-paste0("Layer",
          2:(ncol(cluster_label)+1))
        combined_embedding<-data.frame(
          cluster_label,
          umap_res$combined_embedding)
        if(ncol(umap_res$combined_embedding)!=(3*npcs)){
          stop("label and combined embedding size do not match")
        }
        if(run_adaUMAP){
          if(verbose){
            cat('\n')
            print("Running Adaptive UMAP...")
            setTxtProgressBar(pb = pb, value = 12)
          }
          umap_embedding<-RunAdaUMAP(
            X = combined_embedding,
            metric = distance_metric,
            metric_count = 2,
            py_envir = parent.frame(),
            seed.use = seed.use)
        }
        
        
        
        if(run_adaUMAP){
          if(adjust_UMAP){
            if(verbose){
              cat('\n')
              print("Adjusting UMAP...")
              setTxtProgressBar(pb = pb, value = 18)
            }
            expr_matrix_umap<-umap(
              X = expr_matrix_pca,
              a = 1.8956, 
              b = 0.8006, 
              metric = distance_metric
            )
            umap_embedding<-adjustUMAP(
              pca_embedding = expr_matrix_pca,
              umap_embedding = umap_embedding,
              global_umap_embedding = expr_matrix_umap,
              adjust_method = adjust_method,
              rotate = adjust_rotate,
              seed.use = seed.use)
          }
          if(return_cluster){
            newList<-list("umap_embedding"=umap_embedding,
              "cluster_label"=cluster_label)
            
          }else{
            newList<-umap_embedding
          }
        }else{
          if(return_cluster){
            newList<-list("combined_embedding"=combined_embedding,
              "cluster_label"=cluster_label)
            
          }else{
            newList<-combined_embedding
          }
        }
        
      }else if (ncol(cluster_label) == 3){
        colnames(cluster_label)<-paste0("Layer",
          2:(ncol(cluster_label)+1))
        combined_embedding<-data.frame(
          cluster_label,
          umap_res$combined_embedding)
        if(ncol(umap_res$combined_embedding)!=(4*npcs)){
          stop("label and combined embedding size do not match")
        }
        if(run_adaUMAP){
          if(verbose){
            cat('\n')
            print("Running Adaptive UMAP...")
            setTxtProgressBar(pb = pb, value = 12)
          }
          umap_embedding<-RunAdaUMAP(
            X = combined_embedding,
            metric = distance_metric,
            metric_count = 3,
            py_envir = parent.frame(),
            seed.use = seed.use)
          
        }
        
        
        if(run_adaUMAP){
          if(adjust_UMAP){
            if(verbose){
              cat('\n')
              print("Adjusting UMAP...")
              setTxtProgressBar(pb = pb, value = 18)
            }
            expr_matrix_umap<-umap(
              X = expr_matrix_pca,
              a = 1.8956, 
              b = 0.8006, 
              metric = distance_metric
            )
            umap_embedding<-adjustUMAP(
              pca_embedding = expr_matrix_pca,
              umap_embedding = umap_embedding,
              global_umap_embedding = expr_matrix_umap,
              adjust_method = adjust_method,
              rotate = adjust_rotate,
              seed.use = seed.use)
          }
          if(return_cluster){
            newList<-list("umap_embedding"=umap_embedding,
              "cluster_label"=cluster_label)
            
          }else{
            newList<-umap_embedding
          }
        }else{
          if(return_cluster){
            newList<-list("combined_embedding"=combined_embedding,
              "cluster_label"=cluster_label)
            
          }else{
            newList<-combined_embedding
          }
        }
        
      }else{
        colnames(cluster_label)<-paste0("Layer",
          2:(ncol(cluster_label)+1))
        combined_embedding<-data.frame("Num_Layer"=(ncol(cluster_label)+1),
          cluster_label,
          umap_res$combined_embedding)
        if(run_adaUMAP){
          if(verbose){
            cat('\n')
            print("Running Adaptive UMAP...")
            setTxtProgressBar(pb = pb, value = 12)
          }
          umap_embedding<-RunAdaUMAP(
            X = combined_embedding,
            metric = distance_metric,
            metric_count = 4,
            py_envir = parent.frame(),
            seed.use = seed.use)  
        }
        
        
        if(run_adaUMAP){
          if(adjust_UMAP){
            if(verbose){
              cat('\n')
              print("Adjusting UMAP...")
              setTxtProgressBar(pb = pb, value = 18)
            }
            expr_matrix_umap<-umap(
              X = expr_matrix_pca,
              a = 1.8956, 
              b = 0.8006, 
              metric = distance_metric
            )
            umap_embedding<-adjustUMAP(
              pca_embedding = expr_matrix_pca,
              umap_embedding = umap_embedding,
              global_umap_embedding = expr_matrix_umap,
              adjust_method = adjust_method,
              rotate = adjust_rotate,
              seed.use = seed.use)
          }
          if(return_cluster){
            newList<-list("umap_embedding"=umap_embedding,
              "cluster_label"=cluster_label)
            
          }else{
            newList<-umap_embedding
          }
        }else{
          if(return_cluster){
            newList<-list("combined_embedding"=combined_embedding,
              "cluster_label"=cluster_label)
          }else{
            newList<-combined_embedding
          }
        }
        
      }
    }
    
    
    
  }else{
    
    combined_embedding<-data.frame(
      "cluster_label"=cluster_label,
      combined_embedding)
    
    if(run_adaUMAP){
      if(verbose){
        cat('\n')
        print("Running Adaptive UMAP...")
        setTxtProgressBar(pb = pb, value = 10)
      }
      umap_embedding<-RunAdaUMAP(
        X = combined_embedding,
        metric = distance_metric,
        metric_count = 1,
        py_envir = parent.frame(),
        seed.use = seed.use)
    }
    
    
    if(run_adaUMAP){
      if(adjust_UMAP){
        if(verbose){
          cat('\n')
          print("Adjusting UMAP...")
          setTxtProgressBar(pb = pb, value = 18)
        }
        expr_matrix_umap<-umap(
          X = expr_matrix_pca,
          a = 1.8956, 
          b = 0.8006, 
          metric = distance_metric
        )
        umap_embedding<-adjustUMAP(
          pca_embedding = expr_matrix_pca,
          umap_embedding = umap_embedding,
          global_umap_embedding = expr_matrix_umap,
          adjust_method = adjust_method,
          rotate = adjust_rotate,
          seed.use = seed.use)
      }
      if(return_cluster){
        newList<-list("umap_embedding"=umap_embedding,
          "cluster_label"=cluster_label)
      }else{
        newList<-umap_embedding
      }
    }else{
      if(return_cluster){
        newList<-list("combined_embedding"=combined_embedding,
          "cluster_label"=cluster_label)
        
      }else{
        newList<-combined_embedding
      }
    }
  }
  if(verbose){
    cat('\n')
    print("Finished...")
    setTxtProgressBar(pb = pb, value = 20)
  }
  return(newList)
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
#' @importFrom reticulate py_run_string import import_main py_get_attr py_module_available py_set_seed
#' @import glue
#' @export
#'
#' @examples
#' a<-1
#'
RunAdaUMAP<-function(
  X,
  metric = 'euclidean',
  metric_count = 1,
  py_envir = parent.frame(),
  n.neighbors = 30L,
  n.components = 2L,
  n.epochs = NULL,
  learning.rate = 1.0,
  min.dist = 0.3,
  spread = 1.0,
  set.op.mix.ratio = 1.0,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5,
  a = 1.8956, 
  b = 0.8006,
  uwot.sgd = FALSE,
  seed.use = 42L,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  verbose = FALSE){
  if (!py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip in command line(e.g. pip install umap-learn) or through reticulate package in R (e.g. reticulate::py_install('umap') )")
  }
  py_func_names<-c("adaptive_dist_grad",
    "adaptive_dist2_grad",
    "adaptive_dist3_grad",
    "adaptive_dist_general_grad")
  
  # source the python script into the main python module
  py_run_string(glue(
    "
import numba 
import numpy as np
import warnings

from umap import distances as dist
py_metric='{metric}' 
py_dist = dist.named_distances_with_gradients[py_metric]
warnings.filterwarnings('ignore')
@numba.njit(fastmath=True)
def adaptive_dist_grad(x, y):
    result = 0.0
    npcs = int((len(x)-1)/2)
    if x[0] != y[0]:
        d,grad = py_dist(x[1:(npcs+1)],y[1:(npcs+1)])
    else:
        d,grad = py_dist(x[(npcs+1):(2*npcs+1)],y[(npcs+1):(2*npcs+1)])
    return d, grad

@numba.njit(fastmath=True)
def adaptive_dist2_grad(x, y):
    result = 0.0
    npcs = int((len(x)-2)/3)
    if x[0] != y[0]:
        d,grad = py_dist(x[2:(npcs+2)],y[2:(npcs+2)])
    else:
        if x[1] != y[1] or x[1] == -1:
            d,grad = py_dist(x[(npcs+2):(2*npcs+2)],y[(npcs+2):(2*npcs+2)])
        else:
            d,grad = py_dist(x[(2*npcs+2):(3*npcs+2)],y[(2*npcs+2):(3*npcs+2)])
    return d, grad 

@numba.njit(fastmath=True)
def adaptive_dist3_grad(x, y):
    
    result = 0.0
    npcs = int((len(x)-3)/4)
    if x[0] != y[0]:
        d,grad = py_dist(x[3:(npcs+3)],y[3:(npcs+3)])
    else:
        if x[1] != y[1] or x[1] == -1:
            d,grad = py_dist(x[(npcs+3):(2*npcs+3)],y[(npcs+3):(2*npcs+3)])
        else:
            if x[2] != y[2] or x[2] == -1:
                d,grad = py_dist(x[(2*npcs+3):(3*npcs+3)],y[(2*npcs+3):(3*npcs+3)])
            else:
                d,grad = py_dist(x[(3*npcs+3):(4*npcs+3)],y[(3*npcs+3):(4*npcs+3)])
    return d, grad

@numba.njit(fastmath=True)
def adaptive_dist_general_grad(x, y):
    result = 0.0
    num_layer = x[0]
    npcs = int((len(x)-num_layer)/num_layer)
    processed = False
    for layer in range(1,num_layer):
        if x[layer] != y[layer]:
            print(layer)
            d,grad = py_dist(x[((layer-1)*npcs+num_layer):(layer*npcs+num_layer)],y[((layer-1)*npcs+num_layer):(layer*npcs+num_layer)])
            processed = True
            break
    if not processed:
        d,grad=py_dist(x[((num_layer-1)*npcs+num_layer):(num_layer*npcs+num_layer)],y[((num_layer-1)*npcs+num_layer):(num_layer*npcs+num_layer)])
         
    return d, grad
")  
    ,local = FALSE, convert = TRUE)
  
  
  # copy objects from the main python module into the specified R environment
  py_main <- import_main(convert = TRUE)
  py_main_dict <- py_get_attr(py_main, "__dict__")
  
  Encoding(py_func_names) <- "UTF-8"
  for (py_name in py_func_names){
    py_value <- py_main_dict[[py_name]]
    assign(py_name, py_value, envir = py_envir) 
  }
  if (!is.null(x = seed.use)) {
    py_set_seed(seed = seed.use)
  }
  if (typeof(x = n.epochs) == "double") {
    n.epochs <- as.integer(x = n.epochs)
  }
  
  if (metric_count == 1){
    adaptive_metric<-adaptive_dist_grad
  }
  if (metric_count == 2){
    adaptive_metric<-adaptive_dist2_grad
  }
  if (metric_count == 3){
    adaptive_metric<-adaptive_dist3_grad
  }
  if (metric_count == 4){
    adaptive_metric<-adaptive_dist_general_grad
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



#


get_umap_embedding_adjust<-function(
  pca_embedding,
  pca_center,
  pca_anchor_index,
  pca_dist,
  umap_embedding,
  N_label,
  cluster_,
  label_index,
  adjust_method = "umap",
  distance_metric = "euclidean",
  scale_factor =1,
  rotate = TRUE,
  seed.use = 42
){
  rotation = function(x,y){
    u=x/sqrt(sum(x^2))
    
    v=y-sum(u*y)*u
    v=v/sqrt(sum(v^2))
    
    cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))
    
    sint=sqrt(1-cost^2);
    
    diag(length(x)) - u %*% t(u) - v %*% t(v) + 
      cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v))
  }
  Rotation2to1<-function(umap_center1,umap_center2,pos){
    N_label<-nrow(umap_center1)
    umap_center1_tmp<-t(t(umap_center1[-pos,])-as.numeric(umap_center1[pos,]))
    weight_1<-1/pdist(umap_center1_tmp,c(0,0))@dist
    weight_1<-weight_1/sum(weight_1)
    umap_center2_tmp<-t(t(umap_center2[-pos,])-as.numeric(umap_center2[pos,]))
    #weight_2<-1/pdist(umap_center2_tmp,c(0,0))@dist
    #weight_2<-weight_2/sum(weight_2)
    angles<-sapply(1:(N_label-1), function(i){
      umap1<-umap_center1_tmp[i,]
      umap2<-umap_center2_tmp[i,]
      umap1<-umap1/sqrt(sum(umap1^2))
      umap2<-umap2/sqrt(sum(umap2^2))
      Rumap2toumap1<-rotation(umap2,umap1)
      angle<-acos(Rumap2toumap1[1,1])
      if(Rumap2toumap1[2,1]>=0){
        angle<-acos(Rumap2toumap1[1,1])
      }else{
        angle<- -acos(Rumap2toumap1[1,1])
      }
      angle
    })
    
    #angle2to1<-mean(angles)
    #angle2to1<-median(angles)
    angle2to1<-sum(angles*weight_1)
    R2to1<-diag(cos(angle2to1),2)
    R2to1[2,1]<-sin(angle2to1)
    R2to1[1,2]<- -R2to1[2,1]
    R2to1
  }
  if (adjust_method == "umap"){
    set.seed(seed.use)
    umap_center <-
      umap(
        X = pca_dist,
        n_neighbors = as.integer(x = N_label-1),
        n_components = as.integer(x =2L),
        metric = distance_metric,
        learning_rate = 1.0,
        min_dist = 0.3,
        spread =  1.0,
        set_op_mix_ratio =  1.0,
        local_connectivity =  max(1,N_label-1),
        repulsion_strength = 1,
        negative_sample_rate = 5,
        fast_sgd = FALSE
      )
  }else if (adjust_method == "MDS"){
    umap_center<-isoMDS(pca_dist)$points 
  }else{
    stop("wrong adjust method")
  }
  colnames(umap_center)<-c("UMAP_1","UMAP_2")
  umap_center<-data.frame(umap_center)
  sf1<-(max(umap_embedding[,1])-min(umap_embedding[,1]))/(max(umap_center[,1]) -min(umap_center[,1]))
  sf2<-(max(umap_embedding[,2])-min(umap_embedding[,2]))/(max(umap_center[,2]) -min(umap_center[,2]))
  umap_center[,1]<-umap_center[,1]*sf1*scale_factor
  umap_center[,2]<-umap_center[,2]*sf2*scale_factor
  umap_embedding_mean<-t(sapply(1:N_label, function(i){
    index_i<-which(cluster_ == label_index[i])
    colMeans(as.matrix(umap_embedding[index_i,]))
  }))
  umap_embedding_adjust<-umap_embedding
  
  if (rotate){
    
    umap_center_flip<-umap_center
    umap_center_flip[,1]<-umap_center_flip[,1]*(-1)
    angle_var<-c()
    weight_sample<-c()
    angle_no_flip<-list()
    angle_flip<-list()
    for (i in 1:N_label){
      weight_sample<-c(weight_sample,sum(cluster_ == label_index[i]))
      R2to1<-Rotation2to1(
        umap_center1 = umap_center,
        umap_center2 = pca_center[,c(1,2)],
        pos = i)  
      
      R2to1_flip<-Rotation2to1(
        umap_center1 = umap_center_flip,
        umap_center2 = pca_center[,c(1,2)],
        pos = i) 
      
      angle_vec<-sapply(1:length(pca_anchor_index[[i]]), function(j){
        y<-as.numeric(pca_embedding[pca_anchor_index[[i]][j],c(1,2)]-pca_center[i,c(1,2)])%*%R2to1
        y<-as.numeric(y)
        y<-y/sqrt(sum(y^2))
        x<-umap_embedding[pca_anchor_index[[i]][j],]-umap_embedding_mean[i,]
        x<-as.numeric(x)
        y<-y/sqrt(sum(y^2))
        x<-x/sqrt(sum(x^2))
        
        Rx2y <- rotation(x,y)
        Rx2y <- pmax(Rx2y,-1)
        Rx2y <- pmin(Rx2y,1)
        if(Rx2y[2,1]>=0){
          i
          angle<-acos(Rx2y[1,1])
        }else{
          angle<- - acos(Rx2y[1,1])
        }
        angle
      })
      angle_vec_flip<-sapply(1:length(pca_anchor_index[[i]]), function(j){
        y<-as.numeric(pca_embedding[pca_anchor_index[[i]][j],c(1,2)]-pca_center[i,c(1,2)])%*%R2to1_flip
        y<-as.numeric(y)
        y<-y/sqrt(sum(y^2))
        x<-umap_embedding[pca_anchor_index[[i]][j],]-umap_embedding_mean[i,]
        x<-as.numeric(x)
        y<-y/sqrt(sum(y^2))
        x<-x/sqrt(sum(x^2))
        
        Rx2y <- rotation(x,y)
        
        
        if(Rx2y[2,1]>=0){
          angle<- acos(Rx2y[1,1])
        }else{
          angle<- -1*acos(Rx2y[1,1])
        }
        angle
      })
      
      angle_vec<-angle_vec[angle_vec>min(angle_vec)]
      angle_vec<-angle_vec[angle_vec<max(angle_vec)]
      angle_vec_flip<-angle_vec_flip[angle_vec_flip>min(angle_vec_flip)]
      angle_vec_flip<-angle_vec_flip[angle_vec_flip<max(angle_vec_flip)]
      angle_no_flip[[i]]<- angle_vec
      angle_flip[[i]] <- angle_vec_flip
      angle_var<-rbind(angle_var,(c(var(angle_vec),var(angle_vec_flip))))
    }
    weight_sample<-weight_sample/sum(weight_sample)
    if (sum(angle_var[,1]*weight_sample)>=sum(angle_var[,2]*weight_sample)){
      # use flip
      angle_vec<-angle_flip
    }else{
      # use no flip
      angle_vec<-angle_no_flip
    }
    
    for ( i in 1:N_label){
      anglex2y<-mean(angle_vec[[i]])
      #print(anglex2y)
      Rx2y<-diag(cos(anglex2y),2)
      Rx2y[2,1]<-sin(anglex2y)
      Rx2y[1,2]<- -Rx2y[2,1]
      index_i<-which(cluster_ == label_index[i])
      umap_embedding_adjust[index_i,]<-t(t(t(t(umap_embedding[index_i,])-as.numeric(umap_embedding_mean[i,]))%*%Rx2y)+as.numeric(umap_center[i,]))
    }
    
  }else{
    
    for(i in 1:N_label){
      index_i<-which(cluster_ == label_index[i])
      umap_embedding_adjust[index_i,]<-t(t(umap_embedding[index_i,])-as.numeric(umap_embedding_mean[i,]-umap_center[i,]))
    }
  }
  newList<-list("umap_center"=umap_center,
    "umap_embedding_adjust"=umap_embedding_adjust)
  return(newList)
}


#' adjustUMAP
#'
#' Adjust UMAP to deal with distortion 
#'
#' @details amazing umap
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom pdist pdist
#' @importFrom uwot umap
#' @export
#'
#' @examples
#' a<-1
#'
#'
adjustUMAP<-function(
  pca_embedding,
  umap_embedding,
  global_umap_embedding = NULL,
  adjust_method = "umap",
  distance_metric = "euclidean",
  scale_factor =1,
  rotate = TRUE,
  seed.use = 42,
  min_size = 100,
  maxit_push = NULL
){
  if(is.null(global_umap_embedding)){
    global_umap_embedding<-umap_embedding
  }
  snn_<- FindNeighbors(object = umap_embedding,
    verbose = F)$snn
  cluster_ <- FindClusters(snn_,
    resolution = 0,
    verbose = F)[[1]]
  cluster_ <- as.numeric(as.character(cluster_))
  
  N_label<-length(unique(cluster_))
  label_index<-sort(unique(cluster_))
  if(N_label <= 2){
    return(umap_embedding)
  }
  cluster_size<-sapply(1:N_label, function(i){
    sum(cluster_==label_index[i])
  })
  N_sample<-nrow(pca_embedding)
  
  prop_density<-sapply(1:N_label, function(i){
    index_i<-which(cluster_ == label_index[i])
    set.seed(seed.use)
    sample_index_i<-sample(index_i,min(min_size,length(index_i)) )
    sample_global_dist<-Dist(global_umap_embedding[sample_index_i,])
    sample_local_dist<-Dist(umap_embedding[sample_index_i,])
    mean(c(sample_global_dist))/mean(c(sample_local_dist))
  })
  #print(prop_density)
  #for(i in 1:N_label){
  #  index_i<-which(cluster_ == label_index[i])
  #  cur_umap<-umap_embedding[index_i,]
  #  umap_embedding[index_i,]<-t((t(cur_umap)-as.numeric(colMeans(cur_umap)))*min(2,prop_density[i])+as.numeric(colMeans(cur_umap)))
  #}
  
  pca_center<-t(sapply(1:N_label, function(i){
    
    index_i<-which(cluster_ == label_index[i])
    
    colMeans(as.matrix(pca_embedding[index_i,]))
  }))
  
  
  pca_anchor_index<-lapply(1:N_label, function(i){
    #print(label_index[i])
    index_i<-which(cluster_ == label_index[i])
    set.seed(seed.use)
    sample_index_i<-sample(index_i,min(min_size,length(index_i)) )
    sample_index_dist<-pdist(pca_embedding[sample_index_i,c(1,2)],pca_center[i,c(1,2)])@dist
    nth(x=sample_index_dist,
      k = max(1,min(ceiling(min_size/5),length(index_i))),
      num.of.nths = 2,
      descending = T,
      index.return = T)[,1]
  })
  pca_dist1<-Dist(pca_center)
  
  index_<-which(cluster_size < N_sample*0.1)
  #prop_<-exp(cluster_size/max(cluster_size))/
  #  max(exp(cluster_size/max(cluster_size)))
  #for(i in 1:N_label){
  #  pca_dist1[,i]<-pca_dist1[,i]*prop_[i]
  #  pca_dist1[i,]<-pca_dist1[i,]*prop_[i]
  #}
  
  #pca_dist_vec<-c(pca_dist1[pca_dist1>0])
  #pca_dist_cluster<-Ckmeans.1d.dp(pca_dist_vec,k = 2)
  
  #if(max(pca_dist_vec[pca_dist_cluster$cluster==1])*1.5<
  #    min(pca_dist_vec[pca_dist_cluster$cluster==2])){
  
  #  pca_dist1[pca_dist1 > min(pca_dist_vec[pca_dist_cluster$cluster==2])]<- 
  #    pca_dist1[pca_dist1 > min(pca_dist_vec[pca_dist_cluster$cluster==2])] +
  #    (max(pca_dist_vec[pca_dist_cluster$cluster==1])-
  #    min(pca_dist_vec[pca_dist_cluster$cluster==2]))
  #}
  pca_dist<-as.dist(pca_dist1)
  step1_res<-get_umap_embedding_adjust(
    pca_embedding=pca_embedding,
    pca_center=pca_center,
    pca_anchor_index=pca_anchor_index,
    pca_dist=pca_dist,
    adjust_method = adjust_method,
    distance_metric=distance_metric,
    umap_embedding=umap_embedding,
    N_label=N_label,
    cluster_ = cluster_,
    label_index = label_index,
    scale_factor =scale_factor,
    rotate = rotate,
    seed.use = seed.use)
  umap_center<-step1_res$umap_center
  umap_embedding_adjust<-step1_res$umap_embedding_adjust
  library(Seurat)
  snn_1<- FindNeighbors(object = umap_embedding_adjust,
    verbose = F)$snn
  cluster_1 <- FindClusters(snn_1,
    resolution = 0,
    verbose = F)[[1]]
  cluster_1 <- as.numeric(as.character(cluster_1))
  
  N_label1<-length(unique(cluster_1))
  if ( N_label1 != N_label){
    ## First Adjustment
    label_index1<-sort(unique(cluster_1))
    
    bad_index<-list()
    list_index<-0
    for ( i in 1:N_label1){
      index_i1<-which(cluster_1 == label_index1[i])
      cur_index_len<-length(unique(cluster_[index_i1]))
      if(cur_index_len > 1){
        list_index<-list_index+1
        bad_index[[list_index]]<-c(unique(cluster_[index_i1]))
      }
    }
    for (i in 1:length(bad_index)){
      pos<-min(bad_index[[i]])
      other_pos<-bad_index[[i]][bad_index[[i]]>pos]
      pos<-pos+1
      other_pos<-other_pos+1
      dist_vec<-pca_dist1[,pos]
      all_pos<-dist_vec <= max(dist_vec[other_pos])
      pca_dist1[all_pos,pos]<-min(dist_vec[dist_vec>max(dist_vec[other_pos])])
      pca_dist1[pos,all_pos]<-pca_dist1[all_pos,pos]
    }
    pca_dist<-as.dist(pca_dist1)
    pca_dist1<-as.matrix(pca_dist)
    step2_res<-get_umap_embedding_adjust(
      pca_embedding=pca_embedding,
      pca_center=pca_center,
      pca_anchor_index=pca_anchor_index,
      pca_dist=pca_dist,
      adjust_method = adjust_method,
      distance_metric=distance_metric,
      umap_embedding=umap_embedding,
      N_label=N_label,
      cluster_ = cluster_,
      label_index = label_index,
      scale_factor =scale_factor,
      rotate = rotate,
      seed.use = seed.use)
    umap_center<-step2_res$umap_center
    umap_embedding_adjust<-step2_res$umap_embedding_adjust
    snn_2<- FindNeighbors(object = umap_embedding_adjust,
      verbose = F)$snn
    cluster_2 <- FindClusters(snn_2,
      resolution = 0,
      verbose = F)[[1]]
    cluster_2 <- as.numeric(as.character(cluster_2))
    
    N_label2<-length(unique(cluster_2))
    if(N_label2 != N_label){
      ## Second Adjustment Push Away
      label_index2<-sort(unique(cluster_2))
      
      bad_index<-list()
      list_index<-0
      for ( i in 1:N_label2){
        index_i2<-which(cluster_2 == label_index2[i])
        cur_index_len<-length(unique(cluster_[index_i2]))
        if(cur_index_len > 1){
          list_index<-list_index+1
          bad_index[[list_index]]<-c(unique(cluster_[index_i2]))
        }
      }
      N_label3<-N_label - 1
      cur_iter<-0
      if(is.null(maxit_push)){
        maxit_push<-N_label
      }
      while (N_label3 != N_label & cur_iter < maxit_push) {
        cur_iter<-cur_iter+1
        #print(cur_iter)
        for (i in 1:length(bad_index)){
          pos<-min(bad_index[[i]])
          other_pos<-bad_index[[i]][bad_index[[i]]>pos]
          pos<-pos+1
          other_pos<-other_pos+1
          dist_mat<-pdist(umap_center[pos,],umap_center)@dist
          target_distance<-min(dist_mat[dist_mat>max(dist_mat[other_pos])])
          target_distance<-(target_distance+dist_mat[other_pos])/2
          for (k in 1:length(other_pos)){
            cur_pos<-other_pos[k]
            cur_center<-umap_center[pos,]
            cur_arrow<-umap_center[cur_pos,]
            prop<-target_distance[k]/dist_mat[cur_pos]
            cur_arrow<-(cur_arrow-cur_center)*prop+cur_center
            index_i<-which(cluster_ == label_index[cur_pos])
            umap_embedding_adjust[index_i,]<-t(t(umap_embedding_adjust[index_i,])-as.numeric(umap_center[cur_pos,])+as.numeric(cur_arrow))
          }
        }
        
        snn_3<- FindNeighbors(object = umap_embedding_adjust,
          verbose = F)$snn
        cluster_3 <- FindClusters(snn_3,
          resolution = 0,
          verbose = F)[[1]]
        cluster_3 <- as.numeric(as.character(cluster_3))
        N_label3<-length(unique(cluster_3))
        label_index3<-sort(unique(cluster_3))
        
        bad_index<-list()
        list_index<-0
        for ( i in 1:N_label3){
          index_i3<-which(cluster_3 == label_index3[i])
          cur_index_len<-length(unique(cluster_[index_i3]))
          if(cur_index_len > 1){
            list_index<-list_index+1
            bad_index[[list_index]]<-c(unique(cluster_[index_i3]))
          }
        }
        
      }
      
      if(N_label3 != N_label){
        print("Use Third")
        ## Third Adjustment Rescale
        for (i in 1:length(bad_index)){
          pos<-min(bad_index[[i]])
          other_pos<-bad_index[[i]][bad_index[[i]]>pos]
          pos<-pos+1
          other_pos<-other_pos+1
          dist_mat<-pdist(umap_center[pos,],umap_center)@dist
          target_distance<-min(dist_mat[dist_mat>max(dist_mat[other_pos])])
          re_sf<-target_distance/max(dist_mat[other_pos])
          
          umap_embedding_adjust<-umap_embedding_adjust*re_sf
        }
      }
      
    }}
  return(umap_embedding_adjust)
  snn_<- FindNeighbors(object = umap_embedding_adjust,
    verbose = F)$snn
  cluster_ <- FindClusters(snn_,
    resolution = 0,
    verbose = F)[[1]]
  cluster_ <- as.numeric(as.character(cluster_))
  
  N_label<-length(unique(cluster_))
  label_index<-sort(unique(cluster_))
  
  adjust_umap_center<-t(sapply(1:N_label, function(i){
    index_i<-which(cluster_ == label_index[i])
    colMeans(as.matrix(umap_embedding_adjust[index_i,]))
  }))
  
  umap_dist1<-Dist(adjust_umap_center)
  diag(umap_dist1)<-Inf
  min_dist_vec<-colMins(umap_dist1,value = T)
  update_prop<-1/2
  umap_embedding_adjust4<-umap_embedding_adjust
  if(max(min_dist_vec)<= max(min_dist_vec[min_dist_vec<max(min_dist_vec)])*1.5){
    return(umap_embedding_adjust)
  }else{
    update_index<-which.max(min_dist_vec)[1]
    closed_index<-which.min(umap_dist1[,update_index])
    new_center<-umap_center[update_index,] + update_prop*(umap_center[closed_index,]-umap_center[update_index,])
    index_i<-which(cluster_ == label_index[update_index])
    umap_embedding_adjust4[index_i,]<-t(t(t(t(umap_embedding_adjust[index_i,])-as.numeric(adjust_umap_center[update_index,])))+as.numeric(new_center))
    
    
    snn_4<- FindNeighbors(object = umap_embedding_adjust4,
      verbose = F)$snn
    cluster_4 <- FindClusters(snn_4,
      resolution = 0,
      verbose = F)[[1]]
    cluster_4 <- as.numeric(as.character(cluster_4))
    
    N_label4<-length(unique(cluster_4))
    if(N_label4 != N_label){
      ## Fourth Adjustment
      while (N_label4 != N_label) {
        update_prop<-update_prop*2/3
        umap_embedding_adjust4<-umap_embedding_adjust
        if(max(min_dist_vec)> max(min_dist_vec[min_dist_vec<max(min_dist_vec)])*1.5){
          update_index<-which.max(min_dist_vec)[1]
          closed_index<-which.min(umap_dist1[,update_index])
          new_center<-umap_center[update_index,] + update_prop*(umap_center[closed_index,]-umap_center[update_index,])
          index_i<-which(cluster_ == label_index[update_index])
          umap_embedding_adjust4[index_i,]<-t(t(t(t(umap_embedding_adjust[index_i,])-as.numeric(adjust_umap_center[update_index,])))+as.numeric(new_center))
        }
        
        snn_4<- FindNeighbors(object = umap_embedding_adjust4,
          verbose = F)$snn
        cluster_4 <- FindClusters(snn_4,
          resolution = 0,
          verbose = F)[[1]]
        cluster_4 <- as.numeric(as.character(cluster_4))
        N_label4<-length(unique(cluster_4))
      }
    }
    adjust_umap_center<<-t(sapply(1:N_label, function(i){
      index_i<-which(cluster_ == label_index[i])
      colMeans(as.matrix(umap_embedding_adjust4[index_i,]))
    }))
    return(umap_embedding_adjust4)
  }
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
  scale_factor_separation =3,
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
    min(cluster_center_dist[i,],na.rm = T)/scale_factor_separation/cluster_semi_d[i]
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
  npcs=10,
  nfeatures=2000,
  process_min_size=0
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
        center_method = "mean",
        scale_factor_separation=scale_factor_separation)$combined_embedding
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
    sub_PC_list = sub_PC_list_i)
  
  
  return(combined_embedding)
}

#' FormCombinedEmbedding
#'
#' FormCombinedEmbedding
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
FormCombinedEmbedding<-function(
  expr_matrix,
  expr_matrix_pca,
  cluster_label,
  npcs=10,
  nfeatures =2000,
  center_method = "mean",
  scale_factor_separation=3
){
  sub_PC_list<-SubPCEmbedding(
    expr_matrix = expr_matrix,
    cluster_label = cluster_label,
    npcs = npcs,
    nfeatures = nfeatures)
  combined_embedding<-CombinePC(
    PC_embedding = expr_matrix_pca,
    cluster_label = cluster_label,
    sub_PC_list = sub_PC_list)
  combined_embedding<-ScaleFactor(
    combined_embedding = combined_embedding,
    cluster_label = cluster_label,
    npcs = npcs,
    center_method = center_method,
    scale_factor_separation = scale_factor_separation)$combined_embedding
  combined_embedding<-as.matrix(combined_embedding)
  return(combined_embedding)
}


#' FormAdaptiveCombineList
#'
#' FormAdaptiveCombineList
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
#'
FormAdaptiveCombineList<-function(
  expr_matrix,
  expr_matrix_pca,
  max_stratification,
  stratification_count,
  scale_factor_separation,
  resolution,
  cluster_method,
  npcs,
  nfeatures,
  process_min_size,
  differentail_gene_cutoff = 20,
  do_cluster = TRUE,
  cluster_label = NULL,
  check_differential =TRUE,
  verbose = FALSE
){
  colnames(expr_matrix_pca)<-paste0("Layer",stratification_count,"PC",1:npcs)
  if(nrow(expr_matrix_pca) < process_min_size){
    newList<-list("cluster_label"= -1,
      "combined_embedding"=expr_matrix_pca)
    return(newList)
  }
  if(nrow(expr_matrix_pca)!=ncol(expr_matrix)){
    stop("expr_matrix_pca and expr_matrix do not match")
  }
  N_gene<-nrow(expr_matrix)
  if(is.null(cluster_label[1])){
    do_cluster<-TRUE
  }
  if(do_cluster){
    cluster_label<-DoCluster(
      pc_embedding = expr_matrix_pca,
      method = cluster_method,
      resolution = resolution)$cluster
  }
  # sorted unique cluster label
  label_index<-sort(as.numeric(
    unique(as.character(cluster_label))))
  # number of cluster label
  N_label<-length(label_index)
  
  ## Compare npcs with the cluster size
  size_cluster<-c()
  for ( i in 1:N_label){
    size_cluster<-c(size_cluster,
      sum(cluster_label == label_index[i]))
  }
  if(sum(size_cluster <= npcs) > 0){
    warning("Produce a cluster whose size is less than npcs: Combine it with other cluster, Or please adjust the npcs or resolution")
    cur_index<-which(size_cluster <= npcs)
    pca_center<-t(sapply(1:N_label, function(i){
      index_i<-which(cluster_label == label_index[i])
      colMeans(as.matrix(expr_matrix_pca[index_i,]))
    }))
    sample_index_dist<-pdist(pca_center[-cur_index,],pca_center[cur_index,])@dist
    sample_index_dist<-matrix(sample_index_dist,nrow = sum(size_cluster <= npcs))
    comb_list<-sapply( 1:length(cur_index), function(i){
      c(cur_index[i],Rfast::nth(x=sample_index_dist[i,],
        k = 1,
        num.of.nths = 2,
        descending = F,
        index.return = T)[,1]) 
    })
    for( i in 1:ncol(comb_list)){
      index_1<-which(cluster_label == label_index[comb_list[1,i]])
      cluster_label[index_1]<-label_index[comb_list[2,i]]
    }
    # update label_index
    # sorted unique cluster label
    label_index<-sort(as.numeric(
      unique(as.character(cluster_label))))
    # number of cluster label
    N_label<-length(label_index)
    
    ## Compare npcs with the cluster size
    size_cluster<-c()
    for ( i in 1:N_label){
      size_cluster<-c(size_cluster,
        sum(cluster_label == label_index[i]))
    }
  }
  
  if(check_differential){
    if(N_label == 1){
      newList<-list("cluster_label"= -1,
        "combined_embedding"=expr_matrix_pca)
      return(newList)
    }else{
      ### Begin of else 
      cur_label<- N_label
      while(cur_label>1){
        for (i in 1:(cur_label-1)){
          print(paste0("Current:",label_index[i],"vs",label_index[cur_label]))
          index_1<-which(cluster_label == label_index[i])
          index_2<-which(cluster_label == label_index[cur_label])
          #S<-cbind(expr_matrix[,index_1],expr_matrix[,index_2])
          #S<-CreateSeuratObject(S)
          #label_diff<-c(rep(0,length(index_1)),rep(1,length(index_2)))
          #Idents(S)<-factor(label_diff)
          #.<-capture.output(marker_diff<-FindMarkers(S,
          #  test.use = "t",
          #  ident.1 = unique(S@active.ident)[1],
          # ident.2 =  unique(S@active.ident)[2]))
          #print(sum(marker_diff[,5]<0.05))
          #if (sum(marker_diff[,5]<=0.05) < 5){
          #print(paste0("Don't find differential genes between",label_index[i],"vs",label_index[cur_label]))
          #  cluster_label[index_2]<-label_index[i]
          #  label_index<-sort(as.numeric(
          #    unique(as.character(cluster_label))))
          #  N_label<-length(label_index)
          #  break
          #}
          cur_differential_gene<-0
          cur_gene<-1
          while (cur_differential_gene < differentail_gene_cutoff & 
              cur_gene <= N_gene){
            t_res<-t.test(expr_matrix[cur_gene,index_1],
              expr_matrix[cur_gene,index_2])
            t_rest_res<<-t_res
            if (is.na(t_res$p.value)){
              t_res$p.value<-1
            }
            if (t_res$p.value <0.05){
              cur_differential_gene <- cur_differential_gene+1
            }
            cur_gene<-cur_gene+1
          }
          if(cur_differential_gene < differentail_gene_cutoff){
            print(paste0("Don't find differential genes between",label_index[i],"vs",label_index[cur_label]))
            cluster_label[index_2]<-label_index[i]
            label_index<-sort(as.numeric(
              unique(as.character(cluster_label))))
            N_label<-length(label_index)
            break
          }
        }
        cur_label<-cur_label - 1
        
      }
      ### End of else 
    }
    
    
  }
  
  if(N_label == 1){
    newList<-list("cluster_label"= -1,
      "combined_embedding"=expr_matrix_pca)
    return(newList)
  }
  global_matrix<<-expr_matrix
  global_pca<<-expr_matrix_pca
  global_cluster<<-cluster_label
  
  combined_embedding<-FormCombinedEmbedding(
    expr_matrix=expr_matrix,
    expr_matrix_pca=expr_matrix_pca,
    cluster_label=cluster_label,
    npcs=npcs,
    nfeatures=nfeatures,
    scale_factor_separation=scale_factor_separation
  )
  if(max_stratification == stratification_count + 1){
    newList<-list("cluster_label"= cluster_label,
      "combined_embedding"=combined_embedding)
    return(newList)
  }  
  colnames(combined_embedding)<-c(paste0("Layer",stratification_count,"PC",1:npcs),
    paste0("Layer",(stratification_count+1),"PC",1:npcs))
  
  cluster_label_list<-list()
  combined_embedding_list<-list()
  work_adaptive<-FALSE
  max_ncol_sub<-0
  for ( i in 1:N_label){
    index_i<-which(cluster_label == label_index[i])
    tmp<-FormAdaptiveCombineList(
      expr_matrix = expr_matrix[,index_i],
      expr_matrix_pca = combined_embedding[index_i,(npcs+1):(2*npcs)],
      max_stratification = max_stratification,
      stratification_count = stratification_count + 1,
      scale_factor_separation = scale_factor_separation,
      resolution = resolution,
      cluster_method = cluster_method,
      npcs = npcs,
      nfeatures = nfeatures,
      process_min_size = process_min_size,
      check_differential = check_differential,
      verbose = verbose
    )
    cluster_label_list[[i]]<-tmp$cluster_label
    combined_embedding_list[[i]]<-tmp$combined_embedding
    max_ncol_sub<-max(max_ncol_sub,ncol(tmp$combined_embedding))
    if (sum(tmp$cluster_label!= -1) != 0){
      work_adaptive<-TRUE
    }else{
      if(ncol(tmp$combined_embedding)!=npcs){
        print(ncol(tmp$combined_embedding))
        stop("combined_embedding size is wrong")
      }  
    }
    rm(tmp)
  }
  if (work_adaptive){
    
    if(max_ncol_sub < (2*npcs)){
      stop("combined_embedding size is wrong1")
    }else if (max_ncol_sub == (2*npcs)){
      ## cluster label list is vector
      cluster_label_sub<-rep(-1,length(cluster_label))
      combined_embedding_sub<-matrix(0,
        nrow = nrow(expr_matrix_pca),
        ncol = max_ncol_sub)
      colnames(combined_embedding_sub)<-c(paste0("Layer",(stratification_count+1),"PC",1:npcs),
        paste0("Layer",(stratification_count+2),"PC",1:npcs))
      for ( i in 1:N_label){
        index_i<-which(cluster_label == label_index[i])
        cluster_label_sub[index_i]<-cluster_label_list[[i]]
        
        combined_embedding_sub[index_i,
          1:ncol(combined_embedding_list[[i]])]<-
          combined_embedding_list[[i]]
      }
      newList<-list("cluster_label"= cbind(cluster_label,cluster_label_sub),
        "combined_embedding"=cbind(expr_matrix_pca,combined_embedding_sub))
      return(newList)
      
    }else{
      ## cluster label list is matrix
      max_label_count<-0
      for ( i in 1:N_label){
        if(is.null(dim(cluster_label_list[[i]]))){
          max_label_count<-max(max_label_count,1 )
        }else{
          max_label_count<-max(max_label_count,ncol(cluster_label_list[[i]]))
        }
      }
      cluster_label_sub<-matrix(-1,
        nrow = length(cluster_label),
        ncol = max_label_count)
      combined_embedding_sub<-matrix(0,
        nrow = nrow(expr_matrix_pca),
        ncol = max_ncol_sub)
      
      for ( i in 1:N_label){
        index_i<-which(cluster_label == label_index[i])
        if(is.null(dim(cluster_label_list[[i]]))){
          cluster_label_sub[index_i,1]<-cluster_label_list[[i]]
        }else{
          cluster_label_sub[index_i,1:ncol(cluster_label_list[[i]])
            ]<-cluster_label_list[[i]]
        }
        combined_embedding_sub[index_i,
          1:ncol(combined_embedding_list[[i]])]<-
          combined_embedding_list[[i]]
        if(ncol(combined_embedding_list[[i]]) == max_ncol_sub){
          colnames(combined_embedding_sub)<-colnames(combined_embedding_list[[i]])
        }
      }
      newList<-list("cluster_label"= cbind(cluster_label,cluster_label_sub),
        "combined_embedding"=cbind(expr_matrix_pca,combined_embedding_sub))
      return(newList)
    }
  }else{
    
    newList<-list("cluster_label"= cluster_label,
      "combined_embedding"=combined_embedding)
    return(newList)
  }
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



#################################################################################
######## This is some extra functions
newUmapPlot<-function(expr_pca_combined,
  cell_label= NULL,
  cluster_ = FALSE,
  cluster_label=NULL,
  label_legend =TRUE,
  min.dist = 0.01,
  pt.size=0){
  set.seed(42)
  combinedPC_umap <-
    uwot::umap(
      X = expr_pca_combined,
      a = 1.8956, 
      b = 0.8006, 
      metric = "cosine",
      min_dist = min.dist,
      approx_pow = TRUE, 
    )
  
  combinedPC_umap<-
    data.frame(combinedPC_umap)
  
  colnames(combinedPC_umap) <- 
    c("UMAP_1","UMAP_2")
  
  if(is.null(cell_label)){
    return(combinedPC_umap)
  }
  combinedPC_umap$label<-factor(cell_label)
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category 
    == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, 
    qual_col_pals$maxcolors,rownames(qual_col_pals)))
  
  
  gg<-ggplot(combinedPC_umap)+
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
        col_vector[c(1:length(unique(cell_label)))])
  #+theme(legend.position="none")
  if(!label_legend){
    gg<-gg+theme(legend.position="none")
  }else{
    centers<-dplyr::group_by(combinedPC_umap,label) 
    centers<-dplyr::summarize(centers,x = median(x = UMAP_1),
      y = median(x = UMAP_2),.groups = 'drop')
    
    gg <- gg +
      ggrepel::geom_text_repel(data = centers, 
        mapping = aes(x = x, y = y, 
          label = label), size = 4)
  }
  #combinedPC_umap[,c(1,2,which(colnames(
  #  combinedPC_umap) == "label"))] %>%
  #  dplyr::group_by(label) %>%
  #  summarize(x = median(x = UMAP_1),
  #    y = median(x = UMAP_2)) -> centers
  
  
  if (cluster_){
    combinedPC_umap$cluster_label<-
      factor(cluster_label)
    gg1<-ggplot(combinedPC_umap)+
      geom_point(aes(x = UMAP_1, 
        y = UMAP_2,color = cluster_label),
        size = pt.size)+theme(legend.title =
            element_blank())+
      theme(panel.grid.major = 
          element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = 
            "black"), 
        legend.key=element_blank())+
      guides(color = guide_legend(override.aes = 
          list(size=3)))+
      theme(legend.position="none")
    
    combinedPC_umap[,c(1,2,which(colnames(
      combinedPC_umap) == "cluster_label"))] %>%
      dplyr::group_by(cluster_label) %>%
      summarize(x = mean(x = UMAP_1), 
        y = mean(x = UMAP_2),.groups = 'drop') -> centers
    
    gg1 <- gg1 +
      geom_text(data = centers, 
        mapping = aes(x = x, y = y, 
          label = cluster_label), 
        size = 4)
    newList<-list("umap"=gg,
      "cluster"=gg1,
      "data"=combinedPC_umap[,c(1,2)])
    return(newList)
  }else{
    newList<-list("umap"=gg,
      "data"=combinedPC_umap[,c(1,2)])
    return(newList)
  }
}

RunOrPCA<-function(expr_matrix,count = T,npcs=10,nfeatures=2000){
  if(count){
    expr_matrix<-NormalizeData(
      expr_matrix,
      verbose = F)
  }
  expr_matrix_hvg <- FindVariableFeatures(
    expr_matrix,
    verbose = F)$vst.variance.standardized
  hvg<- nth(x = expr_matrix_hvg,
    k = nfeatures,
    num.of.nths = 2,
    descending = T,
    index.return = T)[,1]
  
  expr_matrix<-expr_matrix[hvg,]
  expr_matrix <- ScaleData(
    expr_matrix,
    verbose = F)
  suppressWarnings(expr_matrix_pca <- RunPCA(
    object = expr_matrix,
    features = rownames(expr_matrix),
    npcs = npcs,
    verbose = F)@cell.embeddings)
  rm(expr_matrix)
  expr_matrix_pca<-data.frame(expr_matrix_pca)
  return(expr_matrix_pca)
}

RunSeuratPCA<-function(expr_matrix,npcs=10,nfeatures=2000){
  expr_matrix <- CreateSeuratObject(counts = expr_matrix)
  
  expr_matrix<-NormalizeData(
    expr_matrix,
    verbose = F)
  expr_matrix <- FindVariableFeatures(expr_matrix,
    nfeatures = nfeatures)
  expr_matrix <- ScaleData(expr_matrix, features = rownames(expr_matrix))
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  suppressWarnings(expr_matrix_pca <- RunPCA(
    object = expr_matrix,
    features = VariableFeatures(object = expr_matrix),
    npcs = npcs,
    verbose = F)@cell.embeddings)
  rm(expr_matrix)
  expr_matrix_pca<-data.frame(expr_matrix_pca)
  return(expr_matrix_pca)
}
