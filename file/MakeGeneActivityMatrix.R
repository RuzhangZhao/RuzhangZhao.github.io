#### MakeGeneActivityMatrix
Extend <- function(
  x,
  upstream = 0,
  downstream = 0,
  from.midpoint = FALSE
) {
  if (any(strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- strand(x = x) == "+" | strand(x = x) == "*"
  if (from.midpoint) {
    midpoints <- start(x = x) + (width(x = x) / 2)
    new_start <- midpoints - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- midpoints + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  } else {
    new_start <- start(x = x) - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- end(x = x) + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  }
  ranges(x = x) <- IRanges(start = new_start, end = new_end)
  x <- trim(x = x)
  return(x)
}

CreateGeneActivityMatrix3 <- function (
  peak.matrix, 
  annotation.file, 
  all.features,
  seq.levels = c(1:22,"X", "Y"), 
  include.body = TRUE, 
  upstream = 2000, 
  downstream = 0, 
  keep.sparse = FALSE, 
  verbose = TRUE) 
{
  #.Deprecated(new = "Signac::GeneActivity", msg = paste("CreateGeneActivityMatrix functionality is being moved to Signac.", 
  #                                                      "Equivalent functionality can be", "achieved via the Signac::GeneActivity function; for", 
  #                                                      "more information on Signac, please see https://github.com/timoast/Signac"))
  # if (!PackageCheck("GenomicRanges", error = FALSE)) {
  #   stop("Please install GenomicRanges from Bioconductor.")
  # }
  # if (!PackageCheck("rtracklayer", error = FALSE)) {
  #   stop("Please install rtracklayer from Bioconductor.")
  # }
  peak.df <- rownames(x = peak.matrix)
  # peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, 
  #                                                           pattern = ":", replacement = "-"), split = "-"))
  peak.df <- do.call(what = rbind, args = strsplit(x = peak.df, split = "-"))
  
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", "start", "end")
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
  BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1
  
  gtf <- annotation.file#rtracklayer::import(con = paste0(dir_data_annotation,"Mus_musculus.GRCm38.75.gtf"))
  #gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, 
  #  pruning.mode = "coarse")
  #if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
  #GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
  #}
  gtf.genes <- gtf#[gtf$type == "gene"]
  if (include.body) {
    suppressWarnings(gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, 
      downstream = downstream))
  }else{
    gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, 
      upstream = upstream, downstream = downstream)
  }
  suppressWarnings(gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom))
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
  gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]
  peak.ids$gene.name <- gene.ids$gene_name
  peak.ids <- as.data.frame(x = peak.ids)
  peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
  annotations <- peak.ids[, c("peak", "gene.name")]
  colnames(x = annotations) <- c("feature", "new_feature")
  if (!keep.sparse) {
    peak.matrix <- as(object = peak.matrix, Class = "matrix")
  }
  
  #if (nbrOfWorkers() > 1) {
  #  mysapply <- future_sapply
  #}
  #else {
  library(pbapply)
    mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
  #}
  all.features<-all.features[all.features%in%annotations$new_feature]
  newmat.list <- mysapply(X = 1:length(x = all.features), FUN = function(x) {
    features.use <- annotations[annotations$new_feature == 
        all.features[x], ]$feature
    submat <- peak.matrix[features.use, ]
    if (length(x = features.use) > 1) {
      submat <- Matrix::colSums(submat)
    }
    if (keep.sparse) {
      return(as(object = as.matrix(x = submat), Class = "dgCMatrix"))
    }else {
      return(as.matrix(x = submat))
    }
  }, simplify = FALSE)
  newmat = do.call(what = cbind, args = newmat.list)
  newmat <- t(x = newmat)
  rownames(x = newmat) <- all.features
  return(as(object = newmat, Class = "dgCMatrix"))
}

CreatePeakActivityMatrix <- function (
  peak.matrix, 
  annotation.file, 
  all.features,
  seq.levels = c(1:22,"X", "Y"), 
  include.body = TRUE, 
  upstream = 2000, 
  downstream = 0, 
  keep.sparse = FALSE, 
  verbose = TRUE) 
{
  #.Deprecated(new = "Signac::GeneActivity", msg = paste("CreateGeneActivityMatrix functionality is being moved to Signac.", 
  #                                                      "Equivalent functionality can be", "achieved via the Signac::GeneActivity function; for", 
  #                                                      "more information on Signac, please see https://github.com/timoast/Signac"))
  # if (!PackageCheck("GenomicRanges", error = FALSE)) {
  #   stop("Please install GenomicRanges from Bioconductor.")
  # }
  # if (!PackageCheck("rtracklayer", error = FALSE)) {
  #   stop("Please install rtracklayer from Bioconductor.")
  # }
  peak.df <- rownames(x = peak.matrix)
  # peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, 
  #                                                           pattern = ":", replacement = "-"), split = "-"))
  peak.df <- do.call(what = rbind, args = strsplit(x = peak.df, split = "-"))
  
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", "start", "end")
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
  BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1
  
  gtf <- annotation.file#rtracklayer::import(con = paste0(dir_data_annotation,"Mus_musculus.GRCm38.75.gtf"))
  #gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, 
  #  pruning.mode = "coarse")
  #if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
  #GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
  #}
  gtf.genes <- gtf#[gtf$type == "gene"]
  if (include.body) {
    suppressWarnings(gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, 
      downstream = downstream))
  }else{
    gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, 
      upstream = upstream, downstream = downstream)
  }
  suppressWarnings(gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom))
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
  gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]
  peak.ids$gene.name <- gene.ids$gene_name
  peak.ids <- as.data.frame(x = peak.ids)
  peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
  annotations <- peak.ids[, c("peak", "gene.name")]
  colnames(x = annotations) <- c("feature", "new_feature")
  if (!keep.sparse) {
    peak.matrix <- as(object = peak.matrix, Class = "matrix")
  }
  
  #if (nbrOfWorkers() > 1) {
  #  mysapply <- future_sapply
  #}
  #else {
  library(pbapply)
  mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
  #}
  all.features<-all.features[all.features%in%annotations$new_feature]
  newmat.list <- mysapply(X = 1:length(x = all.features), FUN = function(x) {
    features.use <- annotations[annotations$new_feature == 
        all.features[x], ]$feature
    #print(dim(peak.matrix[features.use, ]))
    submat <- peak.matrix[features.use, ]
    if (keep.sparse) {
      return(as(object = as.matrix(x = submat), Class = "dgCMatrix"))
    }else {
      return(as.matrix(x = submat))
    }
  }, simplify = FALSE)
  for ( i in 1:length(newmat.list)){
    if(dim(newmat.list[[i]])[2]==1){
      newmat.list[[i]]<-t(newmat.list[[i]])
    }
  }
  newmat = do.call(what = rbind, args = newmat.list)
  #newmat <- t(x = newmat)
  rownames(x = newmat) <- all.features
  return(as(object = newmat, Class = "dgCMatrix"))
}
