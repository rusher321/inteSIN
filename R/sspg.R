#' Sample-specific perturbation of gene interactions
#'
#' @param dat , matrix row is sample id col is feature
#' @param normalid ,vector including normal sample id or normal sample index
#' @param refnet, background net : a matrix include node1 node2
#'
#' @return list including delta rank matrix and edge perturbation matrix
#' @export
#'
#' @importFrom utils combn
#' @references
#' Chen Y, Gu Y, Hu Z, Sun X. Sample-specific perturbation of gene interactions identifies breast cancer subtypes. Brief Bioinform. 2021 Jul 20;22(4):bbaa268. doi: 10.1093/bib/bbaa268. PMID: 33126248; PMCID: PMC8293822.
#' @examples
#' num_samples <- 5
#' num_genes <- 6
#' gene_matrix <- matrix(runif(num_samples * num_genes), ncol = num_samples)
#' x <- as.data.frame(gene_matrix )
#' rownames(x) <- paste0("gene_", 1:num_genes)
#' colnames(x) <- paste0("sample_", 1:num_samples)
#' sspg_res <- sspg(t(x))
sspg <- function(dat, normalid = NULL, refnet = NULL){

  cat("dat: row is sample and column is feature\n")
  # transform the data to rank
  dat_rank <- apply(dat, 1, rank)
  sample_id <- colnames(dat_rank)
  gene_id <- rownames(dat_rank)
  gene_pair <- combn(gene_id, 2)
  if(is.null(refnet)){
    cat("if refnet is null, then use all pair\n")
    delta_rank <- sapply(1:ncol(dat_rank),
           function(x){ delta <- c();
           for(i in 1:ncol(gene_pair)){
             delta[i] <- dat_rank[gene_pair[1,i],x]-dat_rank[gene_pair[2,i], x]
             };
           delta})
    colnames(delta_rank) <- sample_id
    #delta_rank <- cbind(t(gene_pair), delta_rank)
    refnet <- t(gene_pair)
  }else{
    delta_rank <- sapply(1:ncol(dat_rank),
                         function(x){ delta <- c();
                         for(i in 1:nrow(refnet)){
                           delta[i] <- dat_rank[refnet[i,1],x]-dat_rank[refnet[i,2], x]
                         };
                         delta})
    colnames(delta_rank) <- sample_id
    #delta_rank <- cbind(refnet, delta_rank)
  }

  if(is.null(normalid)){
    batch_rank <- apply(delta_rank, 1, mean)
  }else{
    batch_rank <- apply(delta_rank[,normalid], 1, mean)
  }
  epm <- apply(delta_rank, 2, function(x){x-batch_rank})
  out_delta <- data.frame(cbind(refnet, delta_rank))
  out_epm <- data.frame(cbind(refnet, epm))

  out_delta[, 3:ncol(out_delta)] <- apply(out_delta[,3:ncol(out_delta)] , 2, as.numeric)
  out_epm[,3:ncol(out_epm)] <- apply( out_epm[,3:ncol(out_epm)] , 2, as.numeric)

  colnames(out_delta)[1:2] <- colnames(out_epm)[1:2] <- c("feature1", "feature2")
  outlist <- list(out_delta, out_epm)
  names(outlist) <- c("delta_rank", "epm")
  return(outlist)

}
