#' Linear Interpolation to Obtain Network Estimates for Single Samples
#'
#' @param matrix inputdata: row is gene/species, col is sample/cell
#' @param n_method function :correlation method for feature-feature including pearson
#' ,.spearman, .kendall, .sparcc, .MI(mutual information) and .CLR(Context Likelihood)
#' @param alpha significant cutoff, default =0.05
#' @param dat_ref matrix of the reference dataset
#' @param verbose logical
#' @return SNI object:
#' \item{data}{your input table}
#' \item{data_ref}{your reference data}
#' \item{SNI_raw}{the raw value of pertubation matrix for each sample}
#' \item{SNI_z}{the z value of pertubation matrix for each sample}
#' \item{SNI_adj}{the adjcent matrix or network  for each sample}
#' \item{SNI_cen}{the degree/betweeness/closeness for each sample}
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats cor pnorm sd
#' @importFrom SpiecEasi adj2igraph
#' @import igraph
#' @references Kuijjer ML, Tung MG, Yuan G, Quackenbush J, Glass K. Estimating Sample-Specific Regulatory Networks. iScience. 2019 Apr 26;14:226-240. doi: 10.1016/j.isci.2019.03.021.
#' @family SNI
#' @export
#'
#' @examples
#' num_samples <- 100
#' num_genes <- 60
#' gene_matrix <- matrix(runif(num_samples * num_genes), ncol = num_samples)
#' x <- as.data.frame(gene_matrix )
#' rownames(x) <- paste0("gene_", 1:num_genes)
#' colnames(x) <- paste0("sample_", 1:num_samples)
#' lioness_res <- lioness(t(x))
lioness <- function(matrix, n_method = .pearsonF,  alpha = 0.05,
                dat_ref = NULL, verbose = F){
  if(verbose){
    message("row is sample and colmun is feature!\n")
    message("n_method can be one of .pearsonF, .spearmanF, .kendallF, sparccF,
            .mlF, .clrF or you can design one.\n")
  }
  samples <- rownames(matrix)
  nsamp <- nrow(matrix)
  # generate the aggregate network
  if(is.null(dat_ref)){
    aggNet <- n_method(matrix)
  }else{
    aggNet <- n_method(dat_ref)
  }
  cutoff <- -qnorm(alpha, mean = 0, sd = 1)
  # prepare the lioness output
  lionessout_z <- list()
  lionessout_adjNet <- list()
  lionessout_netCent <- list()
  lionessout_raw <- list()
  # run function f and the SS equation
  pb = txtProgressBar(style = 3)
  for(i in 1:nsamp){
    #Sys.sleep(0.5);
    setTxtProgressBar(pb, i/nsamp)
    if(is.null(dat_ref)){
      ss <- n_method(matrix[-i, ]) # apply netFun on all samples minus one
      sin_tmp <- nsamp*(aggNet - ss)+ss # apply LIONESS equation
    }else{
      ss <- n_method(dat_ref)
      aggNet_tmp <- n_method(rbind(matrix[i,], dat_ref))
      sin_tmp <- nsamp*(aggNet - ss)+ss
    }

    # generate the z-score matrix, need to discuss
    z_ss <- (sin_tmp - mean(sin_tmp[lower.tri(sin_tmp)]))/sd(sin_tmp[lower.tri(sin_tmp)])
    lionessout_z[[samples[i]]] <- z_ss
    lionessout_raw[[samples[i]]] <- sin_tmp

    # generate the adjacent matrix
    adj_ss <- z_ss
    adj_ss[abs(adj_ss) < cutoff] <- 0
    adj_ss <- as(sign(adj_ss), "dgCMatrix")
    lionessout_adjNet[[samples[i]]] <- adj_ss
    # generate the network centrality
    lionessout_netCent[[samples[i]]] <- .net_stat(adj_ss)
  }
  Sys.sleep(1);close(pb)

  # combine network centrality
  names(lionessout_raw) <- names(lionessout_z) <- names(lionessout_adjNet) <- names(lionessout_netCent) <- samples

  degree <- do.call("cbind", lapply(lionessout_netCent, function(x){x$degree}) )
  betweenness <- do.call("cbind", lapply(lionessout_netCent, function(x){x$betwenness}))
  closeness <- do.call("cbind", lapply(lionessout_netCent, function(x){x$closeness}))
  colnames(degree) <- colnames(betweenness) <- colnames(closeness) <- samples

  lionessout_netCent_comb <- list(degree, betweenness, closeness)
  names(lionessout_netCent_comb) <- c("degree", "betweenness", "closeness")
  out <- list(data = matrix, data_ref = dat_ref, SNI_raw = lionessout_raw, SNI_z = lionessout_z,
              SNI_adj = lionessout_adjNet,
              SNI_cen = lionessout_netCent_comb)

  class(out) <- "SNI"
  return(out)
}


