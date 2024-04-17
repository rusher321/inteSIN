#' SSN:sample-specific networks
#'
#' @param matrix inputdata: row is gene/species, col is sample/cell
#' @param n_method function :correlation method for feature-feature including pearson
#' spearman, kendall, sparcc, mutual information and Context Likelihood
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
#' @references Liu X, Wang Y, Ji H, Aihara K, Chen L. Personalized characterization of diseases using sample-specific networks. Nucleic Acids Res. 2016 Dec 15;44(22):e164. doi: 10.1093/nar/gkw772
#'
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
#' ssn_res <- ssn(t(x))
ssn <- function(matrix, n_method = .pearsonF,  alpha = 0.05,
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
  # prepare the ssn output
  ssout_z <- list()
  cutoff <- -qnorm(alpha, mean = 0, sd = 1)
  ssout_adjNet <- list()
  ssout_netCent <- list()
  ssout_delta <- list()
  # run function f and the SS equation
  pb = txtProgressBar(style = 3)
  for(i in 1:nsamp){
    #Sys.sleep(0.5);
    setTxtProgressBar(pb, i/nsamp)
    if(is.null(dat_ref)){
      # apply n_method on all samples plus one
      ss <- n_method(rbind(matrix, matrix[i, ]))
    }else{
      ss <- n_method(rbind(dat_ref, matrix[i, ]))
    }
    # apply ssn equation
    delta_ss <- ss-aggNet
    ssout_delta[[samples[i]]] <- delta_ss
    r_ss <- (1-delta_ss*delta_ss)/(nsamp-2)
    z_ss <- delta_ss/r_ss
    ssout_z[[samples[i]]] <- z_ss
    # generate the adjacent matrix
    adj_ss <- z_ss
    adj_ss[abs(adj_ss) < cutoff] <- 0
    adj_ss <- as(sign(adj_ss), "dgCMatrix")
    ssout_adjNet[[samples[i]]] <- adj_ss
    # generate the network centrality
    ssout_netCent[[samples[i]]] <- .net_stat(adj_ss)
  }
  Sys.sleep(1);close(pb)

  # combine network centrality
  names(ssout_z) <- names(ssout_adjNet) <- names(ssout_netCent) <- samples

  degree <- do.call("cbind", lapply(ssout_netCent, function(x){x$degree}) )
  betweenness <- do.call("cbind", lapply(ssout_netCent, function(x){x$betwenness}))
  closeness <- do.call("cbind", lapply(ssout_netCent, function(x){x$closeness}))
  colnames(degree) <- colnames(betweenness) <- colnames(closeness) <- samples

  ssout_netCent_comb <- list(degree, betweenness, closeness)
  names(ssout_netCent_comb) <- c("degree", "betweenness", "closeness")
  out <- list(data = matrix, data_ref = dat_ref, SNI_raw = ssout_delta, SNI_z = ssout_z, SNI_adj = ssout_adjNet,
              SNI_cen = ssout_netCent_comb)

  class(out) <- "SNI"
  return(out)
}
