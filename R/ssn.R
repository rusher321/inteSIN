#' SSN:sample-specific networks
#'
#' @param matrix, row is gene/species, col is smaple/cell
#' @param n_method, function :correlation method for feature-feature
#' @param transF, F means list output; T means matrix
#'
#' @return specific sample network
#'
#' @importFrom reshape melt
#' @importFrom stats cor pnorm sd
#' @references Liu X, Wang Y, Ji H, Aihara K, Chen L. Personalized characterization of diseases using sample-specific networks. Nucleic Acids Res. 2016 Dec 15;44(22):e164. doi: 10.1093/nar/gkw772
#'
#' @export
#'
#' @examples
#' num_samples <- 5
#' num_genes <- 6
#' gene_matrix <- matrix(runif(num_samples * num_genes), ncol = num_samples)
#' x <- as.data.frame(gene_matrix )
#' rownames(x) <- paste0("gene_", 1:num_genes)
#' colnames(x) <- paste0("sample_", 1:num_samples)
#' ssn_res <- ssn(t(x))
ssn <- function(matrix, n_method = .pearsonF, transF = F, alpha = 0.05){

  if(!is.function(n_method)){ stop("please use a function") }
  if(!is.matrix(matrix)) { print("please use a numeric matrix as input") }
  cat("row is sample and colmun is feature!\n")
  cat("n_method need to be one of .pearsonF, .spearmanF, .kendallF and .sparccF!\n")

  samples <- rownames(matrix)
  nsamp <- nrow(matrix)
  # this applies netFun and extracts the aggregate network
  aggNet <- n_method(matrix)
  # prepare the ssn output
  ssout_z <- list()
  ssout_p <- list()
  # run function f and the SS equation
  for(i in 1:nsamp){
    ss <- n_method(rbind(matrix, matrix[i, ])) # apply n_method on all samples plus one
    delta_ss <- ss-aggNet # apply ssn equation
    r_ss <- (1-delta_ss*delta_ss)/(nsamp-2)
    z_ss <- delta_ss/r_ss
    p_ss <- 2*pnorm(-abs(z_ss), 0, 1) # z-score to P
    ssout_z[[samples[i]]] <- z_ss
    ssout_p[[samples[i]]] <- p_ss
  }

  #ssout <- list(ssout_z, ssout_p)
  #names(ssout) <- c("Z_score", "Pvalue")

  if(transF){
    # use the format gene-gene-correlation-sample
    ssout_z <- lapply(ssout_z, function(x){tmp <- as.matrix(x); diag(tmp) <- NA; tmp[lower.tri(tmp)] <- NA;tmp})
    ssout <- melt(ssout_z)
    ssout_p <- melt(ssout_p)
    colnames(ssout) <- c("Feature_1", "Feature_2", "Z_score", "Sample")
    ssout$Pvalue <- ssout_p[,3]
    ssout <- ssout[!is.na(ssout$Z_score), ] # remove the loop and repeat edga
  }else{
    cutoff <- -qnorm(alpha, mean = 0, sd = 1)
    ssout <- lapply(ssout_z, function(x){x[abs(x) < cutoff] <- 0; x <- sign(x); as(x, "dgCMatrix")})
    names(ssout) <- samples
  }
  return(ssout)
}
