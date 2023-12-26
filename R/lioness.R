#' Linear Interpolation to Obtain Network Estimates for Single Samples
#'
#' @param matrix, row is gene/species, col is smaple/cell
#' @param n_method, function correlation method for feature-feature
#' @param transF, list output or matrix
#'
#' @importFrom reshape melt
#'
#' @references
#' Kuijjer ML, Tung MG, Yuan G, Quackenbush J, Glass K. Estimating Sample-Specific Regulatory Networks. iScience. 2019 Apr 26;14:226-240. doi: 10.1016/j.isci.2019.03.021 Add to Citavi project by DOI. Epub 2019 Mar 28. PMID: 30981959 Add to Citavi project by Pubmed ID; PMCID: PMC6463816.
#' @return specific sample network
#' @export
#'
#' @examples
#' num_samples <- 5
#' num_genes <- 6
#' gene_matrix <- matrix(runif(num_samples * num_genes), ncol = num_samples)
#' x <- as.data.frame(gene_matrix )
#' rownames(x) <- paste0("gene_", 1:num_genes)
#' colnames(x) <- paste0("sample_", 1:num_samples)
#' lioness <- lioness(t(x))
lioness <- function(matrix, n_method = .pearsonF, transF = F){

  if(!is.function(n_method)){ stop("please use a function") }
  if(!is.matrix(matrix)) { print("please use a numeric matrix as input") }
  cat("row is sample and colmun is feature!\n")
  cat("n_method need to be one of .pearsonF, .spearmanF, .kendallF and .sparccF!\n")

  samples <- rownames(matrix)
  nsamp <- nrow(matrix)
  # this applies netFun and extracts the aggregate network
  aggNet <- n_method(matrix)
  # prepare the lioness output
  lionessout <- list()
  # run function f and the LIONESS equation
  for(i in 1:nsamp){
    ss <- n_method(matrix[-i, ]) # apply netFun on all samples minus one
    lionessout[[samples[i]]] <- nsamp*(aggNet - ss)+ss # apply LIONESS equation
  }
  if(transF){
    # use the format gene-gene-correlation-sample
    lionessout <- melt(lionessout)
    colnames(lionessout) <- c("Feature_1", "Feature_2", "Rho", "Sample")
    lionessout <- lionessout[lionessout$Feature_1 != lionessout$Feature_2, ]
    lionessout <- lionessout[order(lionessout$Rho), ]
    lionessout <- lionessout[seq(1,nrow(lionessout), 2), ]
  }
  return(lionessout)
}


