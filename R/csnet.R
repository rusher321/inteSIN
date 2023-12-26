#' Construct Cell specific network
#'
#' @param data row is gene/species, col is smaple/cell
#' @param c numeric, the id of cell you intersested, null means to compute all the sample
#' @param alpha numeric, Significant level ;eg. 0.001, 0.01, 0.05 ...
#' @param boxsize numeric, default 0.1, to define the upper and lower
#' @param weighted False or True
#'
#' @return list, the first is single cell network, the second is the degree matrix
#'
#' @importFrom igraph degree
#' @importFrom Matrix sparseMatrix
#' @importFrom stats qnorm
#'
#' @export
#'
#' @references
#' Dai H, Li L, Zeng T, Chen L. Cell-specific network constructed by single-cell RNA sequencing data. Nucleic Acids Res. 2019 Jun 20;47(11):e62. doi: 10.1093/nar/gkz172
#' @examples
#' num_samples <- 5
#' num_genes <- 6
#' gene_matrix <- matrix(runif(num_samples * num_genes), ncol = num_samples)
#' x <- as.data.frame(gene_matrix )
#' rownames(x) <- paste0("gene_", 1:num_genes)
#' colnames(x) <- paste0("sample_", 1:num_samples)
#' csns <- csnet(x, c = NULL, alpha = 0.01, boxsize = 0.1, weighted = FALSE)

csnet <- function(data, c = NULL, alpha = 0.01, boxsize = 0.1, weighted = FALSE) {
  # Construction of cell-specific networks
  sid <- colnames(data)
  gid <- rownames(data)
  # If c is not provided, construct the CSNs for all cells
  if (is.null(c)) {
    c <- seq_len(ncol(data))
  }

  # Default values
  if (is.null(weighted)) {
    weighted <- FALSE
  }
  if (is.null(boxsize)) {
    boxsize <- 0.1
  }
  if (is.null(alpha)) {
    alpha <- 0.01
  }

  n1 <- nrow(data)
  n2 <- ncol(data)

  # Define the neighborhood for each plot
  upper <- lower <- matrix(0, n1, n2)
  for(i in seq_len(n1)){
    s1 <- sort(data[i, ])
    s2 <- order(data[i, ])
    n3 <- n2 - sum(sign(s1))
    h <- round(boxsize/2 * sum(sign(s1)))
    k <- 1
    while(k <= n2){
      s = 0
      while(k + s + 1 <= n2 && s1[k + s + 1] == s1[k]){
         s = s+1
      }
      if(s >= h){
        upper[i, s2[k:(k + s)]] <- data[i, s2[k]]
        lower[i, s2[k:(k + s)]] <- data[i, s2[k]]
      }else{
        upper[i, s2[k:(k + s)]] <- data[i, s2[min(n2, k + s + h)]]
        lower[i, s2[k:(k + s)]] <- data[i, s2[max(n3 * (n3 > h) + 1, k - h)]]
      }
      k = k+s+1
    }
  }
  #print(upper)
  #print(lower)
  # Construction of cell-specific network
  # Initialize csn list
  csn <- vector("list", length = length(c))
  B <- matrix(0, n1, n2)
  p <- -qnorm(alpha, mean = 0, sd = 1)
  #csndm <- matrix(NA, nrow = n1, ncol = n2)

  for(k in c){
    for(j in seq_len(n2)){
      B[, j] <- data[, j] <= upper[, k] & data[, j] >= lower[, k]
    }
    a <- rowSums(B)
    d <- (crossprod(t(B), t(B)) * n2 - outer(a, a)) / sqrt(outer(a, a) * ((n2 - a) %*% t(n2 - a)) / (n2 - 1) + .Machine$double.eps)
    diag(d) <- 0
    #print(d)
    if (weighted) {
      d <- d*(d>0)
      csn[[k]] <- sparseMatrix(i = row(d)[d > p], j = col(d)[d > p], x = rep(1, sum(d > p)), dims = dim(d), dimnames = list(gid, gid))
    } else {
      csn[[k]] <- sparseMatrix(i = row(d)[d > p], j = col(d)[d > p], x = rep(1, sum(d > p)), dims = dim(d), dimnames = list(gid, gid))
    }
    #csndm[,k] <- degree(.adj2igraph(csn[[k]]))
    cat(paste("Sample/Cell ", k, " is completed\n"))

  }
  #names(out) <- c("CSN", "CSNDM")
  names(csn) <- sid
  return(csn)
}

.adj2igraph <- function (Adj, rmEmptyNodes = FALSE, diag = FALSE, edge.attr = list(),
          vertex.attr = list(name = 1:ncol(Adj)))
{
  g <- igraph::graph.adjacency(Adj, mode = "undirected",
                               weighted = TRUE, diag = diag)
  if (length(vertex.attr) > 0) {
    for (i in 1:length(vertex.attr)) {
      attr <- names(vertex.attr)[i]
      g <- igraph::set.vertex.attribute(g, attr, index = igraph::V(g),
                                        vertex.attr[[i]])
    }
  }
  if (length(edge.attr) > 0) {
    for (i in 1:length(edge.attr)) {
      attr <- names(edge.attr)[i]
      g <- igraph::set.edge.attribute(g, attr, index = igraph::E(g),
                                      edge.attr[[i]])
    }
  }
  if (rmEmptyNodes) {
    ind <- igraph::V(g)$name[which(igraph::degree(g) < 1)]
    g <- igraph::delete.vertices(g, ind)
  }
  g
}
