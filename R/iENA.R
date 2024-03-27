#' single sample pearson correlation in iENA
#'
#' @param dat, matrix, row is gene/species, col is smaple/cell
#' @param PID, vector, sample id for specific or all
#' @param transF, F = net output; T = matrix output
#'
#' @return list or matrix
#' @export
#' @importFrom stats var
#' @references
#' Yu X, Zhang J, Sun S, Zhou X, Zeng T, Chen L. Individual-specific edge-network analysis for disease prediction. Nucleic Acids Res. 2017 Nov 16;45(20):e170. doi: 10.1093/nar/gkx787 Add to Citavi project by DOI. PMID: 28981699 Add to Citavi project by Pubmed ID; PMCID: PMC5714249.
#' @examples
#' num_samples <- 5
#' num_genes <- 6
#' gene_matrix <- matrix(runif(num_samples * num_genes), ncol = num_samples)
#' x <- as.data.frame(gene_matrix )
#' rownames(x) <- paste0("gene_", 1:num_genes)
#' colnames(x) <- paste0("sample_", 1:num_samples)
#' iENA_sPCC_res <- iENA_sPCC(t(x))
iENA_sPCC <- function(dat, PID = "all", transF = F, alpha = 0.05, ref = F, refmean = NULL, refvar=NULL){

  if(PID != "all"){
    # remove the id not in the dat
    PID <- rownames(dat)[which(rownames(dat) %in% PID)]
  }else{
    PID <- rownames(dat)
  }
  if(ref){
    ref_mean <- refmean
    ref_var <- refvar
  }else{
    ref_mean <- apply(dat, 2, mean)
    ref_var <- apply(dat, 2, var)
  }
  n_sample <- nrow(dat)
  n_feature <- ncol(dat)
  V1 <- matrix(ref_var, n_feature, n_feature)
  V2 <- matrix(ref_var, n_feature, n_feature, byrow = T)
  G1 <- matrix(ref_mean, n_feature, n_feature)
  G2 <- matrix(ref_mean, n_feature, n_feature, byrow = T)

  feature_id <- colnames(dat)
  sPCC <- list()
  sPCC_z <- list()

  for(i in 1:length(PID)){
    Xi <- as.numeric(dat[i, ])
    X1 <- matrix(Xi, n_feature, n_feature)
    X2 <- matrix(Xi, n_feature, n_feature, byrow = T)
    sPCCi <- (X1-G1)*(X2-G2)/sqrt(V1*V2) # here use the equitation for iENA paper
    rownames(sPCCi) <- colnames(sPCCi) <- feature_id
    sPCC[[PID[i]]] <- sPCCi
    # here to generate the z-score matrix
    z_ss <- (sPCCi-mean(sPCCi[lower.tri(sPCCi)]))/sd(sPCCi[lower.tri(sPCCi)])
    sPCC_z[[PID[i]]] <- z_ss
  }

  if(transF){
    # use the format gene-gene-correlation-sample
    sPCC <- lapply(sPCC, function(x){tmp <- as.matrix(x); diag(tmp) <- NA; tmp[lower.tri(tmp)] <- NA;tmp})
    sPCC <- melt(sPCC)
    colnames(sPCC) <- c("Feature_1", "Feature_2", "Rho", "Sample")
    sPCC <- sPCC[!is.na(sPCC$Rho), ]
  }else{
    cutoff <- -qnorm(alpha, mean = 0, sd = 1)
    sPCC <- lapply(sPCC_z, function(x){x[abs(x) < cutoff] <- 0; x <- sign(x); as(x, "dgCMatrix")})
    names(sPCC) <- PID
  }
  return(sPCC)
}


# iENA_sSpearman <- function(dat, PID = "all", transF = F){
#
#   if(PID != "all"){
#     # remove the id not in the dat
#     PID <- rownames(dat)[which(rownames(dat) %in% PID)]
#   }else{
#     PID <- rownames(dat)
#   }
#
#   ref_mean <- apply(dat, 2, mean)
#   ref_var <- apply(dat, 2, var)
#   n_sample <- nrow(dat)
#   n_feature <- ncol(dat)
#   V1 <- matrix(ref_var, n, n)
#   V2 <- matrix(ref_var, n, n, byrow = T)
#   G1 <- matrix(ref_mean, n, n)
#   G2 <- matrix(ref_mean, n, n, byrow = T)
#
#   feature_id <- colnames(dat)
#   sPCC <- list()
#
#   for(i in 1:length(PID)){
#     Xi <- as.numeric(dat[i, ])
#     X1 <- matrix(Xi, n, n)
#     X2 <- matrix(Xi, n, n, byrow = T)
#     sPCCi <- (X1-G1)*(X2-G2)/sqrt(V1*V2) # here use the equotation for iENA paper
#     rownames(sPCCi) <- colnames(sPCCi) <- feature_id
#     sPCC[[PID[i]]] <- sPCCi
#   }
#
#   if(transF){
#     # use the format gene-gene-correlation-sample
#     sPCC <- melt(sPCC)
#     colnames(sPCC) <- c("Feature_1", "Feature_2", "Rho", "Sample")
#     sPCC <- sPCC[sPCC$Feature_1 != sPCC$Feature_2, ]
#     sPCC <- sPCC[order(sPCC$Rho), ]
#     sPCC <- sPCC[seq(1,nrow(sPCC), 2), ]
#   }
#   return(sPCC)
# }
