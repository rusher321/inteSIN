#' Sample-specific weighted correlation network
#'
#' @param dat, row is gene/species, col is smaple/cell
#' @param s_method, function correlation method for sample-sample
#' @param n_method, function correlation method for feature-feature
#' @param balance, numeric,
#' @param PID, sample id for specific or all
#' @param alpha, significant cutoff
#' @param out_z, logistic, F = output is zscore matrix, T = adj matrix
#' @return outlist, a list including specific network
#' @export
#' @references
#' Chen HH, Hsueh CW, Lee CH, Hao TY, Tu TY, Chang LY, Lee JC, Lin CY. SWEET: a single-sample network inference method for deciphering individual features in disease. Brief Bioinform. 2023 Mar 19;24(2):bbad032. doi: 10.1093/bib/bbad032 Add to Citavi project by DOI. PMID: 36719112 Add to Citavi project by Pubmed ID; PMCID: PMC10025435.
#' @examples
#' num_samples <- 5
#' num_genes <- 6
#' gene_matrix <- matrix(runif(num_samples * num_genes), ncol = num_samples)
#' x <- as.data.frame(gene_matrix )
#' rownames(x) <- paste0("gene_", 1:num_genes)
#' colnames(x) <- paste0("sample_", 1:num_samples)
#' sweet_res <- Sweet(t(x))

Sweet <- function(dat, s_method = .pearsonF, n_method = .pearsonF, balance = 0.1, PID = "all", alpha = 0.05, out_z = F){

  # calculate the weight
  cat("row is sample and colmun is feature!\n")
  cat("s_method and n_method need to be one of .pearsonF, .spearmanF, .kendallF, .miF, .clrF and .sparccF!\n")

  corM <- s_method(t(dat))
  pat_N <- nrow(dat)
  value <- (apply(corM, 1, sum)-1)/(pat_N-1)
  rmax <- max(value)
  rmin <- min(value)
  diff <- rmax - rmin +0.01
  value <- (value - rmin+0.01)/diff
  value <- value * balance * pat_N
  value_d <- data.frame(weight = value, row.names = rownames(dat))

  # calculate the coefficient of raw edge for each sample
  # remove feature with sd =0
  dat <- dat[, apply(dat, 2, sd)!=0]
  corN <- n_method(dat)

  if(PID != "all"){
    # remove the id not in the dat
    PID <- rownames(dat)[which(rownames(dat) %in% PID)]
  }else{
    PID <- rownames(dat)
  }
  netlist <- lapply(PID, function(x){tmp_dat <- rbind(dat, dat[x,]);
  p_cor <- n_method(tmp_dat);
  p_cor <- value_d[x,1]*(p_cor-corN)+corN;
  diag(p_cor) <- 0; p_cor})

  # calculate the z-score of edge for each sample
  all_value <- lapply(netlist, function(x){x[lower.tri(x)]})
  all_value <- unlist(all_value)
  mean_n <- mean(all_value)
  sd_n <- sd(all_value)
  outlist_z <- lapply(netlist, function(x){p_adj <- (x-mean_n)/sd_n;
  diag(p_adj) <- 0; p_adj})

  # calculate the adjacent matrix
  cutoff <- -qnorm(alpha, mean = 0, sd = 1)
  outlist <- lapply(outlist_z, function(x){x[abs(x) < cutoff] <- 0; x <- sign(x); as(x, "dgCMatrix")})
  names(outlist) <- names(outlist_z) <- PID
  if(out_z){
    return(outlist_z)
  }else{
    return(outlist)
  }

}
