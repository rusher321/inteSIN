.pearsonF <- function(x){
  cor_matrix <- cor(x, method = "p")
  return(cor_matrix)
}

.spearmanF <- function(x){
  cor_matrix <- cor(x, method = "spearman")
  return(cor_matrix)
}

.kendelF <- function(x){
  cor_matrix <- cor(x,method = "kendel")
  return(cor_matrix)
}

.sparccF <- function(x){
  #message("Sparcc need use the read cound, specific for microbiome data")
  cor_matrix <- SpiecEasi::sparcc(x)
  out <- cor_matrix$Cor
  colnames(out) <- rownames(out) <- colnames(x)
  return(out)
}

.miF <- function(x){
  #cat("mutual information matrix for nolinear correlation\n")
  mu_matrix <- minet::build.mim(dataset = x, estimator = "spearman")
  return(mu_matrix)
}

.clrF <- function(x){
  #cat("Context Likelihood matrix for nolinear correlation\n")
  mu_matrix <- minet::build.mim(dataset = x, estimator = "spearman")
  return(minet::clr(mu_matrix))

}
