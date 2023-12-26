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
  print("Sparcc need use the read cound, specific for microbiome data")
  cor_matrix <- SpiecEasi::sparcc(x)
  return(cor_matrix)
}

.miF <- function(x){
  cat("mutual information matrix for nolinear correlation\n")
  mu_matrix <- minet::build.mim(dataset = x, estimator = "spearman")
  return(mu_matrix)
}

.clrF <- function(x){
  cat("Context Likelihood matrix for nolinear correlation\n")
  mu_matrix <- minet::build.mim(dataset = x, estimator = "spearman")
  return(minet::clr(mu_matrix))

}
