#' Search Candidate DNB genes
#'
#' Cluster by correlation, caculate DNB scores of all subtrees, and get best DNB
#' score for each time point
#' @param dnb a DNB object
#' @param min_size minimum gene number in subtree. Default 2.
#' @param max_size maximum gene number in subtree. Default Inf.
#' @param included_genes genes to be included in subtree. Default NULL.
#' @param all should all \@code{included_genes} be in subtree. Default False.
#' @param with_ctrl if consider control group. Default T.
#' @param verbose output progress.
#' @return a DNB object
#' @export
#' @references
#' 1. Liu, X., Chang, X., Liu, R., Yu, X., Chen, L., & Aihara, K. (2017). Quantifying critical states of complex diseases using single-sample dynamic network biomarkers. PLoS computational biology, 13(7), e1005633. https://doi.org/10.1371/journal.pcbi.1005633
search_dnb <- function(sni_obj, sid, type = "raw", min_size = 2,
                              max_size = Inf, included_genes = NULL, all = F,
                              verbose = T, method = "SSN", count=F) {

  stopifnot(class(sni_obj) == "SNI")
  candidates_dnb <- list()
  for (id in sid) {
    if (verbose) cat("sample id", id, "\n")
    dnb_lite <- NULL
    if(type == "raw"){
      corM <- sni_obj$SNI_raw[[id]]
    }else if(type == "adj"){
      corM <- sni_obj$SNI_adj[[id]]
    }
    # Hierarchical clustering genes ..
    #dend <- as.dendrogram(hclust(as.dist(1 - dnb$correlation[[tp]])))
    if(method == "SSN"){
      stopifnot("error: ssn method need use the raw type!"= type == "raw")
      dend <- as.dendrogram(hclust(as.dist(2-abs(corM)))) # ref1
    }

    modules <- dendextend::partition_leaves(dend)
    member_nums <- sapply(modules, length)
    modules <- modules[member_nums >= min_size & member_nums <= max_size]
    if (!is.null(included_genes)) {
      if (all == T) {
        modules <- modules[sapply(modules, function(x) {all(included_genes %in% x)})]
      }
      else {
        modules <- modules[sapply(modules, function(x) {any(included_genes %in% x)})]
      }
    }
    # compute the single Expression Deviation for every feature
    if(!is.null(sni_obj$data_ref)){
      feature_mean <- apply(sni_obj$data_ref, 2, dplR::tbrm)
      feature_s <- sni_obj$data[id, ]
      if(any(names(feature_mean) == names(feature_s))){
        if(count){
          sED <- abs(feature_s - feature_mean)/10^7
        }else{
          sED <- abs(feature_s - feature_mean)
        }
      }else{
        stop("error: data and data_ref need have same feature id!")
      }
    }else{
      feature_mean <- apply(sni_obj$data, 2, dplR::tbrm)
      feature_s <- sni_obj$data[id, ]
      if(count){
        sED <- abs(feature_s - feature_mean)/10^7
      }else{
        sED <- abs(feature_s - feature_mean)
      }
    }
    # compute the DNB score
    module_dnbs <- lapply(modules, function(m) .get_DNB_attr(m, corM, sED))
    candidates_dnb[[id]] <- module_dnbs[[which.max(sapply(module_dnbs, "[[", "score"))]]
  }
  return(candidates_dnb)
}



.get_DNB_attr <- function(genes_in, corM, sED) {

  correlation <- corM
  cv <- sED
  genes_out <- setdiff(names(cv), genes_in)
  cor_in <- mean(as.dist(abs(correlation[genes_in, genes_in])))
  cor_out <- mean(abs(correlation[genes_in, genes_out]))
  cv_in <- mean(cv[genes_in])
  # score
  score <- cv_in * cor_in / cor_out
  out <- list(  sED = sED,
                genes = genes_in,
                cv_in = cv_in,
                cor_in = cor_in,
                cor_out = cor_out,
                corM = corM,
                score = score)
  class(out) <- "dnb"
  return(out)

}

.get_final <- function(dnb_res) {
  id <- names(dnb_res)
  max_id <- which.max(sapply(id, function(x){dnb_res[[x]]$score}))
  dnb_genes <- dnb_res[[max_id]]$genes
  dnb_final <- lapply(id, function(x) {
    .get_DNB_attr(dnb_genes, dnb_res[[x]]$corM, dnb_res[[x]]$sED)
  })
  names(dnb_final) <- id
  return(dnb_final)
}

#' Plot DNB
#'
#' Plot DNB attibutes from data frame
#' @param dnb_df a data frame of DNB result.
#' @param ... additional arguments to ggplot2::ggsave
#' @export
plot_DNB <- function(dnb_df,  ...) {
  # score
  p_score <- plot_DNB_attr(dnb_df, "score")

  # cor_in
  p_cor_in <- plot_DNB_attr(dnb_df, "cor_in")

  # cor_out
  p_cor_out <- plot_DNB_attr(dnb_df, "cor_out")

  # cv
  p_cv <- plot_DNB_attr(dnb_df, "cv_in")

  p_merge <- cowplot::plot_grid(p_score, p_cor_in, p_cor_out, p_cv)

  p_merge
}

plot_DNB_attr <- function(dnb_df, attr_name) {
  id <- names(dnb_df)
  value <- sapply(id, function(x){dnb_df[[x]][[attr_name]]})
  dnb_df <- data.frame(time = paste0("T_", 1:length(id)), value)
  colnames(dnb_df)[2] <- attr_name
  dnb_df$time <- as.factor(dnb_df$time)
  ggplot2::ggplot(dnb_df, ggplot2::aes_string("time", attr_name, group =1)) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::geom_point(shape = 21, size = 3) +
    ggplot2::theme_classic()
}




