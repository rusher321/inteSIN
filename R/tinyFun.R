.splitMetaphlan <- function (data, prefix){

  data <- data[grep("k__Bacteria", rownames(data)), ]
  index <- sapply(strsplit(rownames(data), perl = T, split = "\\|"),
                  function(x){rev(x)[1]})

  datalist <- list()
  classify <- c("kingdom","phylum", "class",
                "order", "family", "genus", "species")
  classify2 <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")

  for (i in 1:7) {
    index2 <- grep(classify2[i], index)
    outdat <- data[index2, ]
    rname <- sapply(strsplit(rownames(outdat), perl = T,
                             split = "\\|"), function(x) {
                               y <- rev(x)[1]
                               return(y)
                             })
    rownames(outdat) <- rname
    datalist[[i]] <- outdat
  }
  names(datalist) <- paste0(prefix, "_", classify[1:7])
  return(datalist)
}

.topTax <- function(metadata, K = 20, rmUnclass = F, sample_order =NULL ,
                   tax_order = NULL, tax_colour=NULL){

  # generate the top tax profile
  if(rmUnclass){
    # rm the unclass & renorm
    unclassindex <- which(rownames(metadata) %in% "unclassed")
    metadata <- metadata[-unclassindex, ]
    metadata <- apply(metadata, 2, function(x){y <- x/sum(x); return(y)})
  }

  # order the tax
  if(!is.null(tax_order)){
    data_sub <- metadata[tax_order, ]
    lessdata <- metadata[-which(rownames(metadata) %in% tax_order), ]
    otherdata <- t(data.frame(apply(lessdata, 2, sum)))
    rownames(otherdata) <- "Others"
    qdat <- rbind(data_sub, otherdata)
    naOrder <- c("Others", tax_order)
  }else{
    top <- names(head(sort(apply(metadata , 1, mean), decreasing = T), K))
    topdata <- metadata[top, ]
    lessdata <- metadata[-which(rownames(metadata) %in% top), ]
    otherdata <- t(data.frame(apply(lessdata, 2, sum)))
    rownames(otherdata) <- "Others"
    qdat <- rbind(topdata, otherdata)
    naOrder <- rownames(qdat)
  }


  # order the sample
  if(!is.null(sample_order)){
    idOrder = sample_order
  }else{
    idOrder <- colnames(qdat)[order(qdat[1,], decreasing = T)]
  }
  qdat <- as.data.frame(qdat, check.names=F)
  qdat$sample <- rownames(qdat)
  qdat2 <- melt(qdat)
  colnames(qdat2) <- c("Tax", "Sample", "value")
  qdat2$Tax <- factor(qdat2$Tax, levels = rev(naOrder))
  qdat2$Sample <- factor(qdat2$Sample, levels = idOrder)
  # ggplot
  p = ggplot(qdat2, aes(Sample, value, fill=Tax))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 1)+
    theme_classic()+
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank())+
    xlab("")+ylab("Relative Abundance")+
    scale_y_continuous(expand = c(0,0),breaks = c(0,.2,.4,.6,.8,1))
  if(!is.null(tax_colour)){
    p <- p+scale_fill_manual(values = tax_colour)
  }

  return(p)
}


.net_stat <- function(net_adj){

  cor_matrix <- as.matrix(net_adj)
  feature_name <- rownames(cor_matrix)
  # compute the degree between closeness
  edge_pos <-   sum(cor_matrix[lower.tri(cor_matrix)] == 1)
  edge_neg <-   sum(cor_matrix[lower.tri(cor_matrix)] == -1)

  net_adj <- adj2igraph(net_adj)
  degree_n <-    data.frame(degree = igraph::degree(net_adj), row.names = feature_name)
  E(net_adj)$weight <- abs(E(net_adj)$weight) # translate weight to postive
  betwenness_n <- data.frame(betweeness = igraph::betweenness(net_adj), row.names= feature_name)
  closeness_n <- data.frame(closeness = igraph::closeness(net_adj), row.names = feature_name)
  out <- list(edge_pos, edge_neg, degree_n, betwenness_n,
              closeness_n)
  names(out) <- c("edge_pos", "edge_net", "degree", "betwenness", "closeness")

  return(out)
}

.net_stat_com <- function(netlist){

  all_res <- lapply(netlist, .net_stat)
  pID <- names(netlist)
  edge_num <- data.frame(pos_edge = sapply(all_res, function(x){x[[1]]}),
                         neg_edge = sapply(all_res, function(x){x[[2]]}), row.names = pID)

  degree_stat <- do.call("cbind", lapply(all_res, function(x){x[[3]]}))
  betw_stat <- do.call("cbind", lapply(all_res, function(x){x[[4]]}))
  clos_stat <- do.call("cbind", lapply(all_res, function(x){x[[5]]}))

  colnames(degree_stat) <- colnames(betw_stat) <- colnames(clos_stat) <- pID

  out <- list(edge_num, degree_stat, betw_stat, clos_stat)
  names(out) <- c("edge_number", "degree", "betweeness", "closeness")
  return(out)

}

.mywilcox_2t <- function (datamatrix, configdata, group, ratio = "zero",
                            ...)
{
  config <- .matchpairID(configdata, ...)
  data <- datamatrix[rownames(config), , drop = F]
  out <- matrix(NA, nrow = ncol(data), ncol = 9)
  config <- as.factor(as.character(config[, group]))
  nlevels = levels(droplevels(config))
  for (i in 1:ncol(data)) {
    tmp <- as.numeric(data[, i])
    if(sum(tmp) == 0){
      out[i, 1:9] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
    }else{
      g1 <- tmp[config == nlevels[1]]
      g2 <- tmp[config == nlevels[2]]
      if(sum(g1-g2)==0){
        out[i, 1:9] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
      }else{
        wilcox_sign <- pvalue(wilcoxsign_test(g1 ~ g2))
        effect <- wilcoxonPairedR(x <- tmp, g <- config)
        if (ratio == "zero") {
          or <- tapply(tmp, config, function(x) {
            sum(x != 0, na.rm = T)/length(x)
          })
        }
        else {
          or <- tapply(tmp, config, function(x) {
            sum(!is.na(x))/length(x)
          })
        }
        mean_abun <- tapply(tmp, config, mean, na.rm = T)
        median_abun <- tapply(tmp, config, median, na.rm = T)
        z_score <- statistic(wilcoxsign_test(g1 ~ g2))
        out[i, 1:9] <- c(wilcox_sign, or, mean_abun, median_abun,
                         z_score, effect)
      }
    }
  }
  rownames(out) <- colnames(data)
  colnames(out) <- c("sign_p.value", paste0(rep(nlevels,
                                                3), rep(c("_ratio", "_mean", "_median"),
                                                        each = 2)), "z_score", "effect_size")
  out <- as.data.frame(out)
  out$p.adjust <- p.adjust(out$sign_p.value, method = "BH")
  out$enrich <- ifelse(out$p.adjust < 0.05, ifelse(out$z_score >
                                                     0, nlevels[1], nlevels[2]), "none")
  return(out)
}


.matchpairID <- function(configdat, ID, Time, num = 2){

  names(table(configdat[,ID]))[table(configdat[, ID]) == num] -> matchname
  configdat[!is.na(match(configdat[, ID], matchname)),] -> outconfig
  # sort by ID
  outconfig[order(outconfig[,ID]), ] -> matchdat
  matchdat[order(matchdat[,Time]), ] -> matchdat

  return(matchdat)

}

