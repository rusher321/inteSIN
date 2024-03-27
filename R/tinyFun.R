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


.wilcox_all <- function (datamatrix, configdata, ratio = "zero")
{

  id <- intersect(rownames(datamatrix), rownames(configdata))
  datamatrix <- datamatrix[id, ]
  configdata <- configdata[id, ]
  out <- matrix(NA, nrow = ncol(datamatrix), ncol = 9)
  config <- as.factor(as.character(configdata[, 1]))
  nlevels = levels(droplevels(config))

  for (i in 1:ncol(datamatrix)) {
    tmp <- as.numeric(datamatrix[, i])
    g1 <- tmp[which(config == nlevels[1])]
    g2 <- tmp[which(config == nlevels[2])]
    wilcox_sign <- wilcox.test(g1, g2)$p.value
    data_tmp <- data.frame(tmp = tmp, config = config)
    effect <- as.numeric(rstatix::wilcox_effsize(tmp~config, data = data_tmp)[1,4])

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
    z_score <- statistic(coin::wilcox_test(tmp~config, data_tmp))
    out[i, 1:9] <- c(wilcox_sign, or, mean_abun, median_abun,
                     z_score, effect)
  }
  rownames(out) <- colnames(datamatrix)
  colnames(out) <- c("wilcox_p.value", paste0(rep(nlevels,
                                                  3), rep(c("_ratio", "_mean", "_median"),
                                                          each = 2)), "z_score", "effect_size")
  out <- as.data.frame(out)
  out$p.adjust <- p.adjust(out[,1], method = "BH")
  out$enrich <- ifelse(out$p.adjust < 0.05, ifelse(out$z_score >
                                                     0, nlevels[1], nlevels[2]), "none")
  return(out)

}


.Compare_meta <- function(dat, config, prep = 0.2, method = "all", C = NULL, ...){

  id <- intersect(rownames(dat), rownames(config))
  data_match <- dat[id, ]
  Y <- config[id, 1]
  data_match <- data_match[, apply(data_match, 2, function(x){(sum(x!=0)/length(x)) >= prep})]

  outlist <- list()
  ############ locom compare
  if(method == "all" | method == "locom"){
    res.locom <- LOCOM::locom(otu.table = data_match, Y = Y, C = C,
                              ref.otu = NULL, fdr.nominal = 0.2,
                              n.perm.max = 50000, prev.cut = prep, n.rej.stop = 100, ...)

    out.locom <- data.frame(effect = t(res.locom$effect.size)[,1],
                            pvalue = t(res.locom$p.otu)[,1],
                            qvalue = t(res.locom$q.otu)[,1])

    rownames(out.locom) <- colnames(res.locom$effect.size)
    outlist[["locom"]] <- out.locom
  }

  ########### ancom-bc method
  sam.ID <- rownames(config)
  if (is.null(C)){
    meta.data <- data.frame(group = Y, row.names = id, stringsAsFactors = FALSE)
    formula <- "group"
  } else {
    meta.data <- data.frame(group = Y, C, row.names = id, stringsAsFactors = FALSE)
    colnames(meta.data)[2:ncol(meta.data)] <- colnames(C)
    formula <- paste("group + ",  paste(colnames(C), collapse = " + "), sep = "")
  }

  otu.table.ancombc <- t(data_match)
  otus.keep <- rownames(otu.table.ancombc)
  tax.table <- matrix(sample(letters, 7*length(otus.keep), replace = TRUE),
                      nrow = length(otus.keep), ncol = 7)
  rownames(tax.table) <- rownames(otu.table.ancombc)
  colnames(tax.table) <- c("Kingdom", "Phylum", "Class", "Order", "Family",
                           "Genus", "Species")

  otu.phylo <- phyloseq::otu_table(otu.table.ancombc, taxa_are_rows = TRUE)
  meta.phylo <- phyloseq::sample_data(meta.data)
  tax.phylo <- phyloseq::tax_table(tax.table)
  phylo.obj <- phyloseq::phyloseq(otu.phylo, meta.phylo, tax.phylo)

  if(method == "all" | method == "ancombc"){
    res.ancombc <- ANCOMBC::ancombc(phyloseq = phylo.obj, formula = formula,
                                    p_adj_method = "BH",
                                    zero_cut = 1-prep,
                                    lib_cut = 0,
                                    group = "group", struc_zero = TRUE, neg_lb = FALSE,
                                    tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05,
                                    global = TRUE)

    otu.ancombc <- data.frame(effect = res.ancombc$res$W[,1],
                              pvalue = res.ancombc$res$p_val[,1],
                              qvalue = res.ancombc$res$q_val[,1])

    rownames(otu.ancombc) <- rownames(res.ancombc$res$W)
    outlist[["ancombc"]] <- otu.ancombc
  }
  ##################### wilcox-raw
  if(method == "all" | method == "wilcox"){
    wilcox_res <- wilcox_all(datamatrix = data_match, configdata = meta.data)
    outlist[["wilcox_raw"]] <- wilcox_res
  }
  #################### wilcox-clr-psudo
  if(method == "all" | method == "clr"){
    data_match_p <- data_match
    data_match_p[data_match_p == 0] <- min(data_match[data_match!=0])
    data_match_clr <- robCompositions::cenLR(data_match_p)$x.clr
    data_match_clr_f <- data_match_clr[, colnames(data_match)]
    wilcox_res_clr <- wilcox_all(datamatrix = data_match_clr_f, configdata = meta.data)
    outlist[["wilcox_clr"]] <- wilcox_res_clr
  }

  ################### wilcox-clr-zcomposition
  if(method == "all" |method == "zclr"){
    data_match_z <- zCompositions::cmultRepl(data_match)
    data_match__zclr <- cenLR(data_match_z)$x.clr
    data_match_zclr_f <- data_match__zclr[, colnames(data_match)]
    wilcox_res_zclr <- wilcox_all(datamatrix = data_match_zclr_f, configdata = meta.data)
    outlist[["wilcox_res_zclr"]] <- wilcox_res_zclr
  }

  return(outlist)

}

