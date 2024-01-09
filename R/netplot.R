#' network plot for single sample network
#'
#' @param matrix sparse adjacent matrix or adjacent matrix including gene id or species name
#' @param layout layout of network: one of auto, igraph (from layout_tbl_graph_igraph), dendrogram, manual, linear, matrix, treemap, circlepack, partition or hive. the detail in help of ggraph.
#' @param node_size numeric or character
#' @param node_color color or one variable from metadata
#' @param edge_size  numeric
#' @param edge_color color or "type"
#' @param label_size  numeric or "weight_abs"
#' @param label_color color or one variable from metadata
#' @param metadata  metadata information of node
#'
#' @return ggplot class
#' @importFrom igraph graph_from_data_frame
#' @importFrom  ggraph ggraph
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' num_samples <- 50
#' num_genes <- 100
#' gene_matrix <- matrix(runif(num_samples * num_genes), ncol = num_samples)
#' x <- as.data.frame(gene_matrix )
#' rownames(x) <- paste0("gene_", 1:num_genes)
#' colnames(x) <- paste0("sample_", 1:num_samples)
#' Sweet_res <- Sweet(t(x))
#' netplot_sin(sweet_res$sample_1,layout = "circle", node_size = "degree", label_color = "black", edge_color = "type")
netplot_sin <- function(matrix, layout = "sphere", node_size = 3, node_color = "orange",
                        edge_size =1, edge_color= "gray", label_size = 3, label_color = "grey50",
                        metadata = NULL){

  ## this function based on ggraph

  links <- .translink(matrix) # get the link
  links$Weight_abs <- abs(links$Weight)
  degree <- apply(as.matrix(matrix), 1, function(x){sum(abs(x))})
  id_name <- rownames(as.matrix(matrix))
  nodes <- data.frame(id = id_name, degree = degree, id_name = id_name)
  nodes <- nodes[nodes$degree != 0, ]

  if(!is.null(metadata)){
    nodes <- rbind(nodes, metadata[nodes$id, ,drop = F])
  }
  net <- graph_from_data_frame(d=links, vertices=nodes, directed = F) # translate to igraph class

  raw_plot <- ggraph(graph = net, layout = layout)

  if(!is.numeric(edge_size)){
    raw_plot <- raw_plot+geom_edge_arc(aes(size = links[, edge_size], color = links[, edge_color]))
  }else{
    raw_plot <- raw_plot+geom_edge_arc(edge_width = edge_size, aes(color = links[,edge_color]))
  }

  if(!is.numeric(node_size)){
    raw_plot <- raw_plot+geom_node_point(aes(size = nodes[,node_size], color = node_color))
  }else{
    raw_plot <- raw_plot+geom_node_point(size = node_size, color = node_color)
  }
  out <- raw_plot+geom_node_text(aes(label = id_name), size = label_size, color = label_color,repel = T)+theme_void()
  return(out)
}

.translink <- function(matrix){
  tmp <- as.matrix(matrix);
  diag(tmp) <- NA;
  tmp[lower.tri(tmp)] <- NA;
  tmp_res <- reshape2::melt(tmp)
  colnames(tmp_res) <- c("From", "To", "Weight")
  tmp_res <- tmp_res[!is.na(tmp_res$Weight), ]
  tmp_res <- tmp_res[tmp_res$Weight !=0, ]
  if(nrow(tmp_res) == 0){
    stop("no decet edge!")
  }else{
    tmp_res$type <- ifelse(tmp_res$Weight >0, "+", "-")
    return(tmp_res)
  }

}
