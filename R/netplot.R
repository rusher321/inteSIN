#' network plot for single sample network
#'
#' @param matrix
#' @param layout
#' @param node_size
#' @param node_color
#' @param edge_size
#' @param edge_color
#' @param label_size
#' @param label_color
#' @param metadata
#'
#' @return
#' @export
#'
#' @examples
netplot_sin <- function(matrix, layout = "sphere", node_size = 3, node_color = "orange",
                        edge_size =1, edge_color= "gray", label_size = 3, label_color = "grey50",
                        metadata = NULL){

  ## this function based on ggraph

  links <- .translink(matrix) # get the link
  links$Weight2 <- abs(links$Weight)
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
    raw_plot <- raw_plot+geom_edge_arc(aes(size = Weight2, color = links[, edge_color]))
  }else{
    raw_plot <- raw_plot+geom_edge_arc(size = edge_size, aes(color = links[,edge_color]))
  }

  if(!is.numeric(node_size)){
    raw_plot <- raw_plot+geom_node_point(aes(size = nodes[,node_size], color = node_color))
  }else{
    raw_plot <- raw_plot+geom_node_point( size = node_size, color = node_color)
  }


  out <- raw_plot+geom_node_text(aes(label = id_name), size = label_size, color = label_color,repel = T)+theme_void()
  return(out)
}

.translink <- function(matrix){
  tmp <- as.matrix(matrix);
  diag(tmp) <- NA;
  tmp[lower.tri(tmp)] <- NA;
  tmp_res <- melt(tmp)
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
