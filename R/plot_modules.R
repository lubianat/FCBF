#' plot_one_module
#'
#' Takes as input an adjacency matrix from the find_modules function
#' (output = "matrix") an a list of hubs (FCBF - selected)
#' (Networks are inspired in the CEMiTool package networks)
#' And returns a pdf plotting the modules with their connections
#' @param adjacency_matrix An adjacency matrix from the find modules function
#' (dataframe with a 'genes' column)
#' @param genes A character vector with the FCBF selected genes
#' @param name The name of the plot.
#' @param color The color of the plot. Defaults to a shade of red ("#B30000FF").
#' @param n The max number of gene labels to show.  Defaults to 10.
#' @return A gg object
plot_one_module <- function(adjacency_matrix,
                             name,
                             color = "#B30000FF",
                             n = 10) {
  adj <- as.matrix(adjacency_matrix)

  ig_obj <- graph.adjacency(adj, weighted = TRUE)

  degrees <- igraph::degree(ig_obj, normalized = FALSE)

  ig_obj <-
    igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
  net_obj <- intergraph::asNetwork(ig_obj)
  m <-
    network::as.matrix.network.adjacency(net_obj) # get sociomatrix
  # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
  plotcord <-
    data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))
  # or get it them from Kamada-Kawai's algorithm:
  # plotcord <- data.frame(sna::gplot.layout.kamadakawai(m, NULL))
  colnames(plotcord) <- c("X1", "X2")
  edglist <- network::as.matrix.network.edgelist(net_obj)
  edges <-
    data.frame(plotcord[edglist[, 1], ], plotcord[edglist[, 2], ])
  plotcord$vertex.names <-
    as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
  plotcord$Degree <-
    network::get.vertex.attribute(net_obj, "degree")
  plotcord[, "shouldLabel"] <- FALSE

  max_n <- min(n, length(degrees))
  int_hubs <- names(sort(degrees, decreasing = TRUE))[seq_len(max_n)]
  int_bool <- plotcord[, "vertex.names"] %in% int_hubs
  sel_vertex <- int_hubs
  colnames(edges) <-  c("X1", "Y1", "X2", "Y2")
  #edges$midX  <- (edges$X1 + edges$X2) / 2
  #edges$midY  <- (edges$Y1 + edges$Y2) / 2
  plotcord[which(plotcord[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <-
    TRUE
  plotcord$Degree_cut <-
    cut(plotcord$Degree,
        breaks = 3,
        labels = FALSE)
  plotcord$in_mod <- TRUE
  pl <- ggplot(plotcord)  +
    geom_segment(
      data = edges,
      aes_(
        x =  ~ X1,
        y =  ~ Y1,
        xend =  ~ X2,
        yend =  ~ Y2
      ),
      size = 0.5,
      alpha = 0.7,
      colour = "gray25"
    ) +
    geom_point(aes_(
      x =  ~ X1,
      y =  ~ X2,
      size =  ~ Degree,
      alpha =  ~ Degree
    ),
    color = color) +
    geom_label_repel(
      aes_(
        x =  ~ X1,
        y =  ~ X2,
        label =  ~ vertex.names
      ),
      box.padding = unit(1, "lines"),
      data = function(x) {
        x[x$shouldLabel,]
      }
    ) +
    scale_colour_manual(values = c("#005E87")) +
    labs(title = name) +
    ggplot2::theme_bw(base_size = 12, base_family = "") +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white",
                                               colour = NA),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  return(pl)
}



