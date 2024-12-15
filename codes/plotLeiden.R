

plotLeiden <- function(result, Alus=FALSE, alu_vec=NULL) {
  # Ensure gdf is an igraph object
  if (!inherits(gdf, "igraph")) {
    stop("Input `gdf` must be an igraph object.")
  }
  
  # Assign community memberships to nodes in the graph
  V(gdf)$membership <- result$membership
  
  E(gdf)$scaled_weight <- scales::rescale(E(gdf)$weight, to = c(1, 15))
  
  # Create colormap for the communities
  num_communities <- length(unique(result$membership))
  community_colors <- rainbow(num_communities)
  colors <- community_colors[as.factor(V(gdf)$membership)]
  
  # Plot logic
  plot_ <- function() {
    plot(
      gdf,
      layout = layout_with_fr(gdf),
      vertex.color = colors,
      vertex.size = 4,
      vertex.frame.color = colors,
      edge.color = "grey",  # Default edge color
      edge.width = E(gdf)$scaled_weight,  # Use scaled weights
      vertex.label = NA,
      edge.arrow.size = 0,
      main = ifelse(Alus, "Candidate Chromatin Modulating Alu Communities (Chr22)", "Our Rcpp-based Clustering")
    )
    if (!Alus) {
      legend("topright", legend = unique(result$membership),
             col = community_colors, pch = 19,
             title = "Communities", cex = 0.8)
    }
  }
  
  # Execute the plot
  png("alu_plot.png", width = 2000, height = 2000, res = 300)
  plot_()
  dev.off()
}