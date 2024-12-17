setwd("/home/xuyuan/Desktop/2024 fall/project_data/original_data")
current_dir <- getwd()
R_source_dir <- gsub("original_data", "codes", current_dir)
source(file.path(R_source_dir, "clean_alu.R"))
source(file.path(R_source_dir, "merge_hic_graph.R"))

# store the R file to one folder and store the working directory to another folder
library(data.table)
library(tidyverse)
library(readxl)
library(R.utils)
# we first clean out the ATAC data and prepare the data for the analysis

# we write a function to clean the data and to prepare the data for the analysis, we will construct by the following operations
atac_file <- "atac_data.xlsx"
atac_data <- read_excel(atac_file)
output_csv_name <- "atac_data.csv"
current_dir <- getwd()
new_dir <- gsub("original_data", "temp", current_dir)
output_csv_name <- file.path(new_dir, "atac_data.csv")
write.csv(atac_data, output_csv_name, row.names = FALSE)

# we will use the following code to unzip the file
input_file <- "hg38.hits.gz"
new_dir <- gsub("original_data", "temp", current_dir)
unzip_output_file <- file.path(new_dir, "hg38.hits")
gunzip(input_file, destname = unzip_output_file, remove = FALSE)

# now we unzip the file of the fitHiC data
input_file <- "GSE148434_NS_all.5kb.1e-2.FitHiC.all.txt.gz"
new_dir <- gsub("original_data", "temp", current_dir)
unzip_output_file <- file.path(new_dir, "hic_data.txt")
gunzip(input_file, destname = unzip_output_file, remove = FALSE)

# now we have all the data in the temp folder, we will now process the data
data_dir <- gsub("original_data", "temp", current_dir)
# be careful that the running time is large because we have a lot of data to process
# after calcualtion, we will processed all files by chromosome
returned_path <- clean_alu(
  paste0(data_dir, "/hic_data.txt"),
  paste0(data_dir, "/hg38.hits"),
  data_dir
)

# now we can generate the edges
output_directory <- gsub("original_data", "cleaned", current_dir)

cleaning_dir <- merge_hic_with_alu(data_dir, output_directory)







dataframe <- read.csv("edges_data_frame.csv")
dataframe$group <- NULL

# now we consider the group of the data by using the Leiden algorithm
# we first start using the chromosome 21
data <- dataframe[which(dataframe$chromosome == "chr21"), ]

# print(nrow(dataframe))
# print(summary(dataframe$node_weights))
## we first try to give each weight a very small value to let every data in the model
# dataframe$node_weights <- dataframe$node_weights + 0.01


runIgraph <- function(gdf, iterations, gamma, theta) {
  ldc <- cluster_leiden(
    gdf,
    resolution = gamma,  # Correct argument
    objective_function = "CPM",
    weights = E(gdf)$weight,  # Use edge weights from the graph
    beta = theta,
    n_iterations = iterations,
    vertex_weights = rep(1, vcount(gdf))  # Use vertex count from the graph
  )
  return(ldc)
}

library(igraph)

gdf <- graph_from_data_frame(data, directed = FALSE)
E(gdf)$weight <- as.numeric(data$node_weights)
set.seed(1)
result <- runIgraph(gdf, 50, 0.0005, 0.0001)
print(length(unique(result$membership)))
# we can now plot the graph
# Match data to graph nodes
matched_data <- data[match(V(gdf)$name, data$X), ]

# Assign the dominant family or "Other" if counts are equal
matched_data$dominant_family <- apply(
  matched_data[, c("AluJ_count", "AluS_count", "AluY_count")], 
  1, 
  function(row) {
    if (length(unique(row[!is.na(row)])) == 1 || all(is.na(row)) || all(row == 0)) {
      return(4)  # Assign to "Other" if all counts are equal
    } else {
      return(which.max(row))  # Assign to the max family otherwise
    }
  }
)
matched_data$dominant_family <- as.numeric(matched_data$dominant_family)
# Add the dominant family to the graph
V(gdf)$alu_family <- factor(
  matched_data$dominant_family, 
  levels = 1:4, 
  labels = c("AluJ", "AluS", "AluY", "Other")
)

print(unique(V(gdf)$alu_family))
plotLeiden <- function(result, gdf, Alus = FALSE) {
  # Ensure gdf is an igraph object
  if (!inherits(gdf, "igraph")) {
    stop("Input `gdf` must be an igraph object.")
  }
  
  # Assign community memberships to nodes in the graph
  V(gdf)$membership <- result$membership
  
  # Normalize edge weights for better visualization
  if (!is.null(E(gdf)$weight)) {
    E(gdf)$scaled_weight <- scales::rescale(E(gdf)$weight, to = c(1, 15))
  } else {
    E(gdf)$scaled_weight <- rep(1, ecount(gdf))  # Default weight if none exists
  }
  
  # Create colormap for the communities
  num_communities <- length(unique(result$membership))
  community_colors <- rainbow(num_communities)
  colors <- community_colors[as.factor(V(gdf)$membership)]
  
  # Ensure alu_family attribute exists in the graph
  if (!"alu_family" %in% vertex_attr_names(gdf)) {
    stop("The graph does not contain the 'alu_family' attribute. Please assign it before plotting.")
  }
  
  # Map Alu families to colors
  alu_family_colors <- c("red", "blue", "green", "purple")  # Assign specific colors to AluJ, AluS, AluY
  alu_family_labels <- c("AluJ", "AluS", "AluY", "Other")
  alu_colors <- alu_family_colors[V(gdf)$alu_family]
  
  # Plot logic
  plot_ <- function() {
    plot(
      gdf,
      layout = layout_with_fr(gdf),
      vertex.color = alu_colors,  # Use Alu family colors
      vertex.size = 4,
      vertex.frame.color = alu_colors,
      edge.color = "grey",  # Default edge color
      edge.width = E(gdf)$scaled_weight,  # Use scaled weights
      vertex.label = NA,  # Suppress labels
      edge.arrow.size = 0,
      main = ifelse(Alus, "Candidate Chromatin Modulating Alu Communities (Chr22)", "Our Rcpp-based Clustering")
    )
    if (!Alus) {
      legend("topright", legend = alu_family_labels,
             col = alu_family_colors, pch = 19,
             title = "Alu Families", cex = 0.8)
    } else {
      # Add a legend to the top-right
      legend(
        "topright", 
        legend = alu_family_labels, 
        col = alu_family_colors, 
        pch = 19, 
        title = "Alu Families", 
        cex = 0.8
      )
    }
  }
  
  # Save the plot to a PNG file with 300 DPI
  png("alu_plot.png", width = 2000, height = 2000, res = 300)
  plot_()
  dev.off()
}

plotLeiden(result, gdf, Alus = TRUE)

# 拉普拉斯矩阵
laplacian_matrix <- graph.laplacian(gdf, normalized = FALSE)  # 非归一化
laplacian_matrix_normalized <- graph.laplacian(gdf, normalized = TRUE)  # 归一化

# 查看矩阵
print(laplacian_matrix)

# 将矩阵保存到文件
write.csv(as.matrix(laplacian_matrix), "laplacian_matrix.csv", row.names = FALSE)

# 计算拉普拉斯矩阵
laplacian_matrix <- graph.laplacian(gdf, normalized = TRUE)

# 计算特征值和特征向量
eigen_result <- eigen(as.matrix(laplacian_matrix))

# 提取二阶特征向量
fiedler_vector <- eigen_result$vectors[, 2]

# 将节点划分为两组
node_groups <- ifelse(fiedler_vector > 0, 1, 2)

# 打印节点分组
print(node_groups)

# 热力图可视化拉普拉斯矩阵
library(ggplot2)
library(reshape2)

laplacian_matrix <- graph.laplacian(gdf, normalized = FALSE)
laplacian_df <- as.data.frame(as.table(as.matrix(laplacian_matrix)))

colnames(laplacian_df) <- c("Row", "Column", "Value")
ggplot(laplacian_df, aes(x = Row, y = Column, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Laplacian Matrix Heatmap", x = "Nodes", y = "Nodes")


# 计算拉普拉斯矩阵和特征值/特征向量
laplacian_matrix <- graph.laplacian(gdf, normalized = TRUE)
eigen_result <- eigen(as.matrix(laplacian_matrix))
fiedler_vector <- eigen_result$vectors[, 2]

# 节点分组（根据二阶特征向量分区）
node_groups <- ifelse(fiedler_vector > 0, "Group 1", "Group 2")

# 可视化网络，节点颜色对应分组
library(igraph)
V(gdf)$group <- node_groups
plot(
  gdf,
  vertex.color = V(gdf)$group,
  vertex.size = 5,
  edge.color = "grey",
  main = "Network Partition Based on Fiedler Vector"
)

# 计算 Closeness Centrality
closeness_centrality <- closeness(gdf, mode = "all", weights = E(gdf)$weight)
V(gdf)$closeness <- closeness_centrality

# 可视化网络，节点大小对应中心性
plot(
  gdf,
  vertex.size = scales::rescale(V(gdf)$closeness, to = c(5, 15)),  # 中心性大小
  vertex.color = "skyblue",
  edge.color = "grey",
  vertex.label = NA,
  main = "Network Visualization with Closeness Centrality"
)

# 使用 Leiden 算法计算社区
library(igraph)
leiden_result <- cluster_leiden(gdf, resolution_parameter = 0.01)

# 添加社区分组到图中
V(gdf)$community <- leiden_result$membership

# 可视化网络，节点颜色对应社区
plot(
  gdf,
  vertex.color = V(gdf)$community,
  vertex.size = 5,
  edge.color = "grey",
  main = "Community Detection with Leiden Algorithm"
)

# 拉普拉斯矩阵的特征值分布
laplacian_matrix <- graph.laplacian(gdf, normalized = TRUE)
eigen_result <- eigen(as.matrix(laplacian_matrix))

# 绘制特征值分布
eigenvalues <- eigen_result$values
qplot(eigenvalues, geom = "histogram", bins = 30, fill = I("blue"), alpha = I(0.7)) +
  theme_minimal() +
  labs(title = "Eigenvalue Distribution of Laplacian Matrix", x = "Eigenvalues", y = "Frequency")

num_nodes <- vcount(gdf)
num_edges <- ecount(gdf)
sparsity <- 1 - (num_edges / (num_nodes * (num_nodes - 1) / 2))
print(paste("Network Sparsity:", round(sparsity, 4)))

degree_dist <- degree(gdf)
hist(degree_dist, breaks = 30, main = "Degree Distribution", xlab = "Degree", ylab = "Frequency")

# 计算连通分量
components <- components(gdf)
print(components$csize)  # 每个连通分量的大小

# 提取最大的连通分量
largest_comp <- induced_subgraph(gdf, which(components$membership == which.max(components$csize)))

# 分析最大连通分量
plot(largest_comp, main = "Largest Connected Component")
laplacian_matrix <- graph.laplacian(gdf, normalized = TRUE)
eigenvalues <- eigen(as.matrix(laplacian_matrix))$values
hist(eigenvalues, breaks = 50, main = "Eigenvalue Distribution", xlab = "Eigenvalue", ylab = "Frequency")

# 降维嵌入
laplacian_matrix <- graph.laplacian(gdf, normalized = TRUE)
eigen_result <- eigen(as.matrix(laplacian_matrix))
embedding <- eigen_result$vectors[, 1:2]  # 使用前两个特征向量

plot(embedding, col = V(gdf)$community, pch = 19, main = "Spectral Embedding of Sparse Network")


