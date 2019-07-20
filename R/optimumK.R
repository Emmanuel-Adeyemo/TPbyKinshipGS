
#' June 12, 2019
#' Optimum number of k
#'
#' @description
#' Finds the optimum number of clusters to be used for k-means clustering. The method used id
#' the within sum of squares - wss -
#'
#' @param geno_distance A dataframe of euclidian distance among lines
#' @param intercept_line A number given by user. Might be arbitary at first and the guided after
#'
#' @import ggplot2
#' @import ggfortify
#' @import factoextra
#' @import NbClust
#'
#'
#' @return
#'
#' @export
#'

get.optimum.k <- function(geno_distance, intercept_line){

  optimum_k <- fviz_nbclust(geno_distance, kmeans, method = "wss") +
    geom_vline(xintercept = intercept_line, linetype = 2)+  # include intercept number after optimum number has been determined
    labs(subtitle = "Elbow method")

  return(optimum_k)
}



#' June 14, 2019
#' K-means clustering
#'
#' @description
#' Stratified sampling with k-means method. Divides data into specified group
#' based on the number of k's determined by the algorithm. Number of k can also
#' be specified by user.
#'
#' @param geno_distance A dataframe of euclidian distance among lines
#' @param number_of_clusters A number given by user. The of clusters in the final output
#'
#' @import ggplot2
#' @import ggfortify
#' @import factoextra
#'
#'
#' @return
#'
#' @export
#'

get.kmeans <- function(geno_distance, number_of_clusters){

  kmeans_out <- kmeans(geno_distance, number_of_clusters, nstart=2, iter.max=1000, algorithm="Hartigan-Wong") # options:algorith="Hartigan-Wong", "Lloyd", "Forgy"

  return(kmeans_out)
}



#' June 14, 2019
#' Plot k-means
#'
#' @description
#' Stratified sampling with k-means method. Divides data into specified group
#' based on the number of k's determined by the algorithm. Number of k can also
#' be specified by user.
#'
#' @param cluster_with_pca A dataframe of lines with the corresponding pca and cluster assignment
#' @param cluster_number Column name with cluster assignment to be plotted e.g. cluster 2 is when k-means number = 2
#'
#' @import ggplot2
#' @import ggfortify
#' @import factoextra
#'
#'
#' @return
#'
#' @export
#'

plot.kmeans <- function(cluster_with_pca, cluster_number){

  cluster_with_pca[,cluster_number] <- as.factor(cluster_with_pca[,cluster_number])

  kmeans_plot <- ggplot(cluster_with_pca, aes(x=PC1, y=PC2, shape=cluster_with_pca[,cluster_number], color=cluster_with_pca[,cluster_number])) +
    geom_point()

  kmeans_bw <- kmeans_plot + theme_bw()

  return(kmeans_bw)
}


