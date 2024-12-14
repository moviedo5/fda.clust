
#' @title Functional Data Clustering Measures
#' @description This function computes clustering evaluation measures for functional data objects of class 'fdata'. It supports the calculation of the Dunn, Davies-Bouldin, Calinski-Harabasz, and Silhouette indices. Distances are computed using the metric.lp function from the fda.usc package by default, but users can specify a different distance metric via the 'metric' parameter. Cluster means are calculated using an appropriate function, which can be specified by the user via the 'center_func' argument.
#' 
#' @param X An fdata object containing the functional data. The rows represent observations and the columns represent discrete evaluations of the functional data.
#' @param clusters A vector containing the cluster assignments for each observation in the fdata object.
#' @param index A character string indicating the measure to compute. Possible options are "silhouette" (default), "dunn", "db" (Davies-Bouldin), or "ch" (Calinski-Harabasz).
#' @param metric A function specifying the distance metric to be used. The default is 'metric.lp'.
#' @param par.metric A list of parameters to be passed to the 'metric' function.
#' @param center_func A function to compute cluster means, with possible options like 'func.mean', 'func.trim.mode', or other user-defined functions.
#' 
#' @return The value of the selected clustering measure.
#' 
#' @examples 
#' set.seed(123)
#' t <- seq(0, 2 * pi, length.out = 101)
#' res <- rprocKclust(t, n = c(30, 50, 40), 
#'                    process = c("cos_sin", "sin", "cos"), 
#'                    c = c(-1, 1, 1), k = c(NA, NA, NA), 
#'                    s = c(0.2, 0.3, 0.1))
#' X <- res$X
#' clusters <- res$groups
#' silhouette_val <- fclust.measures(X, clusters, index = "silhouette")
#' dunn_val <- fclust.measures(X, clusters, index = "dunn")
#' db_val <- fclust.measures(X, clusters, index = "db")
#' ch_val <- fclust.measures(X, clusters, index = "ch")
#' print(silhouette_val)
#' print(dunn_val)
#' print(db_val)
#' print(ch_val)
#' 
#' @export
fclust.measures <- function(X, clusters, index = "silhouette", metric = metric.lp, par.metric = list(), center_func = func.mean) {
  if (!inherits(X, "fdata")) {
    stop("The argument 'X' must be an object of class 'fdata'.")
  }
  
  metric <- do.call(metric, c(list(X), par.metric))
  
  if (index == "silhouette") {
    return(silhouette_fdata(X, clusters, metric))
  } else if (index == "dunn") {
    return(dunn_index(X, clusters, metric))
  } else if (index == "db") {
    return(davies_bouldin(X, clusters, metric))
  } else if (index == "ch") {
    return(calinski_harabasz(X, clusters, center_func))
  } else {
    stop("The 'index' argument must be one of 'silhouette', 'dunn', 'db', or 'ch'.")
  }
}

# Dunn Index
dunn_index <- function(X, clusters, metric) {
  intercluster_min <- min(sapply(unique(clusters), function(k) {
    min(metric[clusters == k, clusters != k])
  }))
  
  intracluster_max <- max(sapply(unique(clusters), function(k) {
    max(metric[clusters == k, clusters == k])
  }))
  
  return(intercluster_min / intracluster_max)
}

# Davies-Bouldin Index
davies_bouldin <- function(X, clusters, metric) {
  unique_clusters <- unique(clusters)
  
  db_index <- mean(sapply(unique_clusters, function(k1) {
    s1 <- mean(metric[clusters == k1, clusters == k1])
    max(sapply(setdiff(unique_clusters, k1), function(k2) {
      s2 <- mean(metric[clusters == k2, clusters == k2])
      d12 <- min(metric[clusters == k1, clusters == k2])
      (s1 + s2) / d12
    }))
  }))
  
  return(db_index)
}

# Calinski-Harabasz Index
calinski_harabasz <- function(X, clusters, center_func = func.mean) {
  n <- nrow(X$data)
  k <- length(unique(clusters))
  
  # Calculate the overall mean
  overall_mean <- center_func(X)
  if (inherits(overall_mean, "fdata")) {
    overall_mean <- overall_mean$data
  }
  
  # Calculate the cluster means
  cluster_means <- lapply(unique(clusters), function(k) {
    cluster_data <- X[clusters == k, , drop = FALSE]
    cluster_mean <- center_func(cluster_data)
    if (inherits(cluster_mean, "fdata")) {
      cluster_mean <- cluster_mean$data
    }
    return(cluster_mean)
  })
  cluster_means <- do.call(rbind, cluster_means)
  cluster_sizes <- sapply(unique(clusters), function(k) sum(clusters == k))
  
  # Calculate SSB (sum of squares between clusters)
  SSB <- sum(sapply(1:k, function(j) {
    cluster_size <- cluster_sizes[j]
    cluster_mean <- cluster_means[j, , drop = FALSE]
    sum((cluster_mean - overall_mean)^2) * cluster_size
  }))
  
  # Calculate SSW (sum of squares within clusters)
  SSW <- sum(sapply(unique(clusters), function(k) {
    cluster_data <- X[clusters == k, , drop = FALSE]
    cluster_mean <- cluster_means[which(unique(clusters) == k), , drop = FALSE]
    sum(rowSums((cluster_data$data - matrix(rep(cluster_mean, nrow(cluster_data$data)), 
                                            nrow = nrow(cluster_data$data), 
                                            byrow = TRUE))^2))
  }))
  
  # Calculate the Calinski-Harabasz index
  CH_index <- (SSB / (k - 1)) / (SSW / (n - k))
  
  return(CH_index)
}

# Silhouette Index for Functional Data Clustering
silhouette_fdata <- function(fdata_obj, clusters, metric) {
  n <- nrow(fdata_obj$data)
  silhouette_scores <- numeric(n)
  
  for (i in 1:n) {
    current_cluster <- clusters[i]
    a_i <- mean(metric[i, clusters == current_cluster])
    b_i <- min(sapply(setdiff(unique(clusters), current_cluster), function(c) mean(metric[i, clusters == c])))
    silhouette_scores[i] <- (b_i - a_i) / max(a_i, b_i)
  }
  
  global_silhouette <- mean(silhouette_scores)
  return(global_silhouette)
}


