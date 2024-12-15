#' @title K-Means Clustering for Multivariate Functional Data
#'
#' @description This function performs k-means clustering on multivariate functional datasets using distances calculated with `metric.mfdata`.
#' 
#' @param mfdata A list of `fdata` objects, or an `mfdata` object (a list of `fdata` objects) representing multiple functional datasets.
#' @param ncl Number of clusters.
#' @param max.iter Maximum number of iterations.
#' @param tol Tolerance for convergence.
#' 
#' @return A list with the following components:
#' \describe{
#'   \item{\code{cluster}}{A vector of cluster assignments for each observation.}
#'   \item{\code{centers}}{A list of centroids for each cluster for each functional variable.}
#'   \item{\code{iter}}{Number of iterations until convergence.}
#' }
#' 
#' @examples
#' data(aemet, package = "fda.usc")
#' datos <- mfdata("temp"=aemet$temp,"logprec"=aemet$logprec)
#' dd = metric.mfdata(datos)
#' result <- mfkmeans(datos, ncl = 3)
#' plot(datos,col=result$cluster)
#' # plot(aemet$df[,7:8],col=result$cluster,asp=T)
#' 
#' @rdname mfkmeans 
#' @aliases mfkmeans 
#' @keywords cluster
#' 
#' @export
mfkmeans <- function(mfdata, ncl = 2, max.iter = 100, tol = 1e-4) {
#  mfdata=datos
  #  ncl = 2
  #  max.iter = 100
  #tol = 1e-4
  
  n <- nrow(mfdata[[1]]$data) # Number of observations
  n_vars <- length(mfdata)    # Number of functional variables
  
  # Initialize cluster labels randomly
  cluster_labels <- sample(1:ncl, n, replace = TRUE)
  prev_cluster_labels <- rep(0, n)
  iter <- 0
  disNA <- matrix(NA, nrow = n, ncol = ncl)
  # Main loop for k-means
  while (!all(cluster_labels == prev_cluster_labels) && iter < max.iter) {
    iter <- iter + 1
    prev_cluster_labels <- cluster_labels
    
    # 1. Calculate centroids for each cluster
    centroids <- list()
    for (k in 1:ncl) {
      indices <- which(cluster_labels == k)
      if (length(indices) > 0) {
        # Calculate the centroid for each functional variable
        centroids[[k]] <- lapply(mfdata, function(fdata) func.mean(fdata[indices, ]))
      } else {
        # If no points in cluster, randomly select a point as the centroid
        random_index <- sample(1:n, 1)
        centroids[[k]] <- lapply(mfdata, function(fdata) fdata[random_index, ])
      }
    }
    
    # 2. Calculate distances from each observation to each cluster centroid
    distances <- disNA
    for (k in 1:ncl) {
      temp_centroid_mfdata <- list()
      for (i in seq_along(mfdata)) {
        temp_centroid_mfdata[[i]] <- c(centroids[[k]][[i]],centroids[[k]][[i]])
      }
      names(temp_centroid_mfdata) <- names(mfdata)
      distances[, k] <- metric.mfdata(mfdata, 
                                      temp_centroid_mfdata,method="euclidean")[,1]
    }
    
    # 3. Assign each observation to the closest cluster
    cluster_labels <- max.col(-distances, ties.method = "random")
  }
  
  # Return results
  list(
    cluster = cluster_labels,
    centers = centroids,
    iter = iter
  )
}
