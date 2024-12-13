#' Functional DBSCAN Optimization and Clustering
#'
#' @description Perform DBSCAN clustering on functional data and optimize parameters `eps` and `minPts`.
#' @aliases fdbscan optim.fdbscan
#' @keywords cluster
#' @param fdataobj An object of class `fdata` containing functional data.
#' @param eps Neighborhood parameter (`eps`) for DBSCAN. If NULL, it is estimated automatically.
#' @param minPts Minimum cluster size (`minPts`). If NULL, it is estimated automatically.
#' @param metric Metric function to compute distances. Default is `metric.lp`.
#' @param par.metric List of arguments for the metric function.
#' @return A list with the following elements:
#' \itemize{
#'   \item `optimal`: Data frame with the best parameters (`eps`, `minPts`, `quality`).
#'   \item `model`: DBSCAN clustering model with refined clusters.
#'   \item `results`: Data frame with all combinations of `eps` and `minPts` tested.
#' }
#' @examples
#' \dontrun{
#' t <- seq(0, 2 * pi, length.out = 101)
#' res <- rprocKclust(t, n = c(30, 50, 40), process = c("cos_sin", "sin", "cos"),
#'                    c = c(-1, 1, 1), k = c(NA, NA, NA), s = c(0.2, 0.3, 0.1))
#' opt_results <- optim.fdbscan(res$X, metric = metric.lp, par.metric = list(lp = 2))
#' print(opt_results$optimal)
#' plot(res$X, col = opt_results$model$cluster+1, main = "Optimal DBSCAN Clustering")
#' }
#' @export
fdbscan <- function(fdataobj, eps, minPts = 5, metric = metric.lp, par.metric = list(lp = 2)) {
  # Configuración de parámetros de la métrica
  aaf <- formals(metric)[-c(1, 2)]
  if (length(par.metric) > 0) {
    nam <- names(par.metric)
    for (i in seq_along(par.metric)) {
      aaf[[nam[i]]] <- par.metric[[i]]
    }
  }
  
  # Cálculo de la matriz de distancias
  D <- do.call(metric, c(list(fdataobj), aaf))
  diag(D) <- NA  # Ignorar la diagonal
  n <- nrow(D)
  
  # Inicialización
  cl <- 0
  cluster <- rep(0, n)
  mark <- rep(NA, n)
  vecinos <- apply(D, 1, function(x, e) sum(x <= e, na.rm = TRUE), e = eps)
  mark[vecinos < 1] <- "N"  # Puntos aislados
  mark[vecinos >= minPts] <- "C"  # Puntos núcleo
  mark[vecinos < minPts & vecinos > 0] <- "A"  # Puntos borde
  
  # Verificar si hay puntos núcleo
  nl <- sum(mark == "C" & cluster == 0)
  if (nl == 0) stop("No center points. Consider modifying minPts or eps")
  
  # Algoritmo principal
  while (nl > 0) {
    lori <- which(mark == "C" & cluster == 0)[1]
    cl <- cl + 1
    cluster[lori] <- cl
    
    while (TRUE) {
      nc1 <- sum(cluster == cl)
      if (nc1 == 1) {
        lvec <- which(D[which(cluster == cl & mark == "C"), ] < eps)
      } else {
        lvec <- which(apply(D[which(cluster == cl & mark == "C"), , drop = FALSE], 2,
                            min, na.rm = TRUE) < eps)
      }
      cluster[lvec] <- cl
      
      nc0 <- sum(cluster == 0)
      if (nc0 == 0) break
      
      mDgr <- min(apply(D[which(cluster == cl & mark == "C"), 
                          which(cluster == 0), drop = FALSE], 2, min, na.rm = TRUE))
      if (mDgr >= eps) break
    }
    nl <- sum(mark == "C" & cluster == 0)
  }
  
  # Salida
  return(list(cluster = cluster, mark = mark))
}

#' @export
optim.fdbscan <- function(fdataobj, eps = NULL, minPts = NULL, metric = metric.lp, par.metric = list(lp = 2)) {
  # Get the number of points
  n <- if (is.matrix(fdataobj$data)) nrow(fdataobj$data) else nrow(fdataobj)
  
  # Automatically estimate eps and minPts if not provided
  if (is.null(minPts)) {
    sqrt_n <- sqrt(n)
    minPts <- round(seq(sqrt_n / 4, sqrt_n * 2, length.out = 7))
  }
  if (is.null(eps)) {
    eps <- fkNNdistplot(fdataobj = fdataobj, minPts = sqrt(n), metric = metric, par.metric = par.metric)
    eps <- seq(eps / 3, eps * 2, length.out = 7)
  }
  
  # Create grid of parameter combinations
  grid <- expand.grid(eps = eps, minPts = minPts)
  
  # Evaluate each combination
  results <- lapply(1:nrow(grid), function(i) {
    eps_val <- grid[i, "eps"]
    minPts_val <- grid[i, "minPts"]
    
    # Run DBSCAN with error handling
    result <- tryCatch(
      fdbscan(fdataobj, eps = eps_val, minPts = minPts_val),
      error = function(e) NA
    )
    
    # Evaluate clustering quality
    if (is.na(result)[1]) {
      quality <- NA
    } else {
      quality <- fclust.quality(fdataobj, result$cluster, metric, par.metric)
    }
    
    return(list(eps = eps_val, minPts = minPts_val, quality = quality, model = result))
  })
  
  # Convert results to data frame
  results_df <- do.call(rbind, lapply(results, function(x) data.frame(eps = x$eps, minPts = x$minPts, quality = x$quality)))
  
  # Check for valid results
  if (all(is.na(results_df$quality))) {
    warning("No valid clustering results found.")
    return(list(optimal = NA, model = NA, results = results_df))
  }
  
  # Select the best combination
  best_idx <- which.max(results_df$quality)
  best_result <- results_df[best_idx, ]
  
  # Fit the optimal model
  optimal_model <- results[[best_idx]]$model
  
  # Refine clusters
  refined_model <- refine_clusters(fdataobj, optimal_model, eps = best_result$eps, metric = metric, par.metric = par.metric)
  
  return(list(optimal = best_result, model = refined_model, results = results_df))
}

# Refine clusters to reassign noise points
refine_clusters <- function(fdataobj, model, eps, metric, par.metric = list(lp = 2)) {
  if (is.na(model)[1]) return(model)
  
  # Calculate distance matrix
  D <- do.call(metric, c(list(fdataobj), par.metric))
  diag(D) <- NA
  
  # Reassign noise points to nearest cluster
  for (i in which(model$cluster == 0)) {
    dists <- D[i, ]
    nearest_cluster <- which.min(tapply(dists, model$cluster, mean, na.rm = TRUE))
    if (!is.na(nearest_cluster) && dists[nearest_cluster] <= eps) {
      model$cluster[i] <- nearest_cluster
    }
  }
  return(model)
}

# Evaluate clustering quality
fclust.quality <- function(fdataobj, clusters, metric, par.metric = list(lp = 2)) {
  aaf <- formals(metric)[-c(1, 2)]
  if (length(par.metric) > 0) {
    nam <- names(par.metric)
    for (i in seq_along(par.metric)) {
      aaf[[nam[i]]] <- par.metric[[i]]
    }
  }
  
  # Calculate distance matrix
  D <- do.call(metric, c(list(fdataobj), aaf))
  diag(D) <- 0
  
  # Calculate silhouette scores
  unique_clusters <- unique(clusters[clusters > 0])
  if (length(unique_clusters) < 2) return(NA)
  
  silhouette_scores <- sapply(1:nrow(D), function(i) {
    if (clusters[i] == 0) return(NA)
    
    intra_cluster <- mean(D[i, clusters == clusters[i]], na.rm = TRUE)
    inter_cluster <- sapply(unique_clusters[unique_clusters != clusters[i]], function(c) {
      mean(D[i, clusters == c], na.rm = TRUE)
    })
    inter_cluster <- min(inter_cluster)
    
    (inter_cluster - intra_cluster) / max(intra_cluster, inter_cluster)
  })
  
  mean(silhouette_scores, na.rm = TRUE)
}

# Estimate optimal eps using kNN distances
fkNNdistplot <- function(fdataobj, minPts = NULL, metric = metric.lp, par.metric = list(lp = 2)) {
  if (is.null(minPts)) minPts <- select_minPts(fdataobj)
  
  aaf <- formals(metric)[-c(1, 2)]
  if (length(par.metric) > 0) {
    nam <- names(par.metric)
    for (i in seq_along(par.metric)) {
      aaf[[nam[i]]] <- par.metric[[i]]
    }
  }
  
  # Calculate distance matrix
  D <- do.call(metric, c(list(fdataobj), aaf))
  diag(D) <- NA
  
  # Compute kNN distances
  k_dist <- apply(D, 1, function(x) sort(x, na.last = NA)[minPts])
  eps <- fselect_eps(k_dist)
  
  plot(sort(k_dist), type = "l", main = "kNN Distance Plot", xlab = "Points sorted by k-distance",
       ylab = paste(minPts, "-NN Distance", sep = ""), col = "blue", lwd = 2)
  abline(h = eps, col = "red", lty = 2, lwd = 2)
  return(eps)
}

# Find the optimal eps using the elbow method
fselect_eps <- function(k_dist) {
  k_dist_sorted <- sort(k_dist, decreasing = TRUE)
  n <- length(k_dist_sorted)
  line <- cbind(1:n, k_dist_sorted)
  first_point <- line[1, ]
  last_point <- line[n, ]
  
  distances <- apply(line, 1, function(point) {
    abs((last_point[2] - first_point[2]) * point[1] - (last_point[1] - first_point[1]) * point[2] +
          last_point[1] * first_point[2] - last_point[2] * first_point[1]) /
      sqrt((last_point[2] - first_point[2])^2 + (last_point[1] - first_point[1])^2)
  })
  
  eps_idx <- which.max(distances)
  return(k_dist_sorted[eps_idx])
}

# Automatically select minPts
select_minPts <- function(fdataobj) {
  n <- if (is.matrix(fdataobj$data)) nrow(fdataobj$data) else nrow(fdataobj)
  max(5, round(sqrt(n)))
}




# hacer una funcion silueta y otras medidas

# optics
# kNNdistplot
# extractDBSCAN
