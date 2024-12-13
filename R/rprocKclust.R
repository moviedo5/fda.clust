#' @title Simulate Functional Data for K Clusters
#' 
#' @aliases rprocKclust
#' @keywords cluster
#' @rdname rprocKclust
#' 
#' @description Simulates functional data for K clusters based on specified parameters for 
#' clustering analysis within the `fda.clust` package. It uses `rprocKmu` 
#' to generate the latent processes and `rproc2fdata` from `fda.usc` to generate
#' the functional data.
#' 
#' @param t A numeric vector representing the evaluation points (\code{argvals}) 
#' at which the functional data are computed. 
#' 
#' @param n A vector of integers specifying the number of curves to be generated 
#' for each of the K clusters. 
#' 
#' @param process A character vector of length K specifying the type of each latent process.
#' Options include \code{"poly1"}, \code{"poly2"}, \code{"sin"}, \code{"cos"}, 
#' \code{"sin_cos"}, and \code{"cos_sin"}.
#' 
#' @param c A numeric vector of length K specifying the coefficients for each process.
#' 
#' @param k An optional numeric vector of length K specifying the exponents for the 
#' processes. This parameter is only relevant for \code{"poly1"} and \code{"poly2"}.
#' 
#' @param s A numeric vector specifying the standard deviation of the noise to be 
#' added to each cluster's data. If only one value is provided, it is recycled for all clusters.
#' 
#' @param par.list A list of lists of additional parameters to be 
#' passed to \code{rproc2fdata} for generating the functional data for each cluster. 
#' 
#' @param ... Additional arguments to be passed to \code{rproc2fdata}, allowing 
#' for further customization of the process generation.
#' 
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{X}: The functional data as an \code{fdata} object.
#'   \item \code{X.hat}: The estimated mean functions for each cluster.
#'   \item \code{groups}: The group labels for each curve.
#'   \item \code{color}: The colors associated with each cluster.
#'   \item \code{colors}: The colors assigned to each individual curve.
#' }
#' 
#' @examples
#' t <- seq(0, 2*pi, length.out = 101)
#' res <- rprocKclust(t, 
#'                    n = c(30, 50, 40), 
#'                    process = c("poly1", "sin", "cos"), 
#'                    c = c(10, 1, 1), 
#'                    k = c(2, NA, NA), 
#'                    s = c(0.2, 0.3, 0.1))
#' 
#' plot(res$X, col = res$colors)
#' lines(res$X.hat, lwd = 2, col = res$color)
#' 
#' @export
rprocKclust <- function(t, n, process, c, k = NULL, s = 1, par.list = NULL, ...) {
  
  # Paso 1: Validar la entrada
  K <- length(process) # Número de clusters (K)
  
  # Inicializar y completar k
  if (is.null(k)) {
    k <- rep(0, K) 
  } else if (length(k) != K) {
    stop("The length of 'k' must be equal to the length of 'process'.")
  }
  k[is.na(k)] <- 0  # Reemplazar NA con 0
  
  # Validar longitud de n, c y process
  if (length(n) == 1) n <- rep(n, K) 
  if (length(s) == 1) s <- rep(s, K)
  
  if (!(length(process) == length(c) && length(process) == length(k) && 
        length(process) == length(n) && length(process) == length(s))) {
    stop("Lengths of 'process', 'c', 'k', 'n', and 's' must be the same.")
  }
  
  # Inicializar par.list
  if (is.null(par.list)) {
    par.list <- replicate(K, list(), simplify = FALSE)
  } else if (length(par.list) != K) {
    stop("The length of 'par.list' must be equal to the length of 'process'.")
  }
  
  # Paso 2: Inicializar listas para almacenar los resultados
  X <- NULL
  X.hat <- NULL
  groups <- factor()
  color_palette <- grDevices::rainbow(K)  # Paleta de colores para los K clústeres
  colors <- c()
  
  # Generar los K procesos latentes
  mu_fdata <- rprocKmu(t, process, c, k)  # Devuelve una lista de objetos fdata
  
  # Paso 3: Generar los datos funcionales para cada grupo
  for (i in seq_along(process)) {
    ni <- n[i] # Cantidad de curvas del grupo i
    mui <- mu_fdata$data[i,]  # Extraer el vector de media para el clúster i
  #  if (length(mui) != length(t)) {
  #    stop(paste0("The length of t (", length(t), ") and mu (", length(mui), ") must be the same for cluster ", i, "."))
  #  }
    si <- s[i]  # Desviación estándar para el ruido
    parlisti <- par.list[[i]]  # Parámetros adicionales para rproc2fdata
    
    # Generar los datos funcionales para el grupo i
    if (isFALSE(parlisti) || is.null(parlisti)) {
      Xi <- rproc2fdata(ni, t, mu = mui, sigma = si, ...)
    } else {
      Xi <- rproc2fdata(ni, t, mu = mui, sigma = si,
                        par.list = parlisti, ...)
    }
    
    # Concatenar los resultados
    if (i ==1) {
      X=Xi
      X.hat=fda.usc::func.mean(Xi)
      } else {
      X  <- c(X,Xi)
      X.hat <- c(X.hat,fda.usc::func.mean(Xi))  # Usar la función `func.mean` explícitamente
      }
    groups <- c(groups, rep(i, ni))
    colors <- c(colors, rep(color_palette[i], ni))
  }
  

  # Paso 5: Retornar los resultados
  return(list(
    "X" = X, 
    "X.hat" = X.hat, 
    "groups" = factor(groups), 
    "color" = color_palette, 
    "colors" = colors
  ))
}

# t <- seq(0, 1, length.out = 101)
# res <- rprocKclust(t, 
#                    n = c(50, 40), 
#                    process = c("poly1", "poly2"), 
#                    c = c( 10, 10), 
#                    k = c(1,1.2), 
#                    s = c( 0.3, 0.1))
# names(res)
# plot(res$X,col=res$colors)
# lines(res$X.hat,col=c(1,2))
# 
# 
# t <- seq(0, 2 * pi, length.out = 101)
# res <- rprocKclust(t, 
#                    n = c(50, 40), 
#                    process = c("sin", "cos"), 
#                    c = c( 1, -1), 
#                    k = c(NA, NA), 
#                    s = c( 0.3, 0.1))
# names(res)
# plot(res$X,col=res$colors)
# lines(res$X.hat,col=c(1,2))
