#' @title Generate Latent Processes for Functional Data
#' 
#' @aliases rprocKmu
#' @keywords cluster
#' @rdname rprocKmu
#' 
#' @description Generates multiple latent processes (\code{mu.1}, \code{mu.2}, ..., \code{mu.k}) 
#' based on the specified input parameters, including process type, coefficients, exponents, 
#' and the number of evaluation points.
#' 
#' @param t A numeric vector specifying the grid points (\code{argvals}) at which 
#' the latent processes are evaluated. 
#' 
#' @param process A character vector specifying the types of the latent processes. 
#' Available options include \code{"poly1"}, \code{"poly2"}, \code{"sin"}, 
#' \code{"cos"}, \code{"sin_cos"}, and \code{"cos_sin"}.
#' 
#' @param c A numeric vector of coefficients used to scale each of the latent processes. 
#' 
#' @param k An optional numeric vector of exponents used in the latent processes. 
#' These exponents control the shape and curvature of the processes, 
#' and are only applicable for \code{"poly1"} and \code{"poly2"} process types.
#' 
#' @details 
#' The available latent process types are described as follows: 
#' \itemize{
#'   \item \code{"poly1"}: \eqn{\mu_i(t) = c_i \cdot (1 - t) \cdot t^{k_i}}, 
#'   where \code{c.i} is the scaling coefficient and \code{k.i} controls the curvature.
#'   
#'   \item \code{"poly2"}: \eqn{\mu_i(t) = c_i \cdot (1 - t)^{k_i} \cdot t}, 
#'   where \code{c.i} is the scaling coefficient and \code{k.i} controls the curvature.
#'   
#'   \item \code{"sin"}: \eqn{\mu_i(t) = c_i \cdot \sin(t)}, 
#'   where \code{c.i} is the scaling coefficient that adjusts the amplitude.
#'   
#'   \item \code{"cos"}: \eqn{\mu_i(t) = c_i \cdot \cos(t)}, 
#'   where \code{c.i} is the scaling coefficient that adjusts the amplitude.
#'   
#'   \item \code{"sin_cos"}: \eqn{\mu_i(t) = c_i \cdot (\sin(t) - \cos(t))}, 
#'   where \code{c.i} is the scaling coefficient.
#'   
#'   \item \code{"cos_sin"}: \eqn{\mu_i(t) = c_i \cdot (\cos(t) - \sin(t))}, 
#'   where \code{c.i} is the scaling coefficient.
#' }
#' 
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{mu}: A list of \code{fdata} objects, each representing one of the 
#'   latent processes.
#' }
#' 
#' @examples
#' t <- seq(0, 1, length.out = 101)
#' fproces <- rprocKmu(t, c("poly1", "poly2", "sin", "cos"), 
#'                     c(10, -10, 1, 1), c(1, 1, NA, NA))
#' plot(fproces, main = "Latent Processes")                   
#' 
#' @export
rprocKmu <- function(t, process, c, k = NULL) {
  # Validar la entrada
  k <- if (is.null(k)) rep(NA, length(process)) else k
  if (length(process) != length(c) || length(process) != length(k)) 
    stop("Lengths of 'process', 'c', and 'k' must be the same.")
  
  k[is.na(k)] <- 0  # Si no se proporciona k, se pone a 0
  
  # Inicializar la lista de procesos latentes
  # mu_list <- list()
  
  # Generar los procesos
  for (i in seq_along(process)) {
    mu <- switch(process[i],
                 "poly1" = {
                   c[i] * (1 - t) * t^k[i]
                 },
                 "poly2" = {
                   c[i] * (1 - t)^k[i] * t
                 },
                 "sin" = {
                   c[i] * sin(t)
                 },
                 "cos" = {
                   c[i] * cos(t)
                 },
                 "sin_cos" = {
                   c[i] * (sin(t) - cos(t))
                 },
                 "cos_sin" = {
                   c[i] * (cos(t) - sin(t))
                 },
                 stop(paste("Invalid process type for process", i, ":", 
                            process[i]))
    )
    if (i==1) fmu <- fdata(mu, t,
                names = list(main = process[i], xlab = process[i]))
    else fmu <-  c(fmu,fdata(mu,t))                        
    #mu_list[[i]] <- fdata(mu, t,
    #                      names = list(main = process[i], xlab = process[i]))
  }
  #return(list(mu = fmu))
  fmu
}
