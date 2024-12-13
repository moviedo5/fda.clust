#' @title Generate Latent Processes for Functional Data
#' 
#' @aliases rproc2mu
#' @keywords cluster
#' @rdname rproc2mu
#' 
#' @description Generates two latent processes (\code{mu.1} and \code{mu.2}) 
#' based on the specified input parameters, including process type, coefficients,
#'  exponents, and the number of evaluation points.
#' 
#' @param t A numeric vector specifying the grid points (\code{argvals}) at which 
#' the latent processes are evaluated. 
#' 
#' @param process1,process2 Character strings specifying the types of the first
#'  and second latent processes, respectively. Available options include 
#'  \code{"poly1"} (by default), \code{"poly2"}  (by default), 
#'  \code{"sin"}, and \code{"cos"}.
#' 
#' @param c.1,c.2 Numeric coefficients used to scale the first and second latent
#'  processes, respectively. 
#' 
#' @param k.1,k.2 Optional numeric exponents used in the first and second latent 
#' processes, respectively. These exponents control the shape and curvature of the processes, only applicable for \code{"poly1"}
#'  and \code{"poly2"} process types.
#' 
#' @details 
#' The available latent process types are described as follows: 
#' \itemize{
#'   \item \code{"poly1"}: A polynomial process of the form \eqn{\mu_1(t) = c_1 \cdot (1 - t) \cdot t^{k_1}}, 
#'   where \code{c.1} is the scaling coefficient and \code{k.1} controls the curvature of the process. 
#'   
#'   \item \code{"poly2"}: A polynomial process of the form \eqn{\mu_2(t) = c_2 \cdot (1 - t)^{k_2} \cdot t}, 
#'   where \code{c.2} is the scaling coefficient and \code{k.2} controls the curvature of the process. 
#'   
#'   \item \code{"sin"}: A sinusoidal process of the form \eqn{\mu_i(t) = c_i \cdot \sin(t)}, 
#'   where \code{c.i} is the scaling coefficient that adjusts the amplitude of the sinusoidal wave. 
#'   
#'   \item \code{"cos"}: A cosine process of the form \eqn{\mu_i(t) = c_i \cdot \cos(t)}, 
#'   where \code{c.i} is the scaling coefficient that adjusts the amplitude of the cosine wave. 
#'   
#'   \item \code{"sin_cos"}: A combination of sinusoidal and cosine processes 
#'   of the form \eqn{\mu_1(t) = c_i \cdot (\sin(t) - \cos(t))}, 
#'   where \code{c.i} is the scaling coefficient.
#'      
#'   \item \code{"cos_sin"}: A combination of cosine and sinusoidal processes 
#'   of the form \eqn{\mu_2(t) = c_i \cdot (\cos(t) - \sin(t))}, 
#'   where \code{c.i} is the scaling coefficient.
#' }
#' 
#' 
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{mu.1}: An \code{fdata} object representing the first latent process.
#'   \item \code{mu.2}: An \code{fdata} object representing the second latent process.
#' }
#' 
#' @examples
#' # Example 1: Generate latent processes using poly1 and poly1
#' np <- 101
#' t1 <- seq(0, 1, length.out = np)
#' fproces <- rproc2mu(t1, "poly1", "poly1", 10, -10, 1, 1)
#' plot(fproces$mu.1, col = "blue", 
#'       main = "Latent Processes",
#'       ylim=range(cbind(range(fproces$mu.1),range(fproces$mu.2))))
#' lines(fproces$mu.2, col = "red")
#' 
#' # Example 2: Generate latent processes using poly1 and poly2
#' fproces <- rproc2mu(t1, "poly1", "poly2", 1.1, 1, 1, 1)
#' plot(fproces$mu.1, col = "blue", main = "Latent Processes")
#' lines(fproces$mu.2, col = "red")
#' 
#' # Example 3: Generate latent processes using sin and cos
#' t2 <- seq(0, pi, len = np)
#' fproces <- rproc2mu(t2, "sin", "cos", 1, 1)
#' plot(fproces$mu.1, col = "blue", main = "Latent Processes")
#' lines(fproces$mu.2, col = "red")
#' 
#' # Example 4: Generate latent processes using sin and cos
#' t3 <- seq(0, 2 * pi, len = np)
#' fproces <- rproc2mu(t2, "cos", "cos_sin", 1, 1)
#' plot(fproces$mu.1, col = "blue", main = "Latent Processes")
#' lines(fproces$mu.2, col = "red")


#' @export
rproc2mu <- function(t, process1 = "poly1", process2 ="poly2",
                     c.1, c.2, k.1 = NULL, k.2 = NULL) {
  if (!is.character(process1) || !is.character(process2)) 
    stop("Both 'process1' and 'process2' must be character strings.")
  if (!is.numeric(c.1) || !is.numeric(c.2)) 
    stop("Both 'c.1' and 'c.2' must be numeric.")
  if (!is.numeric(t) )     
    stop("Parameter 't' must be a numeric vector.")
  if (!is.null(k.1) && !is.numeric(k.1)) 
    stop("Parameter 'k.1' must be numeric if provided.")
  if (!is.null(k.2) && !is.numeric(k.2)) 
    stop("Parameter 'k.2' must be numeric if provided.")
  
  # Calculate latent processes based on the specified process type
  # t <- seq(0, 1, length.out = np)
  mu.1 <- switch(process1,
                 "poly1" = {
                   group1 <- bquote(.(c.1)~"(1 - t)t"^.(k.1))
                   c.1 * (1 - t) * t^k.1
                 },
                 "poly2" = {
                   group1 <- bquote(.(c.1)~"(1 - t)"^.(k.1)~"t")
                   c.1 * (1 - t)^k.1 * t
                 },
                 "sin" = {
                   group1= "c.1 * sin(t)"
                   c.1 * sin(t)
                 },
                 "cos" = {
                   group1= "c.1 * cos(t)"
                   c.1 * cos(t)
                 },
                 "sin_cos" = {
                   group1= "c.2 * (sin(t) - cos(t))"
                   c.1 * (sin(t) - cos(t)) # sin(t) - cos(t)
                 },
                 "cos_sin" = {
                   group1= "c.1 * (cos(t) - sin(t))"
                   c.1 * (cos(t) - sin (t)) # cos(x)-sin(x)
                 },
                 stop("Invalid process type for process1")
  )
  
  mu.2 <- switch(process2,
                 "poly1" = {
                   group2 <- bquote(.(c.2)~"(1 - t)t"^.(k.2))
                   c.2 * (1 - t) * t^k.2
                 },
                 "poly2" = {
                   group2 <- bquote(.(c.2)~"(1 - t)"^.(k.2)~"t")
                   c.2 * (1 - t)^k.2 * t
                 },
                 "sin" = {
                  group2= "c.2 * sin(t)"
                   c.2 * sin(t)
                 },
                 "cos" = {
                   group2= "c.2 * cos(t)"
                   c.2 * cos(t)},
                "sin_cos" = {
                   group2= "c.2 * (sin(t) - cos(t))"
                   c.2 * (sin(t) - cos(t)) # sin(t) - cos(t)
                },
                 "cos_sin" = {
                   group2= "c.1 * (cos(t) - sin(t))"
                   c.2 * (cos(t) - sin (t)) # cos(x)-sin(x)
                 },
                 stop("Invalid process type for process2")
  )
  # convertir a mfdata?
  return(list(
    mu.1 = fdata(mu.1, t, names=list(main=process1, "xlab"= group1)),
    mu.2 = fdata(mu.2, t, names=list(main=process2, "xlab"= group2))
  ))
}
