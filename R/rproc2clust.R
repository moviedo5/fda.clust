#' @title Simulate Functional Data for Clustering
#' 
#' @aliases rproc2clust
#' @keywords cluster
#' @rdname rproc2clust
#' 
#' @description Simulates functional data based on specified parameters for 
#' clustering analysis within the `fda.clust` package. It uses `rproc2mu` 
#' to generate the latent processes and `rproc2fdata` from `fda.usc` to generate
#'  the functional data.
#'   
#' @param n.1,n.2 Integers specifying the number of curves to be generated 
#' in the first  and second groups, respectively. 
#' 
#' @param mu.1,mu.2 Objects of class \code{fdata} representing the theoretical 
#' mean functions for the first and second groups, respectively. 
#' 
#' @param s.1,s.2 Numeric values specifying the standard deviation of the noise 
#' to be added to the first and second groups, respectively. 
#' See the \code{sigma} argument of the \code{rproc2fdata} function.
#' 
#' @param par.list.1,par.list.2 Optional lists of additional parameters to be 
#' passed to \code{rproc2fdata} for generating the first and second groups, respectively. 
#' If set to \code{FALSE}, default parameter values will be used.
#' 
#' @param ... Additional arguments to be passed to \code{rproc2fdata}, allowing 
#' for further customization of the process generation.
#' 
#' @return A list containing the following elements:
#' \itemize{
#'   \item \code{X}: The functional data as an \code{fdata} object.
#'   \item \code{X.hat}: The estimated mean functions for each group.
#'   \item \code{group}: The theoretical functions for each group.
#'   \item \code{groups}: The group labels for each curve.
#'   \item \code{color}: The colors associated with each group.
#'   \item \code{colors}: The colors assigned to each individual curve.
#' }
#' 
#' @examples
#' fproces <- rproc2mu( 51, "sin", "cos", 1, 1)
#' res <- rproc2clust (10,50,fproces$mu.1,fproces$mu.2)
#' plot(res$X, col=res$groups)
#' lines(res$X.hat, lwd=2, col=res$color)
#' 
#' fproces <- rproc2mu( 51, "sin", "sin", 1, 2)
#' res <- rproc2clust (10,50,fproces$mu.1,fproces$mu.2)
#' plot(res$X, col=res$groups)
#' lines(res$X.hat, lwd=2, col=res$color)
#' 
#' fproces <- rproc2mu( 101, "poly1", "poly1", 10, -10, 1, 1)
#' res <- rproc2clust (10,50,fproces$mu.1,fproces$mu.2)
#' plot(res$X, col=res$colors)
#' lines(res$X.hat, lwd=2, col=res$color)
#' 
#' fproces <- rproc2mu( 101, "poly1", "poly1", 10, -10, 1, 1)
#' res <- rproc2clust (10,50,fproces$mu.1,fproces$mu.2,"vexponential","vexponential")
#' plot(res$X, col=res$colors)
#' lines(res$X.hat, lwd=2, col=res$color)
#' 
#' fproces <- rproc2mu( 101, "poly1", "poly2", -3, 3, 1, 1.4)
#' res <- rproc2clust (50,50,fproces$mu.1,fproces$mu.2,"brownian","brownian")
#' plot(res$X, col=res$colors)
#' lines(res$X.hat, lwd=2, col=res$color)
#' 
#' par.list.1 <- list(scale = 0.8, theta = 0.3, H = 0.5)
#' par.list.2 <- list(scale = 0.2, theta = 0.3, H = 0.5)
#' fproces <- rproc2mu( 101, "poly1", "poly2", -3, 3, 1, 1.4)
#' res <- rproc2clust (50,50,fproces$mu.1,fproces$mu.2,
#'                     "brownian","brownian", par.list.1, par.list.2)
#' plot(res$X, col=res$colors)
#' lines(res$X.hat, lwd=2, col=res$color)


#' @export
rproc2clust <- function (n.1 = 51, n.2 = n.1, mu.1, mu.2, s.1 = 1, s.2 = 1,
                     par.list.1 = FALSE, par.list.2 = FALSE, ...) {
  t <- mu.1$argvals
  if (is.fdata(mu.1)) mu.1 <- mu.1$data
  if (is.fdata(mu.2)) mu.2 <- mu.2$data
  
  groups <- factor(rep(c(1, 2), c(n.1, n.2)))
  color <- c("darkblue", "darkmagenta")
  colors <- rep(c("blue", "magenta"), c(n.1, n.2))
  
  # Simulación
  if (isFALSE(par.list.1)) {
    X <- rproc2fdata(n.1, t, mu = mu.1, sigma = s.1, ...)
  } else {
    X <- rproc2fdata(n.1, t, mu = mu.1, sigma = s.1, par.list = par.list.1, ...)
  }
  
  if (isFALSE(par.list.2)) {
    X <- c(X, rproc2fdata(n.2, t, mu = mu.2, sigma = s.2, ...))
  } else {
    X <- c(X, rproc2fdata(n.2, t, mu = mu.2, sigma = s.2, par.list = par.list.2, ...))
  }
  
  # Medias funcionales
  X.hat <- func.mean(X[groups == 1])
  X.hat <- c(X.hat, func.mean(X[groups == 2]))
  
  # Salidas
  return (list("X" = X, "X.hat" = X.hat, "groups" = groups,
               "color" = color, "colors" = colors))
}


# ##########################################
# ## Grid of function arguments values
# c.1 <- 25
# c.2 <- rep(c(20, 25, 30, 35), times = 2) * rep(c(-1, 1), each = 4)
# n.1 <- 100
# n.2 <- c(50, 100)
# k.1 <- c(1.0, 1.1)
# k.2 <- 1.0 + 1:4/10
# s.1 <- "vexponential"
# s.2 <- "vexponential"
# par.list.1 <- list(scale = 0.50, theta = 0.3, H = 0.5)
# par.list.2 <- list(scale = 0.25, theta = 0.3, H = 0.5)
# 
# par.args <- expand.grid(s.2, s.1, n.2, n.1, k.2, k.1, c.2, c.1,
#                         stringsAsFactors = FALSE)
# colnames(par.args) <- c("s.2", "s.1", "n.2", "n.1", "k.2", "k.1", "c.2", "c.1")
# 
# # Number of graphs
# n.args <- nrow(par.args) # Number of combinations
# n.graphs <- 16 # Number of graphs per grid
# n.graph <- n.args/n.graphs # Number of images
# 
# t <- seq(0, 1, length.out = 100)
# 
# # Var = 0.5 ####################################################################
# opar <- par(no.readonly = TRUE)
# for (i in 1:n.graph) {
#   # Control sequence
#   seq.g <- 1:n.graphs + n.graphs * (i - 1)
#   for (g in seq.g) {
#     # Simulations
#     set.seed(45)
#     X.sim <- rproc2clust(t,par.args[g, ], par.list.1, par.list.2)
#     ylim <- range(X.sim$X$data) + c(-1.5, 0.0)
#     plot(X.sim$X, main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n",
#          type = "l", ylim = ylim, col = alpha(X.sim$colors, 0.4), lty = 2,
#          lwd = 0.5)
#     lines(X.sim$X.hat, type = "l", col = X.sim$color, lty = 6, lwd = 1)
#     legend('bottom',legend = X.sim$group, col = X.sim$color, pch = 19, lwd = 1,
#            xpd = TRUE, horiz = TRUE, cex = 0.65, seg.len = 1, bty = "n")
#   }
# }
# par(opar)
# 
# 
# # Default ######################################################################
# opar <- par(no.readonly = TRUE)
# 
# for (i in 1:n.graph) {
#   # Secuencia de control
#   seq.g <- 1:n.graphs + n.graphs * (i - 1)
#   par(mfrow = c(4, 4), oma = rep(0.05, 4), mar = rep(0.05, 4))
#   for (g in seq.g) {
#     set.seed(45)
#     X.sim <- rproc2clust(t,par.args[g, ])
#     ylim <- range(X.sim$X$data) + c(-1.5, 0.0)
#     plot(X.sim$X, main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n",
#          type = "l", ylim = ylim, col = alpha(X.sim$colors, 0.4), lty = 2,
#          lwd = 0.5)
#     lines(X.sim$X.hat, type = "l", col = X.sim$color, lty = 6, lwd = 1)
#     legend('bottom',legend = X.sim$group, col = X.sim$color, pch = 19, lwd = 1,
#            xpd = TRUE, horiz = TRUE, cex = 0.65, seg.len = 1, bty = "n")
#   }
# }
# par(opar)

