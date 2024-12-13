#' @title Functional Mean Shift Clustering
#' 
#' @rdname fmeanshift
# @name fmeanshift
#' @aliases fmeanshift
#' @keywords cluster
#' 
#' @description This function applies the mean shift clustering algorithm to a functional data
#' object of class \code{fdata}.  It uses a kernel-based approach to iteratively shift points 
#' towards high-density regions.
#' 
#' @param fdataobj An object of class \code{fdata} containing the functional data to be clustered.
#' @param h The bandwidth parameter. If \code{h < 0}, the bandwidth is estimated using 
#' the \code{h.default} function with \code{prob = abs(h)}. Defaults to \code{h = -0.15}.
#' @param metric A function to compute the distance between elements of \code{fdataobj}. 
#' Defaults to \code{metric.lp}.
#' @param par.metric A list of additional parameters to be passed to the \code{metric} 
#' function. By default, \code{list(lp = 2)} to compute the L2 distance.
#' @param derr A convergence tolerance parameter used to determine when the mean shift
#'  has converged. Defaults to \code{derr = 0.1}.
#' 
#' @details 
#' The \code{fmeanshift} algorithm iteratively shifts each observation towards the mode 
#' of its neighborhood, defined using a kernel 
#' with bandwidth \code{h}. The procedure continues until the shift distance is smaller 
#' than a convergence threshold controlled by \code{derr}.
#' 
#' The distance between functional data is computed using the distance function \code{metric}, 
#' which defaults to the L2 distance provided by \code{metric.lp} 
#' from the \code{fda.usc} package. The bandwidth \code{h} controls the size of the neighborhood c
#' onsidered for the shift.
#' 
#' @return 
#' A list with the following components:
#' \item{cluster}{An integer vector indicating the cluster assignment for each 
#' observation in \code{fdataobj}.}
#' \item{centers}{An \code{fdata} object representing the centers of the clusters.}
#' 
#' @seealso 
#' See \link[fda.usc:metric.lp]{metric.lp} for distance calculations.
#' @examples 
#' \dontrun{
#' set.seed(8:1)
#' t <- seq(0, 2 * pi, length.out = 101)
#' res <- rprocKclust(t, n = c(40, 40), process = c("cos_sin", "sin"),
#'                    c = c(-1, 2), k = c(NA, NA), s = c(0.3, 0.3))
#' # Run mean shift clustering with automatic bandwidth selection
#' result <- fmeanshift(res$X)
#' # Display cluster assignments and centers
#' table(result$cluster,res$groups)
#' plot(result$centers)
#' plot(res$X, col = result$cluster, main = "functional meanshift")
#' }
#' @export
fmeanshift <- function (fdataobj, h = -0.15, metric = metric.lp,
                        par.metric = list(lp = 2), derr = 0.1) {
  aaf = formals(metric)[-c(1, 2)]
  if (length(par.metric) > 0) {
    nam = names(par.metric)
    for (i in 1:length(par.metric)) {
      aaf[[nam[i]]] = par.metric[[i]]
    }
  }
  D = do.call(metric, c(alist(fdataobj), aaf))
  diag(D) = NA
  
  # Mean-Shift
  minlim = min(D, na.rm = TRUE)
  if (minlim == 0) {
    minlim = quantile(D, probs = 0.01, na.rm = TRUE)
  }
  cat(paste("Min. Dist. InterPares:", minlim, "\n"))
  if (h < 0) {
    h0 = h.default(fdataobj, prob = abs(h))
    cat(paste0("Bandwidth:", h0, "\n"))
  } else {
    h0 = h
  }
  iclus = rep(0, nrow(fdataobj))
  
  for (i in 1:nrow(fdataobj)) {
    conv = FALSE
    if (is.fdata(fdataobj)) {
      x0 = fdataobj[i]
    } else {
      x0 = fdataobj[i, , drop = FALSE]
    }
    while (!conv) {
      dd = do.call(metric, c(alist(x0, fdataobj), aaf))
      wd = Ker.norm(dd / h0)
      wd = wd/sum(wd)
      if (is.fdata(fdataobj)) {
        fx = gridfdata(coef = matrix(wd, nrow = 1), fdataobj) 
      } else {
        fx = matrix(wd, nrow = 1) %*% fdataobj
      }
      eps = do.call(metric, c(alist(fx, x0), aaf))
      if (eps < minlim * derr * 0.1) {
        conv = TRUE
      } else {
        x0 = fx
      }
    }
    if (i == 1) {
      mshift = fx
      iclus[i] = 1L
    } else {
      dm = do.call(metric, c(alist(fx, mshift), aaf)) 
      if (any(dm < minlim / 2)) {
        iclus[i] = which(dm < minlim / 2)[1]
      } else {
        if (is.fdata(fdataobj)) {
          mshift = c(mshift, fx)
        } else {
          mshift = rbind(mshift, fx)
        }
        iclus[i] = max(iclus) + 1
      }
    }
  }
  return(list(cluster = iclus, centers = mshift))
}



