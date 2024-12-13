#' @title Hierarchical Clustering for Functional Data
#' 
#' @aliases fhclust
#' @keywords cluster
#' @rdname fhclust
# @name fhclust
#' 
#' @description 
#' Performs hierarchical clustering on functional data using a specified clustering method. 
#' The distance between functional observations is calculated using distance measures 
#' from the \code{fda.usc} package.
#' 
#' @param fdataobj An object of class \code{fdata} representing the functional data to be clustered.
#' Each row corresponds to a functional observation.
#' 
#' @param method A character string specifying the agglomeration method to be used. 
#' Possible values are \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"}, 
#' \code{"mcquitty"}, \code{"median"}, or \code{"centroid"}. Defaults to \code{"ward.D2"}.
#' 
#' @details 
#' The \code{fhclust} function applies hierarchical clustering to functional data, 
#' using distances calculated via the \code{fda.usc::metric.lp} function. 
#' The method for hierarchical clustering can be any of the agglomeration methods available 
#' in \code{hclust}. 
#' 
#' This function is useful for clustering functional data such as time series, curves, 
#' and other functional representations. The function returns an object of class 
#' \code{hclust}, which can be plotted and interpreted as a dendrogram.
#' 
#' @return 
#' An object of class \code{hclust}, which describes the tree produced by the hierarchical clustering process. 
#' The object has the following components:
#' \item{merge}{A numeric matrix describing the merge history.}
#' \item{height}{The height at which the mergers occurred.}
#' \item{order}{A vector giving the order of objects.}
#' \item{labels}{The labels of the objects being clustered.}
#' \item{call}{The call which produced the result.}
#' \item{method}{The agglomeration method used.}
#' 
#' @seealso 
#' \code{\link[stats]{hclust}} for the base R hierarchical clustering function, 
#' and \code{\link[fda.usc]{metric.lp}} for the distance calculation of functional data.
#' 
#' 
#' @examples 
#' \dontrun{
#' t <- seq(0, 2 * pi, length.out = 101)
#' res <- rprocKclust(t, n = c(30, 50, 40), process = c("cos_sin", "sin", "cos"),
#'                    c = c(-1, 1, 1), k = c(NA, NA, NA), s = c(0.2, 0.3, 0.1))
#' # Perform hierarchical clustering using the default method (ward.D2)
#' result <- fhclust(res$X, method = "ward.D2")
#' # Plot the dendrogram
#' plot(result, main = "Dendrogram of Functional Data (Ward's Method)")
#' 
#' # Cut the dendrogram into  clusters
#' groups <- cutree(result, k = 3)
#' print(table(res$groups,groups))
#' }
#' 
#' @export
fhclust <- function(fdataobj, method = "ward.D2") {
  d <- fda.usc::metric.lp(fdataobj)
  hc <- hclust(as.dist(d), method = method)
  return(hc)
}
