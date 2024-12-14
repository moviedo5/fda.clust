# .onUnload <- function(libpath)
# {
#     library.dynam.unload("fda.clust", libpath)
# }
# 
# ## no S4 methodology here; speedup :
# .noGenerics <- TRUE

#' @importFrom graphics lines points abline
#' @importFrom stats hclust as.dist quantile
#' @import fda.usc
#' @importFrom fda.usc is.fdata fdata metric.lp func.mean func.trim.FM

# Declare global variables to avoid warnings about undefined variables
utils::globalVariables(c(
  "metric.lp", "func.mean", "func.trim.FM", 
  "kmeans.center.ini", "Ker.norm", "gridfdata", 
  "h.default", "rproc2fdata"
))

# Optional .onUnload function if you want to unload the package
# .onUnload <- function(libpath) {
#     library.dynam.unload("fda.clust", libpath)
# }

# Optional: Disable S4 method generation to speed up the package
# .noGenerics <- TRUE
