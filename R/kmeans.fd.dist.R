predict_kmeans.fd <- function(object, newdata, ...){
  x = object$centers
  n = NROW(x)
  np <- NCOL(x)
  nn = NROW(newx)
  pp = NCOL(newx)
  if (pp != np) stop("newdata have wrong dimension")
  if (is.null(rownames(newx))) 
    rownames(newx) <- 1:nn
  par.metric <- attr(object$mdist, "par.metric")
  par.metric[["fdata1"]] <- newdata
  par.metric[["fdata2"]] <- x
  a1 <- attr(object$mdist, "call")
  nmdist <- do.call(a1, par.metric)
  #for (i in 1:nr) {  grupo[i] = which.min(d[(nr + 1):(nr + ngroups), i])}  
  grupo <- apply(nmdist,1,which.min)
  return(grupo)
}

