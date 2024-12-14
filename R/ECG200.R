#' @title ECG 200 Dataset
#'
#' @description 
#' Each series traces the electrical activity recorded during one heartbeat. 
#' The two classes are a normal heartbeat and a Myocardial Infarction.
#' 
#' @details 
#' This dataset was formatted by R. Olszewski as part of his thesis 
#' "Generalized feature extraction for structural pattern recognition in time-series data" 
#' at Carnegie Mellon University, 2001.
#' 
#' @format 
#' A list containing the following components:
#' \describe{
#'   \item{\code{df}}{\code{data.frame} with the following variables: 
#'     \describe{
#'       \item{\code{class}}{Corresponding class level of ``CinCECGtorso'' curves with 2 classes.}
#'       \item{\code{sample}}{Factor variable. In the TSC database, the first 100 values (\code{sample="train"}) are used as the training sample, and the remaining 100 (\code{sample="test"}) for testing.}
#'     }
#'   }
#'   \item{\code{x}}{\code{fdata} class object with \code{n = 200} curves (one curve per row) and \code{96} discretization points (one per column).}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name ECG200
#' @usage data(ECG200)
#' @source \url{http://timeseriesclassification.com/description.php?Dataset=CinCECGtorso}
NULL
