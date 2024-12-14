#' @title Growth Data for Berkeley Growth Study
#'
#' @description 
#' The \code{growth_ldata} dataset contains longitudinal growth data from the Berkeley Growth Study. 
#' The data tracks the growth of 39 boys and 54 girls, measuring height at 31 different ages.
#' 
#' @details 
#' The growth data includes information for boys and girls, with the following structure:
#' 
#' @format 
#' A list containing the following components:
#' \describe{
#'   \item{\code{df}}{\code{data.frame} with the following variables: 
#'     \describe{
#'       \item{\code{group}}{A factor variable with 2 levels indicating the group (boys or girls) in the Berkeley Growth Study.}
#'     }
#'   }
#'   \item{\code{x}}{\code{fdata} class object with \code{n = 93} curves (one curve per row) and \code{31} discretization points (one per column). 
#'   It includes 73 curves corresponding to the heights of 39 boys and 54 girls measured in centimeters at 31 different ages.}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name growth_ldata
#' @usage data(growth_ldata)
#' @source \url{https://example-dataset-source.org}
NULL
