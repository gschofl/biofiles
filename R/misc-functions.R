#' @importFrom rmisc %.%
NULL

#' Quickly list all qualifier names
#' 
#' @usage listUniqueQualifs(x)
#' 
#' @param x A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeatureList}},
#' or, \code{\linkS4class{gbFeature}} instance
#' @return A character vector of qualifier names
#' @export
listUniqueQualifs <- unique %.% unlist %.% listQualif