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


#' @usage locusTag(x)
#' @rdname qualif
#' @export
locusTag <- Curry("qualif", which="locus_tag")


#' @usage product(x)
#' @rdname qualif
#' @export
product <- Curry("qualif", which="product")


#' @usage note(x)
#' @rdname qualif
#' @export
note <- Curry("qualif", which="note")


#' @usage proteinId(x)
#' @rdname qualif
#' @export
proteinId <- Curry("qualif", which="protein_id")


#' @usage translation(x)
#' @rdname qualif
#' @export
translation <- AAStringSet %.% Curry("qualif", which="translation")

