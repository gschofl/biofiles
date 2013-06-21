#' @importFrom rmisc compose
#' @importFrom Biostrings AAStringSet
NULL

#' Quickly list all qualifier names
#' 
#' @usage listUniqueQualifs(x)
#' 
#' @param x A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeatureList}},
#' or, \code{\linkS4class{gbFeature}} instance
#' @return A character vector of qualifier names
#' @export
listUniqueQualifs <- compose(unique, unlist, listQualif)


#' @usage locusTag(x)
#' @rdname qualif
#' @export
locusTag <- Curry(qualif, which="locus_tag", use.names=FALSE)


#' @usage product(x)
#' @rdname qualif
#' @export
product <- Curry(qualif, which="product", use.names=FALSE)


#' @usage note(x)
#' @rdname qualif
#' @export
note <- Curry(qualif, which="note", use.names=FALSE)


#' @usage proteinID(x)
#' @rdname qualif
#' @export
proteinID <- Curry(qualif, which="protein_id", use.names=FALSE)


#' @usage translation(x)
#' @rdname qualif
#' @export
translation <- compose(AAStringSet, Curry(qualif, which="translation", use.names=FALSE))

