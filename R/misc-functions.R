#' @importFrom rmisc is.empty
NULL

#' Quickly list all qualifier names
#' 
#' @usage listUniqueQualifs(x)
#' 
#' @param x A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeatureList}},
#' or, \code{\linkS4class{gbFeature}} instance
#' @return A character vector of qualifier names
#' @importFrom rmisc Compose
#' @export
listUniqueQualifs <- Compose(unique, unlist, listQualif)


#' @usage locusTag(x)
#' @rdname qualif
#' @export
locusTag <- Partial(qualif, which="locus_tag", use.names=FALSE)


#' @usage product(x)
#' @rdname qualif
#' @export
product <- Partial(qualif, which="product", use.names=FALSE)


#' @usage note(x)
#' @rdname qualif
#' @export
note <- Partial(qualif, which="note", use.names=FALSE)


#' @usage proteinID(x)
#' @rdname qualif
#' @export
proteinID <- Partial(qualif, which="protein_id", use.names=FALSE)


#' @usage translation(x)
#' @rdname qualif
#' @importFrom Biostrings AAStringSet
#' @export
translation <- Compose(AAStringSet, Partial(qualif, which="translation", use.names=FALSE))

