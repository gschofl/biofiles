#' @include gbLocation-class.R
NULL

#' gbFeature-class
#' 
#' \dQuote{gbFeature} is an S4 class that provides a container
#' for GenBank feature tables.
#' 
#' @slot .seqinfo An \code{environment} containing the genome sequence as
#' an \code{\linkS4class{XStringSet}} object and sequence metadata
#' as a \code{\linkS4class{gbHeader}} object.
#' @slot .id Identifier (index) of the feature in the
#' GenBank record the feature is part of.
#' @slot key The feature key.
#' @slot location A \code{\linkS4class{gbLocation}} object.
#' @slot qualifiers A named character vector. Name attributes
#' correspond to GenBank qualifier tags.
#' @seealso
#'    \code{\linkS4class{gbFeatureList}}    
#' 
#' @name gbFeature-class
#' @rdname gbFeature-class
#' @exportClass gbFeature
setClass("gbFeature",
         representation(.seqinfo="seqinfo",
                        .id="integer",
                        key="character",
                        location="gbLocation",
                        qualifiers="character"))


setValidity2("gbFeature", function (object) {
  TRUE
})

# show -------------------------------------------------------------------


.showGbFeature <- function(object, showInfo=TRUE) {
  op <- options("useFancyQuotes")
  options(useFancyQuotes=FALSE)
  loc <- linebreak(as(location(object), "character"),
                   width=getOption("width") - 4,
                   offset=17, indent=0, split=",", FORCE=TRUE)
  if (all_empty(object@qualifiers)) {
    cat("Feature:         Location/Qualifiers:\n",
        sprintf("%-16s%s\n", key(object), loc))
  } else {
    qua <- names(object@qualifiers)
    val <- linebreak(dQuote(object@qualifiers),
                     width=getOption("width") - 4, offset=17, 
                     indent=-(nchar(qua) + 2), FORCE=TRUE)
    
    cat("Feature:         Location/Qualifiers:\n",
        sprintf("%-16s%s\n", key(object), loc),
        sprintf("%+17s%s=%s\n", "/", qua, val))
  }
  if (showInfo) {
    cat("Seqinfo:\n")
    show(object@.seqinfo)
  }
  options(op)
}


setMethod("show", "gbFeature",
          function (object) {
            .showGbFeature(object, showInfo=TRUE)
          })


# summary ----------------------------------------------------------------


setMethod("summary", "gbFeature",
          function (object, ...) {
            idx  <- c("N", index(object))
            key  <- c("Key", key(object))
            loc  <- c("Location", as(location(object), "character"))
            prod <- c("Product", product(object))
            idx_len <- max(nchar(idx))
            key_len <- max(nchar(key))
            loc_len <- max(nchar(loc))
            idx <- pad(idx, idx_len + 2, "right")
            key <- pad(key, key_len + 3, "right")
            loc <- pad(loc, loc_len + 3, "right")
            showme <- ellipsize(sprintf("%s%s%s%s", idx, key, loc, prod),
                                width=getOption("width") - 1)
            cat(showme, sep="\n")
            return(invisible(NULL))
          })


# getters ----------------------------------------------------------------

setMethod(".sequence", "gbFeature", function (x) .sequence(x@.seqinfo) )

setMethod("getLocus", "gbFeature", function (x) getLocus(x@.seqinfo) )

setMethod("getLength", "gbFeature", function (x) getLength(x@.seqinfo) )

setMethod("getMoltype", "gbFeature", function (x) getMoltype(x@.seqinfo) )

setMethod("getTopology", "gbFeature", function (x) getTopology(x@.seqinfo) )

setMethod("getDivision", "gbFeature", function (x) getDivision(x@.seqinfo) )

setMethod("getDate", "gbFeature", function (x) getDate(x@.seqinfo) )

setMethod("getDefinition", "gbFeature", function (x) getDefinition(x@.seqinfo) )

setMethod("getAccession", "gbFeature", function (x) getAccession(x@.seqinfo) )

setMethod("getVersion", "gbFeature", function (x) getVersion(x@.seqinfo) )

setMethod("getGeneID", "gbFeature", function (x, db='gi') getGeneID(x@.seqinfo, db=db) )

setMethod("getDBLink", "gbFeature", function (x) getDBLink(x@.seqinfo) )

setMethod("getDBSource", "gbFeature", function (x) getDBSource(x@.seqinfo) )

setMethod("getSource", "gbFeature", function (x) getSource(x@.seqinfo) )

setMethod("getOrganism", "gbFeature", function (x) getOrganism(x@.seqinfo) )

setMethod("getTaxonomy", "gbFeature", function (x) getTaxonomy(x@.seqinfo) )

setMethod("getReference", "gbFeature", function (x) getReference(x@.seqinfo) )

setMethod("getKeywords", "gbFeature", function (x) getKeywords(x@.seqinfo) )

setMethod("getComment", "gbFeature", function (x) getComment(x@.seqinfo) )

setMethod("start", "gbFeature",
          function (x, join = FALSE, drop = TRUE) 
            start(x@location, join = join, drop = drop))


setMethod("end", "gbFeature",
          function (x, join = FALSE, drop = TRUE) 
            end(x@location, join = join, drop = drop))


setMethod("strand", "gbFeature",
          function (x, join = FALSE)
            strand(x@location, join = join))


setMethod("width", "gbFeature",
          function (x, join = FALSE)
            width(x@location, join = join))


setMethod("fuzzy", "gbFeature",
          function (x)
            fuzzy(x@location))


setMethod("ranges", "gbFeature",
          function (x, include = "none", exclude = "", join = FALSE) {
            .make_GRanges(x, include = include, exclude = exclude, join = join)
          })


setMethod("location", "gbFeature",
          function (x) x@location)


setMethod("index", "gbFeature",
          function (x) x@.id)


setMethod("key", "gbFeature", 
          function (x) structure(x@key, names=NULL) )


setMethod("qualif", "gbFeature", 
          function (x, which, fixed=FALSE, use.names=TRUE) {
            if (missing(which)) {
              x@qualifiers
            } else {
              .qual_access(x, which, fixed, use.names)
            }
          })


setMethod("dbxref", "gbFeature",
          function (x, db = NULL, ...) {     
            ans <- .qual_access(x, "db_xref")
            if (all(is.na(ans))) {
              return( NA_character_ )
            } else {
              dbs <- strsplitN(unname(ans), ":", 1)
              ids <- strsplitN(unname(ans), ":", 2)
              if (is.null(db)) {
                setNames(ids, dbs)
              } else {
                db_pattern <- paste0(wrap(db, "\\b"), collapse="|")
                db_pos <- grep(db_pattern, dbs, ignore.case=TRUE)
                if (all_empty(db_pos)) {
                  return( NA_character_ )
                } else {
                  setNames(ids[db_pos], dbs[db_pos])
                }
              }
            }
          })


setMethod("getSequence", "gbFeature", function (x) .seq_access(x))

setMethod('.dbSource', 'gbFeature', function (x) parse_dbsource(getDBSource(x)) )

setMethod(".defline", "gbFeature", function (x) {
  paste0("lcl|", key(x), '.', index(x), .dbSource(x), getAccession(x), ' ',
         getDefinition(x))
})

# setters ----------------------------------------------------------------


setReplaceMethod("start", "gbFeature",
                 function(x, check=TRUE, value) {
                   start(x@location, check=check) <- value
                   if (check)
                     validObject(x)
                   x
                 })


setReplaceMethod("end", "gbFeature",
                 function(x, check=TRUE, value) {
                   end(x@location, check=check) <- value
                   if (check)
                     validObject(x)
                   x
                 })


setReplaceMethod("strand", "gbFeature",
                 function(x, value) { 
                   strand(x@location) <- value
                   x
                 })


setReplaceMethod("key", "gbFeature",
                 function (x, check=TRUE, value) {
                   x <- initialize(x, key=value)
                   if (check)
                     validObject(x)
                   x
                 })


setReplaceMethod("qualif", "gbFeature",
                 function (x, which, check=TRUE, value) {
                   x@qualifiers[which] <- value
                   if (check)
                     validObject(x)
                   x
                 })


# listers ----------------------------------------------------------------


setMethod("listQualif", "gbFeature", 
          function (x) {
            names(x@qualifiers)
          })


# testers ----------------------------------------------------------------



setMethod("hasKey", "gbFeature", 
          function (x, key) {
            !is.na(charmatch(key, x@key))
          })


setMethod("hasQualif", "gbFeature",
          function (x, qualifier) {
            !is.na(charmatch(qualifier, names(x@qualifiers)))
          })


# shift ---------------------------------------------------------------


setMethod("shift", "gbFeature",
          function(x, shift=0L, ...) {
            x@location <- shift(x@location, shift)
            x
          })


# subsetting ----------------------------------------------------------


#' @export
setMethod("[[", c("gbFeature", "character", "missing"),
          function(x, i, j) {
            if (i %in% c("key","location",".id")) {
              slot(x, i)
            } else {
              x@qualifiers[i]
            }
          })


#' @export
setMethod("$", "gbFeature",
          function(x, name) {
            if (name %in% c("key","location",".id")) {
              slot(x, name)
            } else {
              x@qualifiers[name]
            }
          })

