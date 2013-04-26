#' @include gbLocation-class.R
NULL


#' gbFeature
#' 
#' \dQuote{gbFeature} is an S4 class that provides a container
#' for GenBank feature tables.
#' 
#' @slot .seqinfo An \code{environment} containing the genome sequence as
#' an \code{\linkS4class{XStringSet}} object and sequence metadata
#' as a \code{\linkS4class{Seqinfo}} object.
#' @slot .id Identifier (index) of the feature in the
#' GenBank record the feature is part of.
#' @slot key The feature key.
#' @slot location A \code{\linkS4class{gbLocation}} object.
#' @slot qualifiers A named character vector. Name attributes
#' correspond to GenBank qualifier tags.       
#' 
#' @rdname gbFeature
#' @export
#' @classHierarchy
#' @classMethods
setClass("gbFeature",
         representation(.seqinfo="environment",
                        .id="integer",
                        key="character",
                        location="gbLocation",
                        qualifiers="character"),
         prototype(.seqinfo=new.env(parent=emptyenv())))


setValidity2("gbFeature", function (object) {
  TRUE
})

# show -------------------------------------------------------------------


#' @autoImports
.showGbFeature <- function(object, showInfo=TRUE) {
  op <- options("useFancyQuotes")
  options(useFancyQuotes=FALSE)
  loc <- linebreak(as(object@location, "character"),
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
    showInfo(seqinfo(object))
  }
  options(op)
}


#' @autoImports
setMethod("show", "gbFeature",
          function (object) {
            .showGbFeature(object, showInfo=TRUE)
          })


# summary ----------------------------------------------------------------


#' @autoImports
setMethod("summary", "gbFeature",
          function (object, ...) {
            idx <- c("N", index(object))
            key <- c("Key", key(object))
            loc <- c("Location", as(location(object), "character"))
            prod <- c("Product", product(object))
            idx_len <-max(nchar(idx))
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


setMethod("seqinfo", "gbFeature",
          function (x) {
            tryCatch(get("seqinfo", x@.seqinfo),
                     error = function (e) Seqinfo() )
          })


#' @autoImports
setMethod("seqlengths", "gbFeature",
          function (x) seqlengths(seqinfo(x)))


#' @autoImports
setMethod("accession", "gbFeature",
          function (x) seqnames(seqinfo(x)))


#' @autoImports
setMethod("definition", "gbFeature",
          function (x) genome(seqinfo(x)))


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
          function (x, which, fixed = FALSE) {
            if (missing(which)) {
              x@qualifiers
            } else {
              .qualAccess(x, which, fixed)
            }
          })


#' @autoImports
setMethod("dbxref", "gbFeature",
          function (x, db = NULL, ...) {     
            ans <- .qualAccess(x, "db_xref")
            if (all(is.na(ans))) {
              return( NA_character_ )
            } else {
              dbs <- strsplitN(ans, ":", 1)
              ids <- strsplitN(ans, ":", 2)
              if (is.null(db)) {
                structure(ids, names = dbs)
              } else {
                db_pattern <- paste(wrap(db, "\\b"), collapse="|")
                db_pos <- grep(db_pattern, dbs, ignore.case=TRUE)
                if (all_empty(db_pos)) {
                  return( NA_character_ )
                } else {
                  structure(ids[db_pos], names = dbs[db_pos])
                }
              }
            }
          })


setMethod("sequence", "gbFeature",
          function (x) .seqAccess(x))


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
                 function(x, check=TRUE, value) { 
                   strand(x@location, check=check) <- value
                   if (check)
                     validObject(x)
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
            if (i %in% c("key","location", ".Id")) {
              slot(x, i)
            } else {
              x@qualifiers[i]
            }
          })

#' @export
setMethod("$", "gbFeature",
          function(x, name) {
            if (name %in% c("key","location",".Id")) {
              slot(x, name)
            } else {
              x@qualifiers[name]
            }
          })

