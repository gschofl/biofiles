
# gbFeature-class -----------------------------------------------------

#' @include gbInfo-class.R
NULL

setClassUnion("charOrNull", c("character", "NULL"))

#' gbFeature
#' 
#' \dQuote{gbFeature} is an S4 class that provides a container
#' for GenBank feature tables.
#' 
#' @slot .Info A \code{\linkS4class{gbInfo}} object.
#' @slot .Id Identifier (index) of the feature in the
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
         representation(.Info="gbInfo",
                        .Id="integer",
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


setMethod("summary", "gbFeature",
          function (object, ...) {
            idx <- pad(index(object), 8, "right")
            key <- pad(key(object), 10, "right")
            loc <- as(location(object)[[1]], "character")
            prod <- ellipsize(product(object), width=getOption("width") - 
                                nchar(idx) - nchar(key) - nchar(loc) - 5)
            showme <- sprintf("%s%s%s  %s\n", idx, key, loc, prod)
            cat(showme)
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


setMethod("partial", "gbFeature",
          function (x)
            partial(x@location))


setMethod("seqinfo", "gbFeature",
          function (x) x@.Info)


setMethod("seqlengths", "gbFeature",
          function (x) seqlengths(seqinfo(x)))


setMethod("accession", "gbFeature",
          function (x) seqnames(seqinfo(x)))


setMethod("definition", "gbFeature",
          function (x) genome(seqinfo(x)))


setMethod("ranges", "gbFeature",
          function (x, with_qual = "none", without_qual = "", join = FALSE) {
            .make_GRanges(x, with_qual = with_qual,
                          without_qual = without_qual,
                          join = join)
          })


setMethod("location", "gbFeature",
          function (x, seqinfo = FALSE) {     
            ans <- x@location
            if (seqinfo) {
              structure(list(ans),
                        accession=accession(x),
                        definition=unname(definition(x)),
                        dir=seqinfo(x)@db@dir)
            } else {
              ans
            }
          })


setMethod("index", "gbFeature",
          function (x, seqinfo = FALSE) {
            ans <- x@.Id
            if (seqinfo) {
              structure(ans,
                        accession=accession(x),
                        definition=unname(definition(x)),
                        dir=seqinfo(x)@db@dir)
            } else {
              ans
            }
          })


setMethod("key", "gbFeature", 
          function (x, seqinfo = FALSE) {
            ans <- structure(x@key, names=NULL)
            if (seqinfo) {
              structure(ans,
                        accession=accession(x),
                        definition=unname(definition(x)),
                        dir=seqinfo(x)@db@dir)
            }
            else {
              ans
            }
          })


setMethod("qualif", "gbFeature", 
          function (x, which, seqinfo = FALSE, fixed = FALSE) {
            if (missing(which)) {
              ans <- x@qualifiers
            } else {
              ans <- .qualAccess(x, which, fixed)
            }
            if (seqinfo) {
              structure(ans,
                        accession=accession(x),
                        definition=unname(definition(x)),
                        dir=seqinfo(x)@db@dir)
            } else {
              ans
            }
          })


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
          function (x) {
            stopifnot(hasValidDb(x))
            db <- slot(seqinfo(x), "db")
            .seqAccess(s=dbFetch(db, "sequence"), x, type=dbFetch(db, "type"))
          })


# setters ----------------------------------------------------------------


setReplaceMethod("start", "gbFeature",
                 function(x, value) {
                   start(x@location) <- value
                   validObject(x)
                   x 
                 })


setReplaceMethod("end", "gbFeature",
                 function(x, value) {
                   end(x@location) <- value
                   validObject(x)
                   x
                 })


setReplaceMethod("strand", "gbFeature",
                 function(x, value) { 
                   strand(x@location) <- value
                   validObject(x)
                   x
                 })


setReplaceMethod("key", "gbFeature",
                 function (x, value, updateDb = FALSE) {
                   x <- initialize(x, key=value)
                   if (updateDb) {
                     db <- slot(seqinfo(x), "db")
                     db$features[x@.Id] <- x
                   }
                   validObject(x)
                   x
                 })


setReplaceMethod("qualif", "gbFeature",
                 function (x, which, value, updateDb = FALSE) {
                   x@qualifiers[which] <- value
                   if (updateDb) {
                     db <- slot(seqinfo(x), "db")
                     db$features[x@.Id] <- x
                   }
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

