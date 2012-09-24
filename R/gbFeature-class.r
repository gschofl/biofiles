
# gbFeature-class -----------------------------------------------------

#' @include gbRange-class.r
NULL

setClassUnion("charOrNull", c("character", "NULL"))

#' gbFeature class
#' 
#' @description
#' \sQuote{\code{gbFeature}} is an S4 class that extends the class
#' \code{\linkS4class{gbFeatureList}}. This class provides a container
#' for feature data retrived from GenBank flat files.
#' 
#' @slot .Dir The path to the database file containing the GenBank
#' record the feature is part of.
#' @slot .ACCN Accession number of the GenBank record that the
#' feature is part of.
#' @slot .DEF The definition line (brief description of the sequence)
#' of the GenBank record the feature is part of.
#' @slot .ID Identifier (sequential index) of the feature in the
#' GenBank record the feature is part of.
#' @slot key
#' @slot location
#' @slot qualifiers Named character vector. Name attributes
#'    correspond to GenBank qualifier tags.       
#' 
#' @rdname gbFeature
#' @exportClass gbFeature
#' @classHierarchy
#' @classMethods
.gbFeature <- setClass("gbFeature",
                       representation(.Dir="character",
                                      .ACCN="character",
                                      .DEF="character",
                                      .ID="integer",
                                      key="character",
                                      location="gbLocation",
                                      qualifiers="charOrNull"))


setValidity("gbFeature", function (object) {
  # at the moment do nothing but the default checks
  TRUE
})


# show -------------------------------------------------------------------


#' @autoImports
setMethod("show", "gbFeature",
          function (object) {
            op <- options("useFancyQuotes")
            options(useFancyQuotes=FALSE)

            loc <- linebreak(as(object@location, "character"),
                             offset=17, indent=0, split=",", FORCE=TRUE)
            
            if (is.null(object@qualifiers)) {
              cat("Feature:         Location/Qualifiers:\n",
                  sprintf("%-16s%s\n", object@key, loc))
            } else {
              qua <- names(object@qualifiers)
              val <- linebreak(dQuote(object@qualifiers), offset=17, 
                               indent=-(nchar(qua) + 2), FORCE=TRUE)
              
              cat("Feature:         Location/Qualifiers:\n",
                  sprintf("%-16s%s\n", object@key, loc),
                  sprintf("%+17s%s=%s\n", "/", qua, val))
            }

            options(op)
            invisible(object)
          })


# summary ----------------------------------------------------------------


setMethod("summary", "gbFeature",
    function (object, ...) {
        showme <- sprintf("%-6s%-20s%-32s\n",
                          object@.ID, object@key, as(object@location, "character"))
        cat(showme)
        return(invisible(TRUE))
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


setMethod("accession", "gbFeature",
          function (x) x@.ACCN)


setMethod("definition", "gbFeature",
          function (x) x@.DEF)


setMethod("range", "gbFeature",
          function (x, join = FALSE)
            range(x@location, join = join))


setMethod("location", "gbFeature",
          function (x, attributes = FALSE, join = FALSE) {     
            ans <- range(x@location, join = join)
            ans@elementMetadata$feature <- x@key
            ans@elementMetadata$id <- x@.ID
            if (attributes) {
              structure(ans,
                        accession=x@.ACCN,
                        definition=x@.DEF,
                        database=x@.Dir)
            } else {
              ans
            }
          })


setMethod("index", "gbFeature",
          function (x, attributes = FALSE) {
            ans <- x@.ID
            if (attributes) {
              structure(ans,
                        accession=x@.ACCN,
                        definition=x@.DEF,
                        database=x@.Dir)
            } else {
              ans
            }
          })


setMethod("key", "gbFeature", 
          function (x, attributes = FALSE) {
            ans <- structure(x@key, names=NULL)
            if (attributes) {
              structure(ans, id=x@.ID,
                        accession=x@.ACCN,
                        definition=x@.DEF,
                        database=x@.Dir)
            }
            else {
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
                if (is_empty(db_pos)) {
                  return( NA_character_ )
                } else {
                  structure(ids[db_pos], names = dbs[db_pos])
                }
              }
            }
          })


setMethod("qualif", "gbFeature", 
          function (x, which, attributes = FALSE, fixed = FALSE) {
            if (missing(which)) {
              ans <- x@qualifiers
            } else {
              ans <- .qualAccess(x, which, fixed)
            }
            if (attributes) {
              structure(ans, id=x@.ID,
                        accession=x@.ACCN,
                        definition=x@.DEF,
                        database=x@.Dir)
            } else {
              ans
            }
          })


setMethod("sequence", "gbFeature",
          function (x) {
            stopifnot(hasValidDb(x))
            db <- initGB(x@.Dir, verbose=FALSE)
            ans <- .seqAccess(dbFetch(db, "sequence"), x, dbFetch(db, "type"))
            ans
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
                     db <- dbInit(x@.Dir, "RDS")
                     db$features[x@.ID] <- x
                   }
                   validObject(x)
                   x
                 })


setReplaceMethod("qualif", "gbFeature",
                 function (x, which, value, updateDb = FALSE) {
                   x@qualifiers[which] <- value
                   if (updateDb) {
                     db <- dbInit(x@.Dir, "RDS")
                     db$features[x@.ID] <- x
                   }
                   validObject(x)
                   x
                 })


# testers ----------------------------------------------------------------



setMethod("hasKey", "gbFeature", 
          function (x, key) 
            not.na(charmatch(key, x@key)))


setMethod("hasQualif", "gbFeature",
          function (x, qualifier)
            not.na(charmatch(qualifier, names(x@qualifiers))))


# shift ---------------------------------------------------------------


setMethod("shift", "gbFeature",
          function(x, shift=0L, ...) {
            x@location <- shift(x@location, shift)
            x
          })


# subsetting ----------------------------------------------------------


#' @export
setMethod("[[", c("gbFeature", "character", "missing"),
          function(x, i, j) slot(x, i))


#' @export
setMethod("$", "gbFeature",
          function(x, name) slot(x, name))

