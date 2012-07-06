
# gbFeature-class -----------------------------------------------------

##' @include gbRange-class.r
NULL

setClassUnion("charOrNull", c("character", "NULL"))

##' gbFeature class
##' 
##' \sQuote{gbFeature} is an S4 class that extends the
##' \code{\link{gbFeatureList-class}}. This class provides a container
##' for feature data retrived from GenBank flat files.
##' 
##' \code{gbFeature} provide the following slots:
##' 
##' \describe{
##'    \item{.Dir}{The path to the database file containing the GenBank
##'    record the feature is part of.}
##'    \item{.ACCN}{Accession number of the GenBank record that the
##'    feature is part of.}
##'    \item{.DEF}{The definition line (brief description of the sequence)
##'    of the GenBank record the feature is part of.}
##'    \item{.ID}{Identifier (sequential index) of the feature in the
##'    GenBank record the feature is part of.}
##'    \item{key}{Feature key (e.g. Source, CDS, gene, etc.)}
##'    \item{location}{An object of \code{\link{gbLocation-class}}}
##'    \item{qualifiers}{Named character vector. Name attributes
##'    correspond to GenBank qualifier tags.}       
##' }
##' 
##' @param ... Slots of gbFeature
##'
##' @name gbFeature-class
##' @rdname gbFeature-class
##' @exportClass gbFeature
##' @aliases show,gbFeature-method
##' @aliases start,gbFeature-method
##' @aliases end,gbFeature-method
##' @aliases strand,gbFeature-method
##' @aliases width,gbFeature-method
##' @aliases range,gbFeature-method
##' @aliases partial,gbFeature-method
##' @aliases getIndex,gbFeature-method
##' @aliases getKey,gbFeature-method
##' @aliases getLocation,gbFeature-method
##' @aliases getQualifier,gbFeature-method
##' @aliases dbXref,gbFeature-method
##' @aliases getSequence,gbFeature-method
##' @aliases hasKey,gbFeature-method
##' @aliases hasQualifier,gbFeature-method
##' @aliases [[,gbFeature-method
##' @aliases $,gbFeature-method
.gbFeature <- 
  setClass("gbFeature",
           representation(.Dir="character",
                          .ACCN="character",
                          .DEF="character",
                          .ID="integer",
                          key="character",
                          location="gbLocation",
                          qualifiers="charOrNull"))


# show-method ---------------------------------------------------------


##' @export
setMethod("show", "gbFeature",
          function (object) {
            op <- options("useFancyQuotes")
            options(useFancyQuotes=FALSE)
            indent <- 16
            pad <- blanks(indent)
            len_feat <- nchar(object@key)
            pad_feat <- blanks(indent - len_feat)
            
            # if necessary wrap the lines
            qua <- names(object@qualifiers)
            loc <- linebreak(as(object@location, "character"),
                             offset=indent+1, indent=0, split=",", FORCE=TRUE)
            val <- linebreak(dQuote(object@qualifiers), offset=indent+1, 
                             indent=-(nchar(qua) + 2), FORCE=TRUE)
            
            cat("Feature:         Location/Qualifiers:\n",
                sprintf("%s%s%s\n", object@key, pad_feat, loc),
                sprintf("%s/%s=%s\n", pad, qua, val))
            
            options(op)
            invisible(object)
          })


# Constructor ---------------------------------------------------------


gbFeature <- function (db_dir, accession, definition, id, key, location, qualifiers) 
{
  .gbFeature(.Dir=as.character(db_dir),
             .ACCN=as.character(accession),
             .DEF=as.character(definition),
             .ID=as.integer(id),
             key=as.character(key),
             location=.getLocation(location),
             qualifiers=qualifiers)
}


# Getter-methods ---------------------------------------------------------


##' @export
setMethod("start", "gbFeature",
          function (x, join = FALSE, drop = TRUE) 
            start(x@location, join = join, drop = drop)
          )

##' @export
setMethod("end", "gbFeature",
          function (x, join = FALSE, drop = TRUE) 
            end(x@location, join = join, drop = drop)
          )

##' @export
setMethod("strand", "gbFeature",
          function (x, join = FALSE)
            strand(x@location, join = join)
          )

##' @export
setMethod("width", "gbFeature",
          function (x, join = FALSE)
            width(x@location, join = join)
          )

##' @export
setMethod("partial", "gbFeature",
          function (x)
            partial(x@location)
          )

##' @export
setMethod("range", "gbFeature",
          function (x, join = FALSE)
            range(x@location, join = join)
          )

##' @export
setMethod("getLocation", "gbFeature",
          function (x, attributes = FALSE, join = FALSE) {     
            ans <- range(x@location, join = join)
            ans@elementMetadata$key <- x@key
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

##' @export
setMethod("getIndex", "gbFeature",
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

##' @export
setMethod("getKey", "gbFeature", 
          function (x, attributes=FALSE) {
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


.qualAccess <- function (x, which = "", fixed = FALSE)
{
  .x <- x@qualifiers
  if (fixed) which <- paste0("\\b", which, "\\b") 
  n_row <- length(.x)
  idx <- matrix(
    vapply(which, grepl, names(.x),
           USE.NAMES=FALSE, FUN.VALUE=logical(n_row)),
    nrow=n_row)
  if (ncol(idx) == 1L) {
    ans <- .x[idx]
    if (length(ans) > 0L) {
      return(ans)
    } else {
      return(structure(NA_character_,
                       names=gsub("\\b", "", which, fixed=TRUE)))
    }
  } 
  else if (ncol(idx) > 1L) {
    ans <- lapply(seq.int(ncol(idx)), function (i) .x[idx[,i]])   
    if (any(na_idx <- vapply(ans, length, FUN.VALUE=integer(1L)) == 0L)) {
      for (na in which(na_idx)) {
        ans[[na]] <- structure(NA_character_,
                               names=gsub("\\b", "", which, fixed=TRUE)[na])
      }
    }
    return(unlist(ans))
  }
}


##' @export
setMethod("dbXref", "gbFeature",
          function (x, db = NULL, ...) {     
            ans <- .qualAccess(x, "db_xref")
            if (all(is.na(ans))) {
              return( NA_character_ )
            } else {
              dbs <- vapply(ans, function (x) strsplit(x, ":")[[1]][1], 
                            character(1), USE.NAMES=FALSE)
              ids <- vapply(ans, function (x) strsplit(x, ":")[[1]][2], 
                            character(1), USE.NAMES=FALSE)
              if (is.null(db)) {
                names(ids) <- dbs
                return( ids )
              } else {
                db_pattern <- paste(sprintf("\\b%s\\b", db), collapse="|")
                db_pos <- grep(db_pattern, dbs, ignore.case=TRUE)
                if (length(db_pos) == 0L) {
                  return( NA_character_ )
                } else {
                  ids <- ids[db_pos]
                  names(ids) <- dbs[db_pos]
                  return( ids )
                }
              }
            }
          })


##' @export
setMethod("getQualifier", "gbFeature", 
          function (x, which = "", attributes = FALSE, fixed = FALSE) {
            if (!any(nzchar(which))) {
              ans <- x@qualifiers
            } else {
              ans <- .qualAccess(x, which, fixed)
            }
            if (attributes) {
              structure(ans, key=x@key, id=x@.ID,
                        accession=x@.ACCN,
                        definition=x@.DEF,
                        database=x@.Dir)
            } else {
              ans
            }
          })


.seqAccess <- function(s, x, type)
{
  ## merge Sequences
  mergeSeq <- function (s, x, type) {
    if (length(start(x)) == 1L) {
      seq <- subseq(s, start=start(x), end=end(x))
    } else {
      seq <- do.call(xscat, Map(subseq, s, start=start(x), end=end(x)))
    }
    seq <- switch(type,
                  DNA=DNAStringSet(seq),
                  AA=AAStringSet(seq),
                  RNA=RNAStringSet(seq))
    seq@ranges@NAMES <- sprintf("%s.%s.%s", x@.ACCN, x@key, x@.ID)
    seq
  }
  
  if (is(x, "gbFeatureList")) {
    ## initiate empty XStringSet
    seq <- switch(type, 
                  DNA=DNAStringSet(),
                  AA=AAStringSet(),
                  RNA=RNAStringSet())
    for (i in seq_along(x)) {
      seq[i] <- mergeSeq(s, x[[i]], type)             
    }
  }
  else if (is(x, "gbFeature")) {
    seq <- mergeSeq(s, x, type)
  }
  
  seq@metadata <- list(definition=x@.DEF, database=x@.Dir)
  seq
}

##' @export
setMethod("getSequence", "gbFeature",
          function (x) {
            stopifnot(hasValidDb(x))
            db <- initGB(x@.Dir, verbose=FALSE)
            ans <- .seqAccess(s=dbFetch(db, "sequence"),
                              x, type=dbFetch(db, "type"))
            ans
          })

##' @export
setMethod("hasKey", "gbFeature", 
          function (x, key) 
            !is.na(charmatch(key, x@key))
          )

##' @export
setMethod("hasQualifier", "gbFeature",
          function (x, qualifier)
            !is.na(charmatch(qualifier, names(x@qualifiers)))
          )


# Replacement methods -------------------------------------------------


##' @export
setMethod("start<-", "gbFeature",
          function(x, value) {
            start(x@location) <- value
            x })

##' @export
setMethod("end<-", "gbFeature",
          function(x, value) {
            end(x@location)  <- value
            x })

##' @export
setMethod("strand<-", "gbFeature",
          function(x, value) { 
            strand(x@location) <- value 
            x})


# Shift ---------------------------------------------------------------


setMethod("shift", "gbFeature",
          function(x, shift=0L, ...) {
            x@location <- shift(x@location, shift)
            x
          })


# Subsetting ----------------------------------------------------------


##' @export
setMethod("[[", c("gbFeature", "character", "missing"),
          function(x, i, j) slot(object, i)
          )

##' @export
setMethod("$", "gbFeature",
          function(x, name) slot(x, name)
          )

