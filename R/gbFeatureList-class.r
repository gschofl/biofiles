#### gbFeatureList-class ####

setOldClass("list")

##' gbFeatureList class
##' 
##' gbFeatureList is an S4 class that extends the
##' \code{\link{gbRecord-class}}. This class provides a container for
##' feature data retrived from GenBank flat files.
##'
##' gbFeatureList objects have four slots in addition to the slots provided
##' by \code{\link{gbRecord-class}} objects:
##' \describe{
##'    \item{.Dir}{The path to the database file containing the GenBank
##'    record the feature list is part of.}
##'    \item{.ACCN}{Accession number of the GenBank record that the
##'    feature list is part of.}
##'    \item{.DEF}{The definition line (brief description of the sequence)
##'    of the GenBank record the feature list is part of.}
##'    \item{.Data}{A list of gbFeature objects}
##' }
##' 
##' @name gbFeatureList-class
##' @rdname gbFeatureList-class
##' @exportClass gbFeatureList
##' @aliases show,gbFeatureList-method
##' @aliases getIndex,gbFeatureList-method
##' @aliases getKey,gbFeatureList-method
##' @aliases getLocation,gbFeatureList-method
##' @aliases getQualifier,gbFeatureList-method
##' @aliases getSequence,gbFeatureList-method
##' @aliases hasKey,gbFeatureList-method
##' @aliases hasQualifier,gbFeatureList-method
##' @aliases [,gbFeatureList-method
##' @aliases [[,gbFeatureList-method
##' @aliases select,select-method,gbFeatureList-method
##' @aliases view,view-method,gbFeatureList-method
.gbFeatureList <- 
  setClass("gbFeatureList", 
           representation(.Dir="character",
                          .ACCN="character",
                          .DEF="character"),
           contains="list")

##' @export
setMethod("show",
          #### show-method ####
          signature="gbFeatureList", 
          function (object) {
            n_f <- length(object)
            cat(sprintf("'%s' with %i features:\n\n", 
                        class(object), n_f))
            if (n_f > 0L) {
              biofiles::show(object[[1L]])
              if (n_f > 1L) {
                cat("\n...\n")
                biofiles::show(object[[n_f]])
              }
            }
            return(invisible(object))
          })


# Constructor ---------------------------------------------------------


#' @keywords internal
gbFeatureList <- function(db_dir, accession, definition, features)
{
  if (!is.list(features))
    stop("'features' must be a list")
  if (!all(vapply(features, is, "gbFeature", FUN.VALUE=TRUE)))
    stop("all elements in 'features' must be gbFeature objects")
  if (!all(vapply(features, function(f) f@.Dir, FUN.VALUE="") == db_dir))
    stop("all elements in 'features' must be from a valid gbData object")
  if (!all(vapply(features, function(f) f@.ACCN, FUN.VALUE="") == accession))
    stop("all elements in 'features' must be from the same gbData object")
  
  .gbFeatureList(.Data=features, .Dir=as.character(db_dir),
                 .ACCN=as.character(accession), .DEF=as.character(definition)) 
}


# Accessors -----------------------------------------------------------


#' @keywords internal
.simplify <- function (x)
{
  if (length(len <- unique(unlist(lapply(x, length)))) > 1L) {
    return(x)
  }
  if (len == 1L) {
    unlist(x, recursive=FALSE)
  } else if (len > 1L) {
    n <- length(x)
    r <- as.vector(unlist(x, recursive=FALSE))
    if (prod(d <- c(len, n)) == length(r)) {
      data.frame(stringsAsFactors=FALSE,
                 matrix(r, nrow=n, byrow=TRUE,
                        dimnames=if (!(is.null(nm <- names(x[[1L]]))))
                          list(NULL, nm)))
    } else {
      x
    }
  }
  else {
    x
  }
}

##' @export
setMethod("getIndex",
          #### getIndex-method ####
          signature(x="gbFeatureList"),
          function (x, attributes=FALSE, simplify=TRUE) {
            ans <- lapply(x, getIndex)
            if (simplify)
              ans <- .simplify(ans)
            if (attributes)
              ans <- structure(ans,
                               accession=x@.ACCN,
                               definition=x@.DEF,
                               database=x@.Dir)
            ans
          })

##' @export
setMethod("getKey",
          #### getKey-method ####
          signature(x="gbFeatureList"),
          function (x, attributes=FALSE, simplify=TRUE) {
            ans <- lapply(x, getKey)
            if (simplify)
              ans <- .simplify(ans)
            if (attributes)
              ans <- structure(ans,
                               id=getIndex(x, attributes=FALSE, simplify=TRUE),
                               accession=x@.ACCN,
                               definition=x@.DEF,
                               database=x@.Dir)
            ans
          })

##' @export
setMethod("getLocation", 
          #### getLocation-method ####
          signature(x="gbFeatureList"),
          function (x, which=c("start", "end", "strand"), 
                    attributes=FALSE, simplify=TRUE, check=TRUE) {
            
            if (check && !all(grepl("start|end|strand", which, ignore.case=TRUE)))
              stop("Invalid location identifier. Use 'start', 'end', or 'strand'")
            
            ans <- lapply(x, getLocation, which=which, attributes=FALSE,
                          check=FALSE)
            
            if (simplify) 
              ans <- .simplify(ans)
            
            if (attributes && simplify)
              ans <- structure(
                data.frame(stringsAsFactors=FALSE,
                           getIndex(x), getKey(x), ans),
                names=c("id", "key", which),
                accession=x@.ACCN,
                definition=x@.DEF,
                database=x@.Dir)
            else if (attributes && !simplify)
              ans <- structure(
                ans,
                names=paste0(getKey(x),".",getIndex(x)),
                accession=x@.ACCN,
                definition=x@.DEF,
                database=x@.Dir)
            
            ans
          })

##' @export
setMethod("getQualifier",
          #### getQualifier-method ####
          signature(x="gbFeatureList"),
          function (x, which=c(""), attributes=FALSE,
                    simplify=TRUE, fixed=FALSE) {
            
            ans <- lapply(x, getQualifier, which=which, attributes=FALSE,
                          fixed=fixed)
            
            if (simplify)
              ans <- .simplify(ans)
            
            if (attributes && simplify)
              ans <- structure(
                data.frame(stringsAsFactors=FALSE,
                           getIndex(x), getKey(x), ans),
                names=c("id", "key", which),
                accession=x@.ACCN,
                definition=x@.DEF,
                database=x@.Dir)
            else if (attributes && !simplify)
              ans <- structure(
                ans,
                names=paste0(getKey(x),".",getIndex(x)),
                accession=x@.ACCN,
                definition=x@.DEF,
                database=x@.Dir)
            
            ans
          })

##' @export
setMethod("getSequence",
          #### getSequence-method ####
          signature("gbFeatureList"),
          function (x) {
            stopifnot(hasValidDb(x))
            db <- initGenBank(x@.Dir, verbose=FALSE)
            .seqAccess(dbFetch(db, "sequence"), x, dbFetch(db, "type"))
          })
          
##' @export
setMethod("hasKey",
          #### hasKey-method ####
          signature("gbFeatureList"), 
          function (x, key) {
            vapply(x, hasKey, key, FUN.VALUE=logical(1L))
          })

##' @export
setMethod("hasQualifier",
          #### hasQualifier-method ####
          signature("gbFeatureList"), 
          function (x, qualifier) {
            vapply(x, hasQualifier, qualifier, FUN.VALUE=logical(1L))
          })


# Subsetting ----------------------------------------------------------


##' @export
setMethod("[", 
          #### [-method ####
          signature(x="gbFeatureList", i="character", j="missing", drop="missing"),
          function(x, i) {
            idx <- vapply(x@.Data, function(f) f@key, character(1L)) == i
            gbFeatureList(x@.Dir, x@.ACCN, x@.DEF, x@.Data[idx])
          })

##' @export
setMethod("[", 
          signature(x="gbFeatureList", i="numeric", j="missing", drop="missing"),
          function(x, i) {
            gbFeatureList(x@.Dir, x@.ACCN, x@.DEF, x@.Data[i])
          })

##' @export
setMethod("[",
          signature(x="gbFeatureList", i="logical", j="missing", drop="missing"),
          function(x, i) {
            gbFeatureList(x@.Dir, x@.ACCN, x@.DEF, x@.Data[i])
          })

##' @export
setMethod("[",
          signature(x="gbFeatureList", i="missing", j="missing", drop="missing"),
          function(x, i) {
            return(x)
          })


# Select --------------------------------------------------------------


##' Select method
##' 
##' Select features from a GenBank Record
##' 
##' @usage select(x, keys=c(""), cols=c(""))
##' 
##' @param x A \sQuote{\code{gbRecord}} or \sQuote{\code{gbFeatureList}}
##' object
##' @param keys keys
##' @param cols cols
##' 
##' @export
##' @docType methods
##' @rdname select-method
setGeneric("select",
           function(x, keys=c(""), cols=c("")) {
             standardGeneric("select")
           })

##' @export
setMethod("select",
          #### select-method ####
          signature(x="gbFeatureList"), 
          function (x, keys=c(""), cols=c("")) {
            ans <- .select(x=x, which=keys)
            ans <- .retrieve(x=ans, what=cols)
            ans
          })

    
# View ----------------------------------------------------------------


##' View method
##' 
##' @usage view(x, n, ...)
##' 
##' @param x x
##' @param n n
##' @param ... Additional parameters.
##' 
##' @export
##' @docType methods
##' @rdname view-method
setGeneric("view",
           function(x, n, ...) {
             standardGeneric("view")
           })

##' @export
setMethod("view",
          #### view-method ####
          signature(x="gbFeatureList"), 
          function (x, n)  {
            for (i in x[seq(if (missing(n)) length(x) else n)]){
              show(i)
              cat("\n")
            }
          })
