#### gbFeatureList objects
setOldClass("list")

##' A list of GenBank features
##' 
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
##' 
##' @examples
##' getSlots("gbFeatureList")
setClass("gbFeatureList", 
         representation(.Dir="character",
                        .ACCN="character",
                        .DEF="character"),
         contains="list")

#### Constructor
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
  
  fl <- new("gbFeatureList", .Data=features, .Dir=as.character(db_dir),
            .ACCN=as.character(accession), .DEF=as.character(definition))
  fl
}

### Accessors ##############################################################

.simplify <- function (x) {
  if (length(len <- unique(unlist(lapply(x, length)))) > 1L)
    return(x)
  if (len == 1L)
    unlist(x, recursive=FALSE)
  else if (len > 1L) {
    n <- length(x)
    r <- as.vector(unlist(x, recursive=FALSE))
    if (prod(d <- c(len, n)) == length(r))
      data.frame(stringsAsFactors=FALSE, matrix(r, nrow=n, byrow=TRUE,
        dimnames=if (!(is.null(nm <- names(x[[1L]])))) list(NULL, nm)))
    else x
  }
  else x
}

##' @aliases index,gbFeatureList,gbFeatureList-method
##' @rdname gbFeatureList-class
setMethod("getIndex",
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

### getKey #################################################################
##' @aliases getKey,gbFeatureList,gbFeatureList-method
##' @rdname accessor-methods
setMethod("getKey",
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


### getLocation ############################################################
##' @aliases getLocation,gbFeatureList,gbFeatureList-method
##' @rdname accessor-methods
setMethod("getLocation", 
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

### getQualifier ###########################################################
##' @aliases getQualifier,gbFeatureList,gbFeatureList-method
##' @rdname accessor-methods
setMethod("getQualifier",
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

### getSequence ############################################################
##' @aliases getSequence,gbFeatureList,gbFeatureList-method
##' @rdname accessor-methods
##' @aliases getSequence,gbFeature,gbFeature-method
##' @rdname accessor-methods
setMethod("getSequence", "gbFeatureList",
          function (x) {
            stopifnot(hasValidDb(x))
            db <- initGenBank(x@.Dir, verbose=FALSE)
            .seqAccess(dbFetch(db, "sequence"), x, dbFetch(db, "type"))
          })
          

##' @aliases hasKey,gbFeatureList,gbFeatureList-method
##' @rdname accessor-methods
setMethod("hasKey", "gbFeatureList", 
          function (x, key) {
            vapply(x, hasKey, key, FUN.VALUE=logical(1L))
          })


##' @aliases hasQualifier,gbFeatureList,gbFeatureList-method
##' @rdname accessor-methods
setMethod("hasQualifier", "gbFeatureList", 
          function (x, qualifier) {
            vapply(x, hasQualifier, qualifier, FUN.VALUE=logical(1L))
          })


### show method ############################################################
##' @export
##' @aliases show,gbFeatureList,gbFeatureList-method
##' @rdname gbFeatureList-class
setMethod("show", signature="gbFeatureList", 
          function (object)  {
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


### Subsetting
##' @export
##' @aliases [,gbFeatureList,gbFeatureList-method
##' @rdname extract-methods
setMethod("[", signature(x="gbFeatureList", i="character", j="missing", drop="missing"),
          function(x, i) {
            idx <- vapply(x@.Data, function(f) f@key, character(1L)) == i
            gbFeatureList(x@.Dir, x@.ACCN, x@.DEF, x@.Data[idx])
          })

##' @export
##' @aliases [,gbFeatureList,gbFeatureList-method
##' @rdname extract-methods
setMethod("[", signature(x="gbFeatureList", i="numeric", j="missing", drop="missing"),
          function(x, i) {
            gbFeatureList(x@.Dir, x@.ACCN, x@.DEF, x@.Data[i])
          })

##' @export
##' @aliases [,gbFeatureList,gbFeatureList-method
##' @rdname extract-methods
setMethod("[", signature(x="gbFeatureList", i="logical", j="missing", drop="missing"),
          function(x, i) {
            gbFeatureList(x@.Dir, x@.ACCN, x@.DEF, x@.Data[i])
          })

##' @export
##' @aliases [,gbFeatureList,gbFeatureList-method
##' @rdname extract-methods
setMethod("[", signature(x="gbFeatureList", i="missing", j="missing", drop="missing"),
          function(x, i) {
            return(x)
          })

### Selecting ##############################################################
##' Select method
##' @export
##' @docType methods
##' @rdname select-methods
setGeneric("select",
           function(x, keys=c(""), cols=c(""))
             standardGeneric("select")
           )

### select method ##########################################################
##' @export
##' @aliases select,gbFeatureList-method
##' @rdname select-methods
setMethod("select",
          signature(x="gbFeatureList"), 
          definition=function (x, keys=c(""), cols=c("")) {
            ans <- .select(x=x, which=keys)
            ans <- .retrieve(x=ans, what=cols)
            ans
          })

    
### viewing ################################################################
##' View method
##' @export
##' @docType methods
##' @rdname view-methods
setGeneric("view",
           function(x, n, ...)
             standardGeneric("view")
           )

### view method ############################################################
##' @export
##' @aliases view,gbFeatureList,gbFeatureList-method
##' @rdname gbFeatureList-class
setMethod("view",
          signature(x="gbFeatureList"), 
          function (x, n)  {
            for (i in x[seq(if (missing(n)) length(x) else n)]){
              show(i)
              cat("\n")
            }
          })
