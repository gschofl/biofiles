#' @include gbRecord-class.R
NULL

#' Class \code{"gbRecordList"}
#' 
#' \dQuote{gbRecordList} is an S4 class that provides a container for
#' \dQuote{\linkS4class{gbRecord}}s retrived from GenBank flat files.
#' For instantiation of a gbRecordList object use the import function
#' \code{\link{gbRecord}} or combine \code{gbRecord} objects using
#' \code{gbRecordList}.
#'
#' @name gbRecordList-class
#' @rdname gbRecordList-class
#' @exportClass gbRecordList
setClass("gbRecordList", contains="list")

#' @param ... \dQuote{\linkS4class{gbRecord}} elements.
#' @return A \dQuote{\linkS4class{gbRecordList}} instance.
#' @rdname gbRecordList-class
#' @export
#' @examples
#' 
#' ###
#' 
gbRecordList <- function(...) {
  listData <- list(...)
  if (length(listData) == 0L) {
    return( new('gbRecordList', .Data = list(new("gbRecord"))) )
  } else {
    if (length(listData) == 1L && is.list(listData[[1L]])) {
      listData <- listData[[1L]]
    }
    if (any(vapply(listData, is, "gbRecordList", FUN.VALUE=FALSE))) {
      listData <- flatten1(listData)
    }
    if (!all(vapply(listData, is, "gbRecord", FUN.VALUE=FALSE))) {
      stop("All elements in '...' must be gbRecord objects")
    }
    names(listData) <- vapply(listData, getAccession, "", USE.NAMES=FALSE)
    return( new('gbRecordList', .Data = listData) )
  }
}


setValidity2("gbRecordList", function (object) {
  if (!all(vapply(object@.Data, is, "gbRecord", FUN.VALUE=logical(1))))
    return("All elements in a gbRecordList must be gbRecord objects")
  
  TRUE
})


setMethod("show", "gbRecordList",
          function (object) { 
            if (all(is.na(getAccession(object)))) {
              cat(sprintf("%s instance with zero records\n", sQuote(class(object))))
            } else {
              cat(sprintf("%s instance with %i records\n", 
                          sQuote(class(object)), length(object)))
              acc <- getAccession(object)
              len <- getLength(object)
              type <- ifelse(getMoltype(object) == 'AA', 'aa', 'bp')
              def <- ellipsize(obj=getDefinition(object),
                               width=getOption("width") - nchar(len) - nchar(type) - 8)
              cat(sprintf("[[%s]]\n  %i %s: %s\n", acc, len, type, def), sep="")
            }
          })


setMethod("summary", "gbRecordList",
          function (object, n=2, ...) {
            x <- lapply(object, summary, n=n, ...=...)
            invisible(NULL)
          })


# getters ----------------------------------------------------------------


setMethod("getLocus", "gbRecordList", function (x) {
  vapply(x, getLocus, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getLength", "gbRecordList", function (x) {
  vapply(x, getLength, FUN.VALUE=0L, USE.NAMES=FALSE)
})

setMethod("getMoltype", "gbRecordList", function (x) {
  vapply(x, getMoltype, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getTopology", "gbRecordList", function (x) {
  vapply(x, getTopology, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getDivision", "gbRecordList", function (x) {
  vapply(x, getDivision, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getDate", "gbRecordList", function (x) {
  Map(getDate, x)
})

setMethod("getDefinition", "gbRecordList", function (x) {
  vapply(x, getDefinition, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getAccession", "gbRecordList", function (x) {
  vapply(x, getAccession, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getVersion", "gbRecordList", function (x) {
  vapply(x, getVersion, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getGeneID", "gbRecordList", function (x, db='gi') {
  vapply(x, getGeneID, db=db, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getDBLink", "gbRecordList", function (x) {
  vapply(x, getDBLink, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getDBSource", "gbRecordList", function (x) {
  vapply(x, getDBSource, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getSource", "gbRecordList", function (x) {
  vapply(x, getSource, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getOrganism", "gbRecordList", function (x) {
  vapply(x, getOrganism, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getTaxonomy", "gbRecordList", function (x) {
  vapply(x, getTaxonomy, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getReference", "gbRecordList", function (x) {
  Map(getReference, x)
})

setMethod("getKeywords", "gbRecordList", function (x) {
  vapply(x, getKeywords, FUN.VALUE="", USE.NAMES=FALSE)
})

setMethod("getComment", "gbRecordList", function (x) {
  vapply(x, getComment, FUN.VALUE="", USE.NAMES=FALSE)
})

#' @export
#' @aliases getFeatures-method,gbRecordList-method
#' @rdname getFeatures-methods
setMethod("getFeatures", "gbRecordList", function(x) {
  .mapply(.features, list(x = x), NULL)
})

#' @export
#' @aliases getFeatures-method,gbRecordList-method
#' @rdname getFeatures-methods
setMethod("ft", "gbRecordList", function(x) {
  .mapply(.features(), list(x = x), NULL)
})
          
#' @export
#' @aliases getSequence-method,gbRecordList-method
#' @rdname getSequence-methods
setMethod("getSequence", "gbRecordList", function (x) {
  Reduce(append, .mapply(.sequence, list(x = x), NULL))
})

# ' @export
# ' @aliases ranges,gbRecordList-method
# ' @rdname ranges-methods
setMethod("ranges", "gbRecordList",
          function (x, join = FALSE, key = TRUE, include = "none", exclude = "") {
            GRangesList(lapply(x, ranges, join = join, key = key,
                               include = include, exclude = exclude))
          })

#' @export
#' @aliases start,gbRecordList-method
#' @rdname start-methods
setMethod("start", "gbRecordList", function (x, join = FALSE, drop = TRUE) {
  lapply(x, start, join = join, drop = drop)
})

#' @export
#' @aliases end,gbRecordList-method
#' @rdname end-methods
setMethod("end", "gbRecordList", function (x, join = FALSE, drop = TRUE) {
  lapply(x, end, join = join, drop = drop)
})

#' @export
#' @aliases strand,gbRecordList-method
#' @rdname strand-methods
setMethod("strand", "gbRecordList", function (x, join = FALSE) {
  lapply(x, strand, join = join)
})

#' @export
#' @aliases width,gbRecordList-method
#' @rdname width-methods
setMethod("width", "gbRecordList", function(x) {
  lapply(x, width)
})

setMethod("joint_width", "gbRecordList", function(x) {
  lapply(x, joint_width)
})

#' @export
#' @aliases fuzzy,gbRecordList-method
#' @rdname fuzzy-methods
setMethod("fuzzy", "gbRecordList", function(x) {
  lapply(x, fuzzy)
})

#' @export
#' @aliases index,gbRecordList-method
#' @rdname index-methods
setMethod("index", "gbRecordList", function(x) {
  lapply(x, index)
})

#' @export
#' @aliases key,gbRecordList-method
#' @rdname key-methods
setMethod("key", "gbRecordList", function(x) {
  lapply(x, key)
})


# listers -------------------------------------------------------------------


##
## listQualif would return a list of lists for gbRecordList
## users should "lapply" over gbRecordLists
##

#' @export
#' @aliases tableQualif,gbRecordList-method
#' @rdname tableQualif-methods
setMethod("tableQualif", "gbRecordList", function(x) {
  lapply(x, tableQualif)
})


# testers ----------------------------------------------------------------


#' @export
#' @aliases hasKey,gbRecordList-method
#' @rdname hasKey-methods
setMethod("hasKey", "gbRecordList", function(x, key) {
  .mapply(hasKey, list(x = x), list(key = key))
})

#' @export
#' @aliases hasQualif,gbRecordList-method
#' @rdname hasQualif-methods
setMethod("hasQualif", "gbRecordList", function(x, qualifier) {
  .mapply(hasQualif, list(x = x), list(qualifier = qualifier))
})


# internal ---------------------------------------------------------------


setMethod('.dbSource', 'gbRecordList', function (x) {
  vapply(x, .dbSource, character(1), USE.NAMES=FALSE)
})

setMethod(".defline", "gbRecordList", function (x) {
  paste0('gi|', getGeneID(x), .dbSource(x), getAccession(x), ' ', getDefinition(x))
})
