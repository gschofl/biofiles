#' @include gbRecord-class.R
NULL

#' gbRecordList
#' 
#' \dQuote{gbRecordList} is an S4 class that provides a container for
#' \dQuote{\linkS4class{gbRecord}}s retrived from GenBank flat files.
#'
#' @rdname gbRecordList
#' @export
#' @classHierarchy
#' @classMethods
setClass("gbRecordList", contains="list")


#' @param ... \dQuote{\linkS4class{gbRecord}} elements.
#' @export
gbRecordList <- function (...) {
  listData <- list(...)
  if (length(listData) == 0L) {
    return( new('gbRecordList', .Data = list(new("gbRecord"))) )
  } else {
    if (length(listData) == 1L && is.list(listData[[1L]])) 
      listData <- listData[[1L]]
    if (!all(vapply(listData, is, "gbRecord", FUN.VALUE=logical(1)))) 
      stop("All elements in '...' must be gbRecord objects")
    names(listData) <- sapply(listData, getAccession)
    return( new('gbRecordList', .Data = listData) )
  }
}


setValidity2("gbRecordList", function (object) {
  if (!all(vapply(object@.Data, is, "gbRecord", FUN.VALUE=logical(1))))
    return("All elements in a gbRecordList must be gbRecord objects")
  
  TRUE
})


#' @autoImports
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
  vapply(x, getLocus, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getLength", "gbRecordList", function (x) {
  vapply(x, getLength, FUN.VALUE=integer(1), USE.NAMES=FALSE)
})

setMethod("getMoltype", "gbRecordList", function (x) {
  vapply(x, getMoltype, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getTopology", "gbRecordList", function (x) {
  vapply(x, getTopology, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getDivision", "gbRecordList", function (x) {
  vapply(x, getDivision, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getDate", "gbRecordList", function (x) {
  Map(getDate, x)
})

setMethod("getDefinition", "gbRecordList", function (x) {
  vapply(x, getDefinition, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getAccession", "gbRecordList", function (x) {
  vapply(x, getAccession, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getVersion", "gbRecordList", function (x) {
  vapply(x, getVersion, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getGeneID", "gbRecordList", function (x, db='gi') {
  vapply(x, getGeneID, db=db, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getDBLink", "gbRecordList", function (x) {
  vapply(x, getDBLink, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getDBSource", "gbRecordList", function (x) {
  vapply(x, getDBSource, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getSource", "gbRecordList", function (x) {
  vapply(x, getSource, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getOrganism", "gbRecordList", function (x) {
  vapply(x, getOrganism, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getTaxonomy", "gbRecordList", function (x) {
  vapply(x, getTaxonomy, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getReference", "gbRecordList", function (x) {
  Map(getReference, x)
})

setMethod("getKeywords", "gbRecordList", function (x) {
  vapply(x, getKeywords, FUN.VALUE=character(1), USE.NAMES=FALSE)
})

setMethod("getComment", "gbRecordList", function (x) {
  vapply(x, getComment, FUN.VALUE=character(1), USE.NAMES=FALSE)
})


setMethod("getFeatures", "gbRecordList", 
          function (x) Map(getFeatures, x))


setMethod("getSequence", "gbRecordList", 
          function (x) {
            Reduce(append, Map(getSequence, x))
            })


setMethod("ranges", "gbRecordList",
          function (x, join = FALSE, key = TRUE, include = "none", exclude = "") {
            GRangesList(lapply(x, ranges, join = join, key = key,
                               include = include, exclude = exclude))
          })


setMethod("start", "gbRecordList",
          function (x, join = FALSE, drop = TRUE) {
            lapply(x, start, join = join, drop = drop)
          })


setMethod("end", "gbRecordList",
          function (x, join = FALSE, drop = TRUE) {
            lapply(x, end, join = join, drop = drop)
          })


setMethod("strand", "gbRecordList",
          function (x, join = FALSE) {
            lapply(x, strand, join = join)
          })


setMethod("width", "gbRecordList",
          function (x, join = FALSE) {
            lapply(x, width, join = join)
          })


# internal ---------------------------------------------------------------

setMethod('.dbSource', 'gbRecordList', function (x) {
  vapply(x, .dbSource, character(1), USE.NAMES=FALSE)
})

setMethod(".defline", "gbRecordList", function (x) {
  paste0('gi|', getGeneID(x), .dbSource(x), getAccession(x), ' ', getDefinition(x))
})
