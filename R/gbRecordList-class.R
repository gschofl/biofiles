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
              si <- seqinfo(object)
              acc <- seqnames(si)
              len <- unname(getLength(si))
              type <- base::ifelse(
                vapply(object, slot, name='moltype', FUN.VALUE=character(1)) == 'AA',
                'aa', 'bp'
              )
              def <- 
                ellipsize(obj=unname(genome(si)),
                          width=getOption("width") - base::nchar(len) - base::nchar(type) - 8)
              cat(sprintf("[[%s]]\n  %i %s: %s\n", acc, len, type, def), sep="")
            }
          })


setMethod("summary", "gbRecordList",
          function (object, n=2, ...) {
            x <- lapply(object, summary, n=n, ...=...)
            invisible(NULL)
          })


# getters ----------------------------------------------------------------


#' @autoImports
setMethod("seqinfo", "gbRecordList",
          function (x) {
            suppressWarnings(
              base::Reduce(
                GenomicRanges::merge, base::lapply(x, seqinfo)
              )
            )
          })


setMethod("getLength", "gbRecordList",
          function (x) unname(seqlengths(seqinfo(x))))


setMethod("getAccession", "gbRecordList", 
          function (x) seqnames(seqinfo(x)))


setMethod("getGeneID", "gbRecordList", function (x) {
            vapply(x, getGeneID, FUN.VALUE=character(1), USE.NAMES=FALSE)
          })


setMethod("getDefinition", "gbRecordList", 
          function (x) unname(genome(seqinfo(x))))


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
