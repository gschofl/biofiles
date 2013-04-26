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
#' @autoImports
gbRecordList <- function (...) {
  listData <- list(...)
  if (length(listData) == 0L) {
    return( new('gbRecordList', .Data = list(new("gbRecord"))) )
  } else {
    if (length(listData) == 1L && is.list(listData[[1L]])) 
      listData <- listData[[1L]]
    if (!all(vapply(listData, is, "gbRecord", FUN.VALUE=logical(1)))) 
      stop("All elements in '...' must be gbRecord objects")
    names(listData) <- sapply(listData, accession)
    return( new('gbRecordList', .Data = listData) )
  }
}


setValidity("gbRecordList", function (object) {
  if (!all(vapply(object@.Data, is, "gbRecord", FUN.VALUE=logical(1))))
    return("All elements in a gbRecordList must be gbRecord objects")
  
  TRUE
})


#' @autoImports
setMethod("show", "gbRecordList",
          function (object) { 
            if (all(is.na(accession(object)))) {
              cat(sprintf("%s instance with zero records\n", sQuote(class(object))))
            } else {
              cat(sprintf("%s instance with %i records\n", 
                          sQuote(class(object)), length(object)))
              si <- seqinfo(object)
              acc <- seqnames(si)
              len <- unname(seqlengths(si))
              type <- ifelse(vapply(object, slot, name='type',
                                    FUN.VALUE=character(1)) == 'AA',
                             'aa', 'bp')
              def <- 
                ellipsize(obj=unname(genome(si)),
                          width=getOption("width") - nchar(len) - nchar(type) - 8)
              cat(sprintf("[[%s]]\n  %i %s: %s\n", acc, len, type, def))
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
            l <- lapply(x, seqinfo)
            suppressWarnings(Reduce(GenomicRanges::merge, l))
          })


setMethod("seqlengths", "gbRecordList",
          function (x) seqlengths(seqinfo(x)))


setMethod("accession", "gbRecordList", 
          function (x) seqnames(seqinfo(x)))


setMethod("definition", "gbRecordList", 
          function (x) genome(seqinfo(x)))


setMethod("features", "gbRecordList", 
          function (x) Map(features, x))


setMethod("sequence", "gbRecordList", 
          function (x) {
            Reduce(append, Map(sequence, x))
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

