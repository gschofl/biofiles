
# gbRange-class ----------------------------------------------------------

##' @include gbLocation-class.r
NULL

##' gbRange class
##' 
##' @param ... Slots of gbRange
##' 
##' @exportClass gbRange
##' @name gbRange-class
##' @rdname gbRange-class
##' @aliases $,gbRange-method
##' @aliases [,gbRange-method
##' @aliases [[,gbRange-method
.gbRange <- setClass("gbRange", contains="IRanges")


# initialize-method ------------------------------------------------------


#' @keywords internal
setMethod("initialize",
          signature(.Object = "gbRange"),
          function (.Object, start, width, strand, ...) {
            if (missing(start) || missing(width) || missing(strand)) {
              stop("Missing arguments")
            }
            if (!all(strand %in% c(1,-1))) {
              stop("Strand must be encoded as 1 (plus strand) or -1 (minus strand)")
            }
            anno=list(...)
            z <- vapply(c(list(start, width, strand), anno), length, numeric(1))
            if (length(unique(z)) != 1L) {
              stop("Arguments have unequal length")
            }
            r <- callNextMethod(.Object, start = start, width = width)
            r@elementMetadata <- if (length(anno) > 0) {
              DataFrame(strand, anno)
            } else {
              DataFrame(strand)
            }
            r
          })


# show-method ------------------------------------------------------------


#' @export
setMethod("show", "gbRange",
          function (object) {
            lo <- length(object)
            cat(sprintf("%s of length %s\n", class(object), lo))
            if (lo == 0L){
              return(NULL)
            } 
            if (lo < 20L) {
              showme <- 
                as.data.frame(cbind(as.data.frame(object),
                                    as.data.frame(object@elementMetadata)),
                              row.names = paste0("[", seq_len(lo), "]"))
            } else {
              n <- 8
              headshow <- as(head(object, n), "data.frame")
              tailshow <- as(tail(object, n), "data.frame")
              rows <- c(paste0("[", c(1:n), "]"),
                        "---",
                        paste0("[", (lo - n + 1):lo, "]"))
              showme <- as.data.frame(rbind(headshow, "---", tailshow),
                                      row.names = rows)
            }
            show(showme)
          })


#' As(gbRange, "data.frame")
#'
#' @name as
setAs("gbRange", "data.frame",
      function (from) {
        data.frame(as.data.frame(from),
                   as.data.frame(stringsAsFactors=FALSE,
                                 from@elementMetadata))
      })


# subsetting-methods -----------------------------------------------------


#' @export
setMethod("$", "gbRange",
          function (x, name) as(x, "data.frame")[[name]]
)

## Ignores the 'drop' argument and behaves as if it was set to FALSE
#' @export
setMethod("[", "gbRange",
          function (x, i, j, ..., drop=TRUE) {
            if (missing(j) && length(list(...)) == 0L)
              callNextMethod(x=x, i=i)
            else
              as(x, "data.frame")[i, j, ..., drop=FALSE]
          })


#' @export
setMethod("[[", "gbRange",
          function (x, i, j, ...)  {
            as(x, "data.frame")[[i]]
          })
