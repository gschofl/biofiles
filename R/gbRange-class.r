
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


##' @keywords internal
setMethod("initialize",
          signature(.Object = "gbRange"),
          function (.Object, start, width, strand, ...) 
          {
            if (missing(start) || missing(width) || missing(strand)) {
              stop("Missing arguments")
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


##' @keywords internal
setMethod("show", "gbRange",
          function (object)
          {
            lo <- length(object)
            cat(class(object), " of length ", lo, "\n", sep = "")
            if (lo == 0L) 
              return(NULL)
            else {
              showme <- as.data.frame(cbind(as.data.frame(object), as.data.frame(object@elementMetadata)),
                                      row.names = paste("[", seq_len(lo), "]", sep = ""))
              show(showme)
            }  
          })


##' @keywords internal
setMethod("coerce", signature(from="gbRange", to="data.frame"),
          function(from, to, strict) {
            data.frame(as.data.frame(x),
                       as.data.frame(stringsAsFactors=FALSE, x@elementMetadata)
            )
          })


##' @export
setMethod("$", "gbRange",
          function (x, name) as(x, "data.frame")[[name]]
)

##' @export
setMethod("[", "gbRange",
          function (x, i, j, ..., drop) {
            if (missing(j))
              callNextMethod(x, i, j, ..., drop)
            else
              as(x, "data.frame")[i, j, drop]
          })

##' @export
setMethod("[[", "gbRange",
          function (x, i, j, ...)  {
            as(x, "data.frame")[[i]]
          })