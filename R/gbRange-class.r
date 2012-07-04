
# gbRange-class ----------------------------------------------------------

##' @include gbLocation-class.r
NULL

##' gbRange class
##' 
##' @exportClass gbRange
##' @name gbRange-class
##' @rdname gbRange-class
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
          function(object)
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

