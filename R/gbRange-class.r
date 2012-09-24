#' @include gbLocation-class.r
NULL

# gbRange-class ----------------------------------------------------------

#' gbRange class
#' 
#' @rdname gbRange
#' @export
#' @classHierarchy
#' @classMethods
.gbRange <- setClass("gbRange", contains="IRanges")


# initialize-method ------------------------------------------------------


#' @keywords internal
#' @autoImports
setMethod("initialize", "gbRange",
          function (.Object, start, width, strand, ...) {
            if (missing(start) || missing(width) || missing(strand)) {
              stop("Missing arguments")
            }
            if (any(strand %ni% c(1,-1))) {
              stop("Strand must be encoded as 1 (plus strand) or -1 (minus strand)")
            }
            anno=list(...)
            z <- vapply(c(list(start, width, strand), anno), length, numeric(1))
            if (length(unique(z)) != 1L) {
              stop("Arguments have unequal length")
            }
            r <- callNextMethod(.Object, start = start, width = width)
            elementMetadata(r) <- if (length(anno) > 0) {
              DataFrame(strand, anno)
            } else {
              DataFrame(strand)
            }
            r
          })


# show-method ------------------------------------------------------------


setMethod("show", "gbRange",
          function (object) {
            lo <- length(object)
            cat(sprintf("%s of length %s\n", class(object), lo))
            if (lo == 0L){
              return(NULL)
            } 
            if (lo < 20L) {
              showme <- 
                as.data.frame(as(object, "data.frame"),
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


setAs("gbRange", "data.frame",
      function (from) {
        data.frame(as.data.frame(from),
                   as.data.frame(stringsAsFactors=FALSE,
                                 from@elementMetadata))
      })


# getters ----------------------------------------------------------------


setMethod("start", "gbRange", function (x) callNextMethod(x))


setMethod("end", "gbRange", function (x) callNextMethod(x))


setMethod("width", "gbRange", function (x) IRanges::width(x))


# shift ------------------------------------------------------------------


setMethod("shift", "gbRange",
    function (x, shift = 0L, use.names = TRUE) {
        IRanges::shift(x=x, shift=shift, use.names=use.names)
    })


# subsetting -------------------------------------------------------------


setMethod("$", "gbRange",
          function (x, name) as(x, "data.frame")[[name]])


## Ignores the 'drop' argument and behaves as if it was set to FALSE
setMethod("[", "gbRange",
          function (x, i, j, ..., drop=TRUE) {
            if (missing(j) && length(list(...)) == 0L)
              callNextMethod(x=x, i=i)
            else
              as(x, "data.frame")[i, j, ..., drop=FALSE]
          })


setMethod("[[", "gbRange",
          function (x, i, j, ...)  {
            as(x, "data.frame")[[i]]
          })
