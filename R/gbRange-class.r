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
              return( .Object )
            }
            if (any(strand %ni% c(1,-1,NA))) {
              stop("Strand must be encoded as 1 (plus strand), -1 (minus strand), or NA")
            }
            
            anno <- as.list(unlist(list(...))) 
            z <- vapply(c(list(start, width, strand), anno), length, numeric(1))
            if (length(unique(z)) != 1L) {
              stop("Arguments have unequal length")
            }
            
            r <- callNextMethod(.Object, start = start, width = width,
                                NAMES = anno[["names"]])
            
            anno[["names"]] <- NULL
            r@elementMetadata <- if (length(anno) > 0) {
              DataFrame(strand, anno)
            } else {
              DataFrame(strand)
            }
            r
          })


# show-method ------------------------------------------------------------


#' @autoImports
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


#' @export
setAs("gbRange", "data.frame",
      function (from) {
        # there is no proper setAs() method defined for coercion from
        # IRanges to data.frames in IRanges, only as.data.frame.
        # So we can't use callNextMethod() here.
        op <- options(stringsAsFactors = FALSE)
        df <- data.frame(IRanges::as.data.frame(from),
                         as(from@elementMetadata, "data.frame"))
        options(op)
        df
      })


# getters ----------------------------------------------------------------


setMethod("start", "gbRange", function (x) callNextMethod(x))


setMethod("end", "gbRange", function (x) callNextMethod(x))


setMethod("width", "gbRange", function (x) IRanges::width(x))


setMethod("strand", "gbRange", 
          function (x) elementMetadata(x)[["strand"]])


setMethod("range", "gbRange", 
          function (x) {
            elementMetadata(x) <- NULL
            x
          })


setMethod("annotation", "gbRange", 
          function (x)  elementMetadata(x))


setMethod("sequence", "gbRange",
          function (x, seq, ...) {
            
            if (is(seq, "gbRecord")) {
              seq <- sequence(seq)
            }
            
            if (is(seq, "DNAStringSet")) {
              if (length(seq) > 1) {
                warning("'seq' contains multiple sequences. Only the first will be used")
              }
              seq <- seq[[1]]
            }
            
            if (!is(seq, "DNAString")) {
              stop("'seq' must be a 'DNAString' object")
            }
            
            x <- sort(x)
            start <- biofiles::start(x)
            end <- biofiles::end(x)
            strand <- biofiles::strand(x)
            seqs <- as(Views(seq, start, end), "DNAStringSet")
            o <- order(c(start[strand == -1], start[strand == 1]))
            
            append(reverseComplement(seqs[strand == -1]),
                   seqs[strand == 1])[o]
          })


# shift ------------------------------------------------------------------


setMethod("shift", "gbRange",
    function (x, shift = 0L, use.names = TRUE) {
        IRanges::shift(x=x, shift=shift, use.names=use.names)
    })


# subsetting -------------------------------------------------------------


#' @export
setMethod("$", "gbRange",
          function (x, name) as(x, "data.frame")[[name]])


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

