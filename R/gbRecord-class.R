#' @include gbFeatureList-class.R
NULL

setClassUnion("gbLocationOrNull", members=c("gbLocation", "NULL"))

#' gbRecord-class
#' 
#' \dQuote{gbRecord} is an S4 class that provides a container for data
#' parsed from a GenBank record.
#'
#' @name gbRecord-class
#' @rdname gbRecord-class
#' @exportClass gbRecord
new_gbRecord <- setClass(
  "gbRecord",
  slots=c(
    seqinfo = "seqinfo",
    features = "gbFeatureList",
    contig = "gbLocationOrNull"
  )
)


setValidity2("gbRecord", function(object) {
  # at the moment do nothing but the default checks
  TRUE
})


# Internal getters ----------------------------------------------------------


setMethod('.seqinfo', 'gbRecord', function(x) {
  x@seqinfo
})

setMethod('.features', 'gbRecord', function(x) {
  x@features
})

setMethod('.contig', 'gbRecord', function(x) {
  x@contig
})

setMethod('.locus', 'gbRecord', function(x) {
  .locus(.seqinfo(x))
})

setMethod('.header', 'gbRecord', function(x) {
  .header(.seqinfo(x))
})

setMethod('.sequence', 'gbRecord', function(x) {
  .sequence(.seqinfo(x))
})

setMethod('.dbSource', 'gbRecord', function(x) {
  parse_dbsource(getDBSource(x))
})

setMethod(".defline", "gbRecord", function(x) {
  paste0('gi|', getGeneID(x), .dbSource(x), getAccession(x), ' ', getDefinition(x))
})


# constructor ------------------------------------------------------------


#' \code{gbRecord} or \code{\linkS4class{gbRecordList}} objects can be
#' construced by parsing GenBank flat files or an \code{\linkS4class{efetch}}
#' object containing one or more GenBank records.
#'
#' @details
#' For a description of the GenBank format see
#' \url{http://www.ncbi.nlm.nih.gov/collab/FT/}
#'
#' @param gb A vector of paths to GenBank record files,
#' an \code{\linkS4class{efetch}} object containing GenBank record(s), or
#' a \code{textConnection} to a character vector that can be parsed as
#' a Genbank record.
#' @param with_sequence Include sequence information if avaliable.
#' @return A \code{\linkS4class{gbRecord}} or
#' \code{\linkS4class{gbRecordList}} 
#' @export
gbRecord <- function(gb, with_sequence = TRUE) {
  ##
  ## Parse 'efetch' objects
  ##
  if (is(gb, "efetch")) {
    stopifnot(require(reutils))
    if (rettype(gb) %ni% c('gb', 'gbwithparts', 'gp') || retmode(gb) != "text") {
      stop("Must use efetch with rettype='gbwithparts','gb', or 'gp' and retmode='text'")
    }
    gb <- usplit(content(gb, "text"), "\n\n")
    n <- length(gb)
    gb_list <- vector("list", n)
    for (i in seq_len(n)) {
      gb_data <- usplit(gb[i], "\n")
      gb_list[[i]] <- .parseGbRecord(gb_data, with_sequence)
    }
  ##
  ## Parse random textConnections
  ##
  } else if (is(gb, "textConnection")) {
    con <- gb
    on.exit(close(con))
    gb_list <- list(.parseGbRecord(gb_data=readLines(con), with_sequence))
  ##
  ## Parse GeneBank flat files  
  ##
  } else if (tryCatch(all(file.exists(gb)), error = function() FALSE)) {
    n <- length(gb)
    gb_list <- vector("list", n)
    for (i in seq_len(n)) {
      con <- file(gb[i], open="rt")
      gb_list[[i]] <- .parseGbRecord(gb_data=readLines(con), with_sequence)
      close(con)
    }
  } else {
    stop(paste0("'gb' must be the path to a GenBank flat file or an 'efetch' ",
                "object containing GenBank records"))
  }
  if (length(gb_list) == 1L) {
    gb_list[[1L]]
  } else {
    gbRecordList(gb_list)
  }
}


# show -------------------------------------------------------------------


#' @importFrom Biostrings width
#' @importFrom XVector toString
#' @importFrom XVector subseq
.show_gbRecord <- function(x) {
  if (length(getAccession(x)) == 0L) {
    showme <- sprintf("%s instance with no features\n", sQuote(class(x)))
  } else {
    S <- .sequence(x)
    W <- getOption("width")
    type <- getMoltype(x)
    showme <- paste0(
      sprintf("%s instance with %i features\n", sQuote(class(x)), length(.features(x))),
      .header(x)$to_string(write_to_file = FALSE),
      if (length(S) != 0) {
        if (Biostrings::width(S) < W - 14) {
          sprintf("ORIGIN      %s\n", XVector::toString(S))
        } else {
          sprintf("ORIGIN      %s\n            ...\n            %s\n",
                  XVector::toString(
                    XVector::subseq(S, start=1, end=W - 14)
                  ),
                  XVector::toString(
                    XVector::subseq(S, start=length(S[[1L]]) -  W + 15,
                                    end=length(S[[1L]]))
                  )
          )
        }
      },
      if (!is.null(x)) {
        sprintf("CONTIG      %s\n", linebreak(as(.contig(x), "character"),
                                              offset=13, split=",", FORCE=TRUE))
      })
  }
  
  cat(showme)
}


setMethod("show", "gbRecord",
          function(object) { 
            .show_gbRecord(object)
          })


# summary ----------------------------------------------------------------


setMethod("summary", "gbRecord",
          function(object, n=7, ...) {
            acc  <- getAccession(object)
            len  <- getLength(object)
            type <- if (getMoltype(object) == 'AA') 'aa' else 'bp'
            def  <- ellipsize(obj=getDefinition(object),
                              width=getOption("width") - nchar(len) - nchar(type) - 8)
            cat(sprintf("[[%s]]\n  %i %s: %s\n", acc, len, type, def), sep="")
            summary(.features(object), n=n, setoff=2)
            
            return(invisible(NULL))
          })


# getters ----------------------------------------------------------------

setMethod("getLocus", "gbRecord", function(x) getLocus(.seqinfo(x)) )

setMethod("getLength", "gbRecord", function(x) getLength(.seqinfo(x)) )

setMethod("getMoltype", "gbRecord", function(x) getMoltype(.seqinfo(x)) )

setMethod("getTopology", "gbRecord", function(x) getTopology(.seqinfo(x)) )

setMethod("getDivision", "gbRecord", function(x) getDivision(.seqinfo(x)) )

setMethod("getDate", "gbRecord", function(x) getDate(.seqinfo(x)) )
  
setMethod("getDefinition", "gbRecord", function(x) getDefinition(.seqinfo(x)) )

setMethod("getAccession", "gbRecord", function(x) getAccession(.seqinfo(x)) )

setMethod("getVersion", "gbRecord", function(x) getVersion(.seqinfo(x)) )

setMethod("getGeneID", "gbRecord", function(x, db='gi') getGeneID(.seqinfo(x), db=db) )

setMethod("getDBLink", "gbRecord", function(x) getDBLink(.seqinfo(x)) )

setMethod("getDBSource", "gbRecord", function(x) getDBSource(.seqinfo(x)) )

setMethod("getSource", "gbRecord", function(x) getSource(.seqinfo(x)) )

setMethod("getOrganism", "gbRecord", function(x) getOrganism(.seqinfo(x)) )

setMethod("getTaxonomy", "gbRecord", function(x) getTaxonomy(.seqinfo(x)) )

setMethod("getReference", "gbRecord", function(x) getReference(.seqinfo(x)) )

setMethod("getKeywords", "gbRecord", function(x) getKeywords(.seqinfo(x)) )

setMethod("getComment", "gbRecord", function(x) getComment(.seqinfo(x)) )

setMethod("getFeatures", "gbRecord", function(x) .features(x))

setMethod("getSequence", "gbRecord",  function(x) .sequence(x))


setMethod("ranges", "gbRecord",
          function(x, join = FALSE, key = TRUE, include = "none", exclude = "") {
            .make_GRanges(.features(x), join = join, include = include,
                          exclude = exclude, key = key)
          })


setMethod("start", "gbRecord",
          function(x, join = FALSE, drop = TRUE) {
            start(.features(x), join = join, drop = drop)
          })


setMethod("end", "gbRecord",
          function(x, join = FALSE, drop = TRUE) {
            end(.features(x), join = join, drop = drop)
          })


setMethod("strand", "gbRecord",
          function(x, join = FALSE) {
            strand(.features(x), join = join)
          })


setMethod("width", "gbRecord",
          function(x, join = FALSE) {
            width(.features(x), join = join)
          })


# listers ----------------------------------------------------------------


setMethod("listQualif", "gbRecord", 
          function(x) {
            lapply(.features(x), listQualif)
          })


# select -----------------------------------------------------------------


setMethod("select", "gbRecord",
          function(x, ..., keys = NULL, cols = NULL) {
            ans <- .features(x)
            ans <- .select(ans, ..., keys = keys)
            ans <- .retrieve(ans, cols = cols)
            ans
          })


# shift ------------------------------------------------------------------


setMethod("shift", "gbRecord", function(x, shift, split=FALSE, order=FALSE) {
  .shift(x=x, shift=shift, split=split, order=order)
})
            

setMethod("revcomp", "gbRecord", function(x, order=TRUE) {
  .revcomp(x=x, order=order)
})


