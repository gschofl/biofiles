#' @include gbFeatureList-class.R
#' @importFrom reutils rettype retmode content
NULL

setClassUnion("gbLocationOrNull", members=c("gbLocation", "NULL"))

#' Class \code{"gbRecord"}
#'
#' \dQuote{gbRecord} is an S4 class that provides a container for data
#' parsed from a GenBank or GenPept records. For instantiation of a gbRecord
#' object use the import function \code{\link{gbRecord}}.
#' 
#' @slot seqinfo A \code{"\linkS4class{seqinfo}"} instance; This is a 
#' reference class holding the sequence as an \code{"\linkS4class{XStringSet}"}
#' instance and header of the file containing metadata as a
#' \code{"\linkS4class{gbHeader}"} object.
#' @slot features A \code{"\linkS4class{gbFeatureList}"} instance.
#' @slot contig If present, a CONTIG record.
#' 
#' @seealso
#'    The constructor, \code{\link{gbRecord}}
#'
#' @section Accessor functions:
#'  \code{\link{header}}, \code{\link{ft}}, \code{\link{ranges}}
#'
#' @name gbRecord-class
#' @rdname gbRecord-class
#' @exportClass gbRecord
new_gbRecord <- setClass(
  "gbRecord",
  slots=c(
    seqinfo  = "seqinfo",
    features = "gbFeatureList",
    contig   = "gbLocationOrNull"
  )
)


setValidity2("gbRecord", function(object) {
  # at the moment do nothing but the default checks
  TRUE
})


# constructor ------------------------------------------------------------


#' Read a GenBank/GenPept-format file.
#' 
#' Import data from GenBank/GenPept files into R, represented as an instance
#' of the \code{\linkS4class{gbRecord}} or \code{\linkS4class{gbRecordList}}
#' classes.
#' 
#' @details
#' For a sample GenBank record see
#' \url{http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html},
#' for a detailed description of the GenBank feature table format see
#' \url{http://www.ncbi.nlm.nih.gov/collab/FT/}
#'
#' @param gbk A vector of paths to GenBank format files,
#' an \code{\link[reutils]{efetch}} object containing GenBank record(s), or
#' a \code{textConnection} to a character vector that can be parsed as
#' a Genbank record.
#' @return An instance of the \code{\linkS4class{gbRecord}} or
#' \code{\linkS4class{gbRecordList}} classes.
#' @rdname gbRecord
#' @export
#' @examples
#' 
#' ### import from file
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' getHeader(x)
#' getFeatures(x)
#' 
#' ### quickly extract features as GRanges
#' ranges(x["CDS"], include=c("product", "note", "protein_id"))
#' 
#' ### import directly from NCBI
#' \dontrun{
#' require(reutils)
#' x <- gbRecord(efetch("139189709", "protein", rettype="gp", retmode="text"))
#' x
#' }
#' 
#' ## Directly subset features
#' x[[1]]
#' 
#' ## import a file containing multiple GenBank records as a
#' ## gbRecordList. With many short records it pays of to
#' ## run the parsing in paralle
#' gss_file <- system.file("extdata", "gss.gbk", package="biofiles")
#' 
#' \dontrun{
#'    require(doParallel)
#'    registerDoParallel(cores=4)
#' }
#' 
#' gss <- gbRecord(gss_file)
#' gss
#' 
gbRecord <- function(gbk) {
  if (is(gbk, "efetch")) {
    ##
    ## Parse'efetch' objects
    ##
    if (rettype(gbk) %ni% c('gb','gbwithparts','gp') || retmode(gbk) != "text") {
      stop("Use efetch with rettype='gbwithparts','gb', or'gp' and retmode='text'")
    }
    gbk_data <- readLines(textConnection(content(gbk, "text")))
    return(parse_gb_record(gbk_data))
  } else if (is(gbk, "textConnection")) {
    ##
    ## Parse textConnections
    ##
    on.exit(close(gbk))
    gbk_data <- readLines(gbk)
    return(parse_gb_record(gbk_data))
  } else if (tryCatch(all(file.exists(gbk)), error = function(e) FALSE)) {
    ##
    ## Parse GeneBank flat files  
    ##
    n <- length(gbk)
    gbk_list <- vector("list", n)
    for (i in seq_len(n)) {
      con <- file(gbk[i], open="rt")
      on.exit(close(con))
      gbk_list[[i]] <- parse_gb_record(gb_record=readLines(con))
    }
    if (length(gbk_list) == 1L) {
      return(gbk_list[[1L]])
    } else {
      return(gbRecordList(gbk_list))
    }
  } else {
    stop(paste0("'gbk' must be the path to GenBank flat file or an'efetch' ",
                "object containing GenBank records"))
  }
}


# show -------------------------------------------------------------------


#' @importFrom XVector toString subseq
.show_gbRecord <- function(x) {
  if (length(getAccession(x)) == 0L) {
    showme <- sprintf("%s instance with no features\n", sQuote(class(x)))
  } else {
    S <- .sequence(x)
    W <- getOption("width")
    showme <- paste0(
      sprintf("%s instance with %i features\n", sQuote(class(x)), length(x)),
      .header(x)$to_string(write_to_file = FALSE),
      if (length(S) != 0) {
        if (width(S) < W - 16) {
          sprintf("ORIGIN      %s\n", XVector::toString(S))
        } else {
          sprintf("ORIGIN      %s\n            ...\n            %s\n",
                  XVector::toString(
                    XVector::subseq(S, start=1, end=W-16)
                  ),
                  XVector::toString(
                    XVector::subseq(S, start=length(S[[1L]])-W+17, end=length(S[[1L]]))
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



# summary, length -----------------------------------------------------------


setMethod("summary", "gbRecord",
          function(object, n=7, ...) {
            acc  <- getAccession(object)
            len  <- getLength(object)
            type <- if (getMoltype(object) =='AA')'aa' else'bp'
            def  <- ellipsize(obj=getDefinition(object),
                              width=getOption("width") - nchar(len) - nchar(type) - 8)
            cat(sprintf("[[%s]]\n  %i %s: %s\n", acc, len, type, def), sep="")
            summary(.features(object), n=n, setoff=2)
            invisible(NULL)
          })


setMethod("length", "gbRecord", function(x) {
  length(.features(x))
})


# Internal getters ----------------------------------------------------------


setMethod('.seqinfo','gbRecord', function(x) {
  x@seqinfo
})

setMethod('.features','gbRecord', function(x) {
  x@features
})

setMethod('.contig','gbRecord', function(x) {
  x@contig
})

setMethod('.locus','gbRecord', function(x) {
  .locus(.seqinfo(x))
})

setMethod('.header','gbRecord', function(x) {
  .header(.seqinfo(x))
})

setMethod('.sequence','gbRecord', function(x) {
  .sequence(.seqinfo(x))
})

setMethod('.dbSource','gbRecord', function(x) {
  parse_dbsource(getDBSource(x))
})

setMethod(".defline", "gbRecord", function(x) {
  paste0('gi|', getGeneID(x), .dbSource(x), getAccession(x),'', getDefinition(x))
})


# getters ----------------------------------------------------------------


setMethod("getLocus", "gbRecord", function(x) getLocus(.seqinfo(x)))

setMethod("getLength", "gbRecord", function(x) getLength(.seqinfo(x)))

setMethod("getMoltype", "gbRecord", function(x) getMoltype(.seqinfo(x)))

setMethod("getTopology", "gbRecord", function(x) getTopology(.seqinfo(x)))

setMethod("getDivision", "gbRecord", function(x) getDivision(.seqinfo(x)))

setMethod("getDate", "gbRecord", function(x) getDate(.seqinfo(x)))
  
setMethod("getDefinition", "gbRecord", function(x) getDefinition(.seqinfo(x)))

setMethod("getAccession", "gbRecord", function(x) getAccession(.seqinfo(x)))

setMethod("getVersion", "gbRecord", function(x) getVersion(.seqinfo(x)))

setMethod("getGeneID", "gbRecord", function(x, db='gi') getGeneID(.seqinfo(x), db=db) )

setMethod("getDBLink", "gbRecord", function(x) getDBLink(.seqinfo(x)))

setMethod("getDBSource", "gbRecord", function(x) getDBSource(.seqinfo(x)))

setMethod("getSource", "gbRecord", function(x) getSource(.seqinfo(x)))

setMethod("getOrganism", "gbRecord", function(x) getOrganism(.seqinfo(x)))

setMethod("getTaxonomy", "gbRecord", function(x) getTaxonomy(.seqinfo(x)))

setMethod("getReference", "gbRecord", function(x) getReference(.seqinfo(x)))

setMethod("getKeywords", "gbRecord", function(x) getKeywords(.seqinfo(x)))

setMethod("getComment", "gbRecord", function(x) getComment(.seqinfo(x)))

#' @export
#' @aliases getHeader,gbRecord-method
#' @rdname getHeader-methods
setMethod("getHeader", "gbRecord", function(x) .header(x))

#' @export
#' @aliases getHeader,gbRecord-method
#' @rdname getHeader-methods
setMethod("header", "gbRecord", function(x) .header(x))

#' @export
#' @aliases getFeatures,gbRecord-method
#' @rdname getFeatures-methods
setMethod("getFeatures", "gbRecord", function(x) .features(x))

#' @export
#' @aliases getFeatures,gbRecord-method
#' @rdname getFeatures-methods
setMethod("ft", "gbRecord", function(x) .features(x))

#' @export
#' @aliases getSequence,gbRecord-method
#' @rdname getSequence-methods
setMethod("getSequence", "gbRecord",  function(x) .sequence(x))

# ' @export
# ' @aliases ranges,gbRecord-method
# ' @rdname ranges-methods
setMethod("ranges", "gbRecord",
          function(x, join = FALSE, key = TRUE, include = "none", exclude = "") {
            .make_GRanges(.features(x), join = join, include = include,
                          exclude = exclude, key = key)
          })

#' @export
#' @aliases start,gbRecord-method
#' @rdname start-methods
setMethod("start", "gbRecord", function(x, join = FALSE, drop = TRUE) {
  start(.features(x), join = join, drop = drop)
})

#' @export
#' @aliases end,gbRecord-method
#' @rdname end-methods
setMethod("end", "gbRecord", function(x, join = FALSE, drop = TRUE) {
  end(.features(x), join = join, drop = drop)
})

#' @export
#' @aliases strand,gbRecord-method
#' @rdname strand-methods
setMethod("strand", "gbRecord", function(x, join = FALSE) {
  strand(.features(x), join = join)
})

#' @export
#' @aliases width,gbRecord-method
#' @rdname width-methods
setMethod("width", "gbRecord", function(x) {
  width(.features(x))
})

#' @export
#' @aliases width,gbRecord-method
#' @rdname width-methods
setMethod("joint_width", "gbRecord", function(x) {
  joint_width(.features(x))
})

#' @export
#' @aliases dbxref,gbRecord-method
#' @rdname dbxref-methods
setMethod("dbxref", "gbRecord", function(x, db = NULL, na.rm = TRUE, ...) {
  dbxref(.features(x), db = db, na.rm = na.rm, ...)
})

#' @export
#' @aliases location,gbRecord-method
#' @rdname location-methods
setMethod("location", "gbRecord", function(x, join = FALSE) {
  lapply(.features(x), location)
})

#' @export
#' @aliases fuzzy,gbRecord-method
#' @rdname fuzzy-methods
setMethod("fuzzy", "gbRecord", function(x) {
  do.call(rbind, lapply(.features(x), fuzzy))
})

#' @export
#' @aliases index,gbRecord-method
#' @rdname index-methods
setMethod("index", "gbRecord", function(x) {
  vapply(.features(x), function(x) x@.id, 0)
})

#' @export
#' @aliases key,gbRecord-method
#' @rdname key-methods
setMethod("key", "gbRecord", function(x) {
  vapply(.features(x), function(f) f@key, "")
})

#' @export
#' @aliases qualif,gbRecord-method
#' @rdname qualif-methods
setMethod("qualif", "gbRecord", function(x, which = "", fixed = FALSE, use.names = TRUE) {
  ans <- .qual_access(.features(x), which, fixed, use.names)
  if (use.names) {
    .simplify(ans, unlist=FALSE)
  } else {
    .simplify(ans, unlist=TRUE)
  }
})


# listers ----------------------------------------------------------------


#' @export
#' @aliases listQualif,gbRecord-method
#' @rdname listQualif-methods
setMethod("listQualif", "gbRecord", function(x) {
  lapply(.features(x), listQualif)
})


# testers ----------------------------------------------------------------


#' @export
#' @aliases hasKey,gbRecord-method
#' @rdname hasKey-methods
setMethod("hasKey", "gbRecord", function(x, key) {
  vapply(.features(x), hasKey, key, FUN.VALUE=FALSE)
})

#' @export
#' @aliases hasQualif,gbRecord-method
#' @rdname hasQualif-methods
setMethod("hasQualif", "gbRecord", function(x, qualifier) {
  vapply(.features(x), hasQualif, qualifier, FUN.VALUE=FALSE)
})


# subsetting ----------------------------------------------------------------


#' Method extensions to extraction operator for gbRecord objects.
#'
#' See the documentation for the \code{\link[base]{Extract}} generic,
#' defined in the R \code{\link[base]{base-package}} for the expected behavior. 
#'
#' @seealso  \code{\link[base]{Extract}}
#' 
#' @export
#' @aliases [,gbRecord,ANY,ANY,ANY-method
#' @docType methods
#' @rdname extract-methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' 
#' ## Extract a gbFeatureList from a gbRecord:
#' x[1:4]
#' 
#' ## Extract a gbFeature
#' x[[1]]
#' 
#' ## Extract ggFeatures by Feature Key
#' x["CDS"]
#' 
setMethod("[", c("gbRecord", "ANY", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            .features(x)[i, ...]
          })


#' @export
#' @aliases [,gbRecord,missing,ANY,ANY-method
#' @rdname extract-methods
setMethod("[", c("gbRecord", "missing", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            .features(x)
          })


#' @export
#' @aliases [[,gbRecord-method
#' @rdname extract-methods
setMethod("[[", "gbRecord",
          function(x, i, j, ...) {
            .features(x)[[i]]
          })


# select, shift, revcomp ----------------------------------------------------


setMethod("select", "gbRecord", function(x, ..., keys = NULL, cols = NULL) {
  .retrieve(.select(.features(x), ..., keys = keys), cols = cols)
})


setMethod("shift", "gbRecord", function(x, shift, split=FALSE, order=FALSE) {
  .shift(x=x, shift=shift, split=split, order=order)
})
            

setMethod("revcomp", "gbRecord", function(x, order=TRUE) {
  .revcomp(x=x, order=order)
})


