#' @include gbFeatureTable-class.R
#' @importFrom reutils rettype retmode content
NULL

setClassUnion("gbLocationOrNull", members = c("gbLocation", "NULL"))

#' Class \code{"gbRecord"}
#'
#' \dQuote{gbRecord} is an S4 class that provides a container for data
#' parsed from GenBank, GenPept, Embl or IMGT records. For instantiation of a
#' gbRecord object use the import function \code{\link{gbRecord}}.
#' 
#' @slot seqinfo A \code{"\linkS4class{seqinfo}"} instance; This is a 
#' reference class holding the sequence as an \code{"\linkS4class{XStringSet}"}
#' instance and header of the file containing metadata as a
#' \code{"\linkS4class{gbHeader}"} object.
#' @slot features A \code{"\linkS4class{gbFeatureTable}"} instance.
#' @slot contig If present, a CONTIG record.
#' 
#' @seealso
#'    The constructor, \code{\link{gbRecord}}
#'
#' @section Accessor functions:
#'  \code{\link{header}}, \code{\link{ft}}, \code{\link{ranges}}
#'
#' @export
new_gbRecord <- setClass(
  "gbRecord",
  slots = c(
    seqinfo  = "seqinfo",
    features = "gbFeatureTable",
    contig   = "gbLocationOrNull"
  )
)


setValidity2("gbRecord", function(object) {
  # at the moment do nothing but the default checks
  TRUE
})


# constructor ------------------------------------------------------------


#' Read a GenBank/GenPept or Embl format file.
#' 
#' Import data from GenBank/GenPept, Embl, or IMGT/HLA flat files into R,
#' represented as an instance of the \code{\linkS4class{gbRecord}} or
#' \code{\linkS4class{gbRecordList}} classes.
#' 
#' @details
#' For a sample GenBank record see
#' \url{http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html},
#' for a detailed description of the GenBank feature table format see
#' \url{http://www.ncbi.nlm.nih.gov/collab/FT/}.
#' 
#' For a description of the EMBL flat file format see
#' \url{ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt}.
#' 
#' For a description of the format and conventions of IMGT/HLA flat files
#' see \url{http://www.ebi.ac.uk/ipd/imgt/hla/docs/manual.html}.
#'
#' @note
#' The \code{\linkS4class{gbRecord}} class is modelled after the Genbank flat
#' file format. Both Embl and IMGT/HLA files do not fit this model perfectly,
#' so some pretty arbitrary choices were made to make the data from these files
#' fitr the model.
#'
#' @param rcd A vector of paths to GenBank/Embl format records,
#' an \code{\link[reutils]{efetch}} object containing GenBank record(s), or
#' a \code{textConnection} to a character vector that can be parsed as
#' a Genbank or Embl record.
#' @param progress Print a nice progress bar if parsing multiple Genbank records.
#' (This will not work if you process the records in parallel.)
#' @return An instance of the \code{\linkS4class{gbRecord}} or
#' \code{\linkS4class{gbRecordList}} classes.
#' @seealso
#'  \code{\link{genomeRecordFromNCBI}}
#' @export
#' @examples
#' \dontrun{
#' ### import from file
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package = "biofiles")
#' x <- gbRecord(gbk_file)
#' }
#' 
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' getHeader(x)
#' getFeatures(x)
#' 
#' ### quickly extract features as GRanges
#' ranges(x["CDS"], include = c("product", "note", "protein_id"))
#' 
#' ## Directly subset features
#' x[[1]]
#' 
#' ### import directly from NCBI
#' \dontrun{
#' require(reutils)
#' x <- gbRecord(efetch("139189709", "protein", rettype = "gp", retmode = "text"))
#' x
#' }
#' 
#' ## import a file containing multiple GenBank records as a
#' ## gbRecordList. With many short records it pays of to
#' ## run the parsing in parallel
#' \dontrun{
#' gss_file <- system.file("extdata", "gss.gbk", package = "biofiles")
#' library(doParallel)
#' registerDoParallel(cores = 4)
#' gss <- gbRecord(gss_file)
#' gss
#' }
gbRecord <- function(rcd, progress = FALSE) {
  if (missing(rcd)) {
    ## instantiate an empty gbRecord
    return(new_gbRecord())
  }
  if (is(rcd, "efetch")) {
    assert_that(retmode(rcd) == "text", rettype(rcd) %is_in% c('gb','gbwithparts','gp'))
    rcd <- content(rcd, "textConnection")
  }
  if (is(rcd, "textConnection")) {
    on.exit(close(rcd))
    return(parse_record(rcd = readLines(rcd), progress = progress))
  } else if (tryCatch(all(file.exists(rcd)), error = function(e) FALSE)) {
    n <- length(rcd)
    rcd_list <- vector("list", n)
    for (i in seq_len(n)) {
      con <- file(rcd[i], open = "rt")
      on.exit(close(con))
      rcd_list[[i]] <- parse_record(rcd = readLines(con), progress = progress)
    }
    if (length(rcd_list) == 1L) {
      return(rcd_list[[1L]])
    } else {
      return(gbRecordList(rcd_list))
    }
  } else {
    stop(paste0("'rcd' must be the path to GenBank or Embl flat file or an 'efetch' ",
                "object containing GenBank records"))
  }
}


# show -------------------------------------------------------------------


#' @importFrom XVector toString subseq
show_gbRecord <- function(x) {
  if (.header(x)$is_empty()) {
    showme <- sprintf("An empty object of class %s.\n", sQuote(class(x)))
  } else {
    S <- .sequence(x)
    W <- getOption("width")
    showme <- paste0(
      sprintf("An object of class %s, with %i features\n", sQuote(class(x)), length(x)),
      .header(x)$to_string(write_to_file = FALSE),
      if (length(S) != 0) {
        if (width(S) < W - 16) {
          sprintf("ORIGIN      %s\n", XVector::toString(S))
        } else {
          sprintf("ORIGIN      %s\n            ...\n            %s\n",
                  toString(subseq(S, start = 1, end = W-16)),
                  toString(subseq(S, start = length(S[[1L]])-W+17, end = length(S[[1L]]))))
        }
      },
      sprintf("CONTIG      %s\n", linebreak(as(.contig(x), "character"), 
                                            offset = 13, split = ",", FORCE = TRUE))
    )
  }
  cat(showme, sep = "\n")
}

setMethod("show", "gbRecord", function(object) { 
  show_gbRecord(object)
})


# summary, length -----------------------------------------------------------


#' @param n How many elements should be summarized in head and tail.
#' @rdname summary-methods
setMethod("summary", "gbRecord",
          function(object, n = 7, ...) {
            acc  <- getAccession(object)
            len  <- getLength(object)
            type <- if (getMoltype(object) == 'AA')'aa' else'bp'
            def  <- ellipsize(obj = getDefinition(object),
                              width = getOption("width") - nchar(len) - nchar(type) - 8)
            cat(sprintf("[[%s]]\n  %i %s: %s\n", acc, len, type, def), sep = "")
            summary(.features(object), n = n, setoff = 2)
            invisible(NULL)
          })


#' Get the number of gbFeatures.
#'
#' @param x A \code{gbRecord}
#' @return An integer
#' @export
#' @docType methods
#' @rdname length-methods
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
  paste0('gi|', getGeneID(x), .dbSource(x), getAccession(x), ' ', getDefinition(x))
})


# getters ----------------------------------------------------------------

#' @rdname accessors
setMethod("getLocus", "gbRecord", function(x) getLocus(.seqinfo(x)))
#' @rdname accessors
setMethod("getLength", "gbRecord", function(x) getLength(.seqinfo(x)))
#' @rdname accessors
setMethod("getMoltype", "gbRecord", function(x) getMoltype(.seqinfo(x)))
#' @rdname accessors
setMethod("getTopology", "gbRecord", function(x) getTopology(.seqinfo(x)))
#' @rdname accessors
setMethod("getDivision", "gbRecord", function(x) getDivision(.seqinfo(x)))
#' @rdname accessors
setMethod("getDate", "gbRecord", function(x) getDate(.seqinfo(x)))
#' @rdname accessors
setMethod("getDefinition", "gbRecord", function(x) getDefinition(.seqinfo(x)))
#' @rdname accessors
setMethod("getAccession", "gbRecord", function(x) getAccession(.seqinfo(x)))
#' @rdname accessors
setMethod("getVersion", "gbRecord", function(x) getVersion(.seqinfo(x)))
#' @rdname accessors
setMethod("getGeneID", "gbRecord", function(x, db = 'gi') getGeneID(.seqinfo(x), db = db) )
#' @rdname accessors
setMethod("getDBLink", "gbRecord", function(x) getDBLink(.seqinfo(x)))
#' @rdname accessors
setMethod("getDBSource", "gbRecord", function(x) getDBSource(.seqinfo(x)))
#' @rdname accessors
setMethod("getSource", "gbRecord", function(x) getSource(.seqinfo(x)))
#' @rdname accessors
setMethod("getOrganism", "gbRecord", function(x) getOrganism(.seqinfo(x)))
#' @rdname accessors
setMethod("getTaxonomy", "gbRecord", function(x) getTaxonomy(.seqinfo(x)))
#' @rdname accessors
setMethod("getReference", "gbRecord", function(x) getReference(.seqinfo(x)))
#' @rdname accessors
setMethod("getKeywords", "gbRecord", function(x) getKeywords(.seqinfo(x)))
#' @rdname accessors
setMethod("getComment", "gbRecord", function(x) getComment(.seqinfo(x)))

#' @rdname getHeader-methods
setMethod("getHeader", "gbRecord", function(x) .header(x))

#' @rdname getHeader-methods
setMethod("header", "gbRecord", function(x) .header(x))

#' @rdname getFeatures-methods
setMethod("getFeatures", "gbRecord", function(x) .features(x))

#' @rdname getFeatures-methods
setMethod("ft", "gbRecord", function(x) .features(x))

#' @rdname getSequence-methods
setMethod("getSequence", "gbRecord",  function(x) .sequence(x))

#' @describeIn ranges
setMethod("ranges", "gbRecord",
          function(x, join = FALSE, key = TRUE, include = "none", exclude = "") {
            .GRanges(x = .features(x), join = join, include = include,
                     exclude = exclude, key = key)
          })

#' @describeIn start
setMethod("start", "gbRecord", function(x, join = FALSE) {
  start(.features(x), join = join)
})

#' @describeIn end
setMethod("end", "gbRecord", function(x, join = FALSE) {
  end(.features(x), join = join)
})

#' @describeIn strand
setMethod("strand", "gbRecord", function(x, join = FALSE) {
  strand(.features(x), join = join)
})

#' @describeIn span
setMethod("span", "gbRecord", function(x, join = FALSE) {
  span(.features(x), join = join)
})

#' @rdname dbxref-methods
setMethod("dbxref", "gbRecord", function(x, db = NULL, ...) {
  dbxref(.features(x), db = db, ...)
})

#' @rdname location-methods
setMethod("location", "gbRecord", function(x, join = FALSE) {
  lapply(.features(x), location)
})

#' @describeIn fuzzy
setMethod("fuzzy", "gbRecord", function(x) {
  do.call(rbind, lapply(.features(x), fuzzy))
})

#' @rdname index-methods
setMethod("index", "gbRecord", function(x) {
  index(.features(x))
})

#' @rdname key-methods
setMethod("key", "gbRecord", function(x) {
  vapply(.features(x), slot, name = 'key', FUN.VALUE = "")
})

#' @rdname qualif-methods
setMethod("qualif", "gbRecord", function(x, which = "", fixed = FALSE, use.names = TRUE) {
  ans <- .qual_access(.features(x), which, fixed, use.names)
  if (use.names) {
    .simplify(ans, unlist = FALSE)
  } else {
    .simplify(ans, unlist = TRUE)
  }
})


# listers ----------------------------------------------------------------


#' @rdname qualifList-methods
setMethod("qualifList", "gbRecord", function(x) {
  lapply(.features(x), qualifList)
})

#' @rdname qualifTable-methods
setMethod("qualifTable", "gbRecord", function(x) {
  tbls <- tbl_qual(.features(x))
  Reduce(tbl_merge, tbls)
})

#' @rdname featureTable-methods
setMethod("featureTable", "gbRecord", function(x) {
  table(key(.features(x)))
})

# testers ----------------------------------------------------------------


#' @rdname hasKey-methods
setMethod("hasKey", "gbRecord", function(x, key) {
  vapply(.features(x), hasKey, key, FUN.VALUE = FALSE)
})

#' @rdname hasQualif-methods
setMethod("hasQualif", "gbRecord", function(x, qualifier) {
  vapply(.features(x), hasQualif, qualifier, FUN.VALUE = FALSE)
})


# subsetting ----------------------------------------------------------------


#' Method extensions to extraction operator for gbRecord objects.
#'
#' See the documentation for the \code{\link[base]{Extract}} generic,
#' defined in the R \code{\link[base]{base-package}} for the expected behavior.
#'  
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' or \code{\linkS4class{gbRecord}} object.
#' @param i indices specifying elements to extract. With \code{gbFeatureTable}s
#' and \code{gbRecord}s, a \code{character} index is matched against feature keys;
#' \code{gbFeature}s a \code{character} index is matched against qualifiers.
#' @param j Not used.
#' @param ... Not used.
#' @param drop Not used.
#' @return A \code{\linkS4class{gbFeatureTable}} object or elements of a
#' \code{\linkS4class{gbFeature}} object.
#' @seealso  \code{\link[base]{Extract}}
#' @export
#' @docType methods
#' @rdname extract-methods
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' 
#' ## Extract a gbFeatureTable from a gbRecord:
#' x[1:4]
#' 
#' ## Extract a gbFeature
#' x[[1]]
#' 
#' ## Extract ggFeatures by Feature Key
#' x["CDS"]
#' 
setMethod("[", c("gbRecord", "ANY", "missing", "ANY"), function(x, i, j, ..., drop = TRUE) {
  .features(x)[i, ...]
})

#' @rdname extract-methods
setMethod("[", c("gbRecord", "missing", "missing", "ANY"), function(x, i, j, ..., drop = TRUE) {
  .features(x)
})

#' @rdname extract-methods
setMethod("[[", "gbRecord", function(x, i, j, ...) {
  .features(x)[[i]]
})


# select, shift, revcomp ----------------------------------------------------


#' @rdname manip-methods
setMethod("filter", "gbRecord", function(x, ..., .cols = NULL) {
  ans <- .filter(.features(x), ..., .cols = .cols)
  if (is.null(.cols)) {
    ans <- as(ans, 'gbRecord')
  }
  ans
})

#' @rdname manip-methods
setMethod("select", "gbRecord", function(x, ..., .cols = NULL) {
  .select(.features(x), ..., .cols = .cols)
})

#' @describeIn shift
setMethod("shift", "gbRecord", function(x, shift, split = FALSE, order = TRUE) {
  .shift(x = x, shift = shift, split = split, order = order)
})
            
#' @rdname revcomp-methods
setMethod("revcomp", "gbRecord", function(x, order = TRUE) {
  .revcomp(x = x, order = order)
})


