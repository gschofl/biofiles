#' @include gbFeatureList-class.R
NULL

setClassUnion("XStringSetOrNull", members=c("XStringSet", "NULL"))
setClassUnion("gbLocationOrNull", members=c("gbLocation", "NULL"))

#' gbRecord
#' 
#' \dQuote{gbRecord} is an S4 class that provides a container for data
#' parsed from a GenBank record.
#'
#' @rdname gbRecord
#' @export
#' @classHierarchy
#' @classMethods
setClass("gbRecord", representation(seqinfo = "Seqinfo",
                                    locus = "character",
                                    type = "character",
                                    topology = "character",
                                    division = "character",
                                    date = "POSIXlt",
                                    version = "character",
                                    GI = "character",
                                    dblink = "character",
                                    dbsource = "character",
                                    keywords = "character",
                                    source = "character",
                                    organism = "character",
                                    lineage = "character",
                                    references = "character",
                                    comment = "character",
                                    features = "gbFeatureList",
                                    sequence = "XStringSetOrNull",
                                    contig = "gbLocationOrNull"))


setValidity("gbRecord", function (object) {
  # at the moment do nothing but the default checks
  TRUE
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
#' @autoImports
gbRecord <- function (gb, with_sequence = TRUE) {
  ##
  ## Parse 'efetch' objects
  ##
  if (is(gb, "efetch")) {
    if (gb@rettype %ni% c("gb", "gp") || gb@retmode != "text")
      stop("Must use efetch with rettype='gbwithparts','gb', or 'gp' and retmode='text'")
    
    split_gb <- unlist(strsplit(content(gb, "text"), "\n\n"))
    n <- length(split_gb)
    parsed_data <- vector("list", n)
    for (i in seq_len(n)) {
      gb_data <- unlist(strsplit(split_gb[i], "\n"))
      parsed_data[[i]] <- .parseGB(gb_data, with_sequence=with_sequence)
    }
  ##
  ## Parse random textConnections
  ##
  } else if (is(gb, "textConnection")) {
    con <- gb
    on.exit(close(con))
    parsed_data <- list(.parseGB(readLines(con), with_sequence))
  ##
  ## Parse GeneBank flat files  
  ##
  } else if (tryCatch(all(file.exists(gb)), error = function() FALSE)) {
    n <- length(gb)
    parsed_data <- vector("list", n)
    for (i in seq_len(n)) {
      con <- file(gb[i], open="rt")
      parsed_data[[i]] <- .parseGB(gb_data=readLines(con), with_sequence)
      close(con)
    }
  } else {
    stop("'gb' must be the path to a GenBank flat file or an 'efetch' object containing GenBank records")
  } 

  gbr_list <- list()
  for (gbk in parsed_data) {
    gbr <- with(gbk, 
                new("gbRecord",
                    seqinfo = header[["seqinfo"]],
                    locus = header[["locus"]],
                    type = header[["type"]],
                    topology = header[["topology"]],
                    division = header[["division"]],
                    date = header[["date"]],
                    version = header[["version"]],
                    GI = header[["GI"]],
                    dblink = header[["dblink"]],
                    dbsource = header[["dbsource"]],
                    keywords = header[["keywords"]],
                    source = header[["source"]],
                    organism = header[["organism"]],
                    lineage = header[["lineage"]],
                    references = header[["references"]],
                    comment = header[["comment"]],
                    features = features,
                    sequence = sequence,
                    contig = contig)
    )
    gbr_list <- c(gbr_list, gbr)
  }
  
  if (length(gbr_list) == 1L) 
    gbr_list[[1L]]
  else
    gbRecordList(gbr_list)
}


# show -------------------------------------------------------------------


setMethod("show", "gbRecord",
          function (object) { 
            if (is.na(accession(object))) {
              showme <- sprintf("%s instance with no features\n", sQuote(class(object)))
            } else {
              S <- sequence(object)
              W <- getOption("width")
              showme <- paste0(
                sprintf("%s instance with %i features\n", 
                        sQuote(class(object)), length(features(object))),
                sprintf("LOCUS       %s\n",
                        linebreak(paste(object@locus, seqlengths(object),
                                        if (object@type == 'AA') 'aa' else 'bp',
                                        object@type, object@topology,
                                        object@division, object@date),
                                  offset=13, FORCE=TRUE)),
                sprintf("DEFINITION  %s\n",
                        linebreak(definition(object), offset=13, FORCE=TRUE)),
                sprintf("ACCESSION   %s\n", accession(object)),
                sprintf("VERSION     %s GI:%s\n", object@version, object@GI),
                sprintf("DBLINK      Project: %s\n", object@dblink),
                if (object@type == "AA") {
                  sprintf("DBSOURCE    %s\n",
                          linebreak(object@dbsource, offset=13, FORCE=TRUE))
                },
                sprintf("KEYWORDS    %s\n",
                        linebreak(object@keywords, offset=13, FORCE=TRUE)),
                sprintf("SOURCE      %s\n",
                        linebreak(object@source, offset=13, FORCE=TRUE)),
                sprintf("  ORGANISM  %s\n",
                        linebreak(object@organism, offset=13, FORCE=TRUE)),
                sprintf("            %s\n",
                        linebreak(object@lineage, offset=13, FORCE=TRUE)),
                sprintf("REFERENCE   %s\n", object@references),
                sprintf("COMMENT     %s\n",
                        linebreak(object@comment, offset=13, FORCE=TRUE)),
                if (!is.null(S)) {
                  if (S@ranges@width[1L] < W - 14)
                    sprintf("ORIGIN      %s\n", toString(S))
                  else
                    sprintf("ORIGIN      %s\n            ...\n            %s\n",
                            toString(subseq(S, start=1, end=W - 14)),
                            toString(subseq(S, start=length(S[[1L]]) -  W + 15,
                                            end=length(S[[1L]]))))
                },
                if (!is.null(object@contig)) {
                  sprintf("CONTIG      %s\n", linebreak(as(object@contig, "character"),
                                                        offset=13, split=",", FORCE=TRUE))
                })
            }
            
            cat(showme)
          })


# summary ----------------------------------------------------------------


setMethod("summary", "gbRecord",
          function (object, n=7, ...) {
            si <- seqinfo(object)
            acc <- seqnames(si)
            len <- unname(seqlengths(si))
            type <- if (object@type == 'AA') 'aa' else 'bp'
            def <- 
              ellipsize(obj=unname(genome(si)),
                        width=getOption("width") - nchar(len) - nchar(type) - 8)
            cat(sprintf("[[%s]]\n  %i %s: %s\n", acc, len, type, def))
            obj.features <- features(object)
            if (length(obj.features) > 2*n) {
              head <- head(obj.features, n=n)
              tail <- tail(obj.features, n=n)
              x <- lapply(head, summary)
              cat("...\n")
              x <- lapply(tail, summary)  
            } else  {
              x <- lapply(obj.features, summary)
            }
            
            return(invisible(NULL))
          })


# getters ----------------------------------------------------------------


setMethod("seqinfo", "gbRecord",
          function (x) x@seqinfo)


setMethod("seqlengths", "gbRecord",
          function (x) seqlengths(seqinfo(x)))


setMethod("accession", "gbRecord", 
          function (x) seqnames(x@seqinfo))


setMethod("definition", "gbRecord", 
          function (x) genome(x@seqinfo))


setMethod("features", "gbRecord", 
          function (x) x@features)


setMethod("sequence", "gbRecord", 
          function (x) x@sequence)


setMethod("ranges", "gbRecord",
          function (x, join = FALSE, key = TRUE, include = "none", exclude = "") {
            .make_GRanges(x@features, join = join, include = include,
                          exclude = exclude, key = key)
          })


setMethod("start", "gbRecord",
          function (x, join = FALSE, drop = TRUE) {
            start(features(x), join = join, drop = drop)
          })


setMethod("end", "gbRecord",
          function (x, join = FALSE, drop = TRUE) {
            end(features(x), join = join, drop = drop)
          })


setMethod("strand", "gbRecord",
          function (x, join = FALSE, drop = TRUE) {
            strand(features(x), join = join, drop = drop)
          })


# listers ----------------------------------------------------------------


setMethod("listQualif", "gbRecord", 
          function (x) {
            lapply(x@features, listQualif)
          })


# subsetting ----------------------------------------------------------


setMethod("[[", c("gbRecord", "character", "missing"),
          function(x, i) slot(object=x, name=i))


setMethod("$", "gbRecord",
          function(x, name) slot(object=x, name=name))


# select -----------------------------------------------------------------


setMethod("select", "gbRecord",
          function (x, ..., keys = NULL, cols = NULL) {
            ans <- x@features
            ans <- .select(ans, ..., keys = keys)
            ans <- .retrieve(ans, cols = cols)
            ans
          })


# shift ------------------------------------------------------------------


setMethod("shift", "gbRecord",
          function(x, shift, split=FALSE, order=FALSE, updateDb=FALSE)
            .shift_features(x=x, shift=shift, split=split, order=order))

