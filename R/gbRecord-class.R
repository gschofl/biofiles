#' @include gbFeatureList-class.R
NULL

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
new_gbRecord <- setClass(
  "gbRecord",
  slots=c(
    seqinfo = "seqinfo",
    features = "gbFeatureList",
    contig = "gbLocationOrNull"
  )
)



setValidity2("gbRecord", function (object) {
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
gbRecord <- function (gb, with_sequence = TRUE)
{
  ##
  ## Parse 'efetch' objects
  ##
  if (is(gb, "efetch"))
  {
    if (gb@rettype %ni% c('gb','gbwithparts','gp') || gb@retmode != "text")
      stop("Must use efetch with rettype='gbwithparts','gb', or 'gp' and retmode='text'")
    
    split_gb <- unlist(strsplit(Rentrez::content(gb, "text"), "\n\n"))
    n <- length(split_gb)
    gbr_list <- vector("list", n)
    for (i in seq_len(n))
    {
      gb_data <- unlist(strsplit(split_gb[i], "\n"))
      gbr_list[[i]] <- .parseGbRecord( gb_data, with_sequence )
    }
  ##
  ## Parse random textConnections
  ##
  }
  else if (is(gb, "textConnection"))
  {
    con <- gb
    on.exit(close(con))
    gbr_list <- list( .parseGbRecord(readLines(con), with_sequence) )
  ##
  ## Parse GeneBank flat files  
  ##
  }
  else if (tryCatch(all(file.exists(gb)), error = function() FALSE))
  {
    n <- length(gb)
    gbr_list <- vector("list", n)
    for (i in seq_len(n))
    {
      con <- file(gb[i], open="rt")
      gbr_list[[i]] <- .parseGbRecord( gb_data=readLines(con), with_sequence )
      close(con)
    }
  }
  else
    stop(paste0("'gb' must be the path to a GenBank flat file or an 'efetch' ",
                "object containing GenBank records"))
  
  if (length(gbr_list) == 1L) 
    gbr_list[[1L]]
  else
    gbRecordList( gbr_list )
}


# show -------------------------------------------------------------------

#' @importFrom XVector toString  subseq
.show_gbRecord <- function (x) {
  if (is.na(getAccession(x))) {
    showme <- sprintf("%s instance with no features\n", sQuote(class(x)))
  } else {
    S <- getSequence(x)
    W <- getOption("width")
    type <- getMoltype(x)
    showme <- paste0(
      sprintf("%s instance with %i features\n", sQuote(class(x)), length(getFeatures(x))),
      sprintf("LOCUS       %s\n",
              linebreak(paste(getLocus(x), getLength(x),
                              if (type == 'AA') 'aa' else 'bp',
                              type, getTopology(x), getDivision(x),
                              unname(getDate(x)[2])),
                        offset=13, FORCE=TRUE)),
      sprintf("DEFINITION  %s\n",
              linebreak(getDefinition(x), offset=13, FORCE=TRUE)),
      sprintf("ACCESSION   %s\n", getAccession(x)),
      sprintf("VERSION     %s GI:%s\n", getVersion(x), getGeneID(x)),
      sprintf("DBLINK      Project: %s\n", getDBLink(x)),
      if (type == "AA") {
        sprintf("DBSOURCE    %s\n",
                linebreak(getDBSource(x), offset=13, FORCE=TRUE))
      },
      sprintf("KEYWORDS    %s\n",
              linebreak(getKeywords(x), offset=13, FORCE=TRUE)),
      sprintf("SOURCE      %s\n",
              linebreak(getSource(x), offset=13, FORCE=TRUE)),
      sprintf("  ORGANISM  %s\n",
              linebreak(getOrganism(x), offset=13, FORCE=TRUE)),
      sprintf("            %s\n",
              linebreak(getTaxonomy(x), offset=13, FORCE=TRUE)),
      sprintf("REFERENCE   %s\n", getReference(x)),
      sprintf("COMMENT     %s\n",
              linebreak(getComment(x), offset=13, FORCE=TRUE)),
      if (length(S) != 0) {
        if (S@ranges@width[1L] < W - 14)
          sprintf("ORIGIN      %s\n", XVector::toString(S))
        else
          sprintf("ORIGIN      %s\n            ...\n            %s\n",
                  XVector::toString(
                    XVector::subseq(S, start=1, end=W - 14)
                  ),
                  XVector::toString(
                    subseq(S, start=length(S[[1L]]) -  W + 15,
                           end=length(S[[1L]]))
                  )
          )
      },
      if (!is.null(x)) {
        sprintf("CONTIG      %s\n", linebreak(as(x@contig, "character"),
                                              offset=13, split=",", FORCE=TRUE))
      })
  }
  
  cat(showme)
}


setMethod("show", "gbRecord",
          function (object) { 
            .show_gbRecord(object)
          })


# summary ----------------------------------------------------------------


#' @autoImports
setMethod("summary", "gbRecord",
          function (object, n=7, ...) {
            acc <- getAccession(object)
            len <- getLength(object)
            type <- if (getMoltype(object) == 'AA') 'aa' else 'bp'
            def <- 
              ellipsize(obj=getDefinition(object),
                        width=getOption("width") - base::nchar(len) - base::nchar(type) - 8)
            cat(sprintf("[[%s]]\n  %i %s: %s\n", acc, len, type, def), sep="")
            summary(getFeatures(object), n=n, setoff=2)
            
            return(invisible(NULL))
          })


# getters ----------------------------------------------------------------

setMethod("getLocus", "gbRecord", function (x) getLocus(x@seqinfo) )

setMethod("getLength", "gbRecord", function (x) getLength(x@seqinfo) )

setMethod("getMoltype", "gbRecord", function (x) getMoltype(x@seqinfo) )

setMethod("getTopology", "gbRecord", function (x) getTopology(x@seqinfo) )

setMethod("getDivision", "gbRecord", function (x) getDivision(x@seqinfo) )

setMethod("getDate", "gbRecord", function (x) getDate(x@seqinfo) )
  
setMethod("getDefinition", "gbRecord", function (x) getDefinition(x@seqinfo) )

setMethod("getAccession", "gbRecord", function (x) getAccession(x@seqinfo) )

setMethod("getVersion", "gbRecord", function (x) getVersion(x@seqinfo) )

setMethod("getGeneID", "gbRecord", function (x, db='gi') getGeneID(x@seqinfo, db=db) )

setMethod("getDBLink", "gbRecord", function (x) getDBLink(x@seqinfo) )

setMethod("getDBSource", "gbRecord", function (x) getDBSource(x@seqinfo) )

setMethod("getSource", "gbRecord", function (x) getSource(x@seqinfo) )

setMethod("getOrganism", "gbRecord", function (x) getOrganism(x@seqinfo) )

setMethod("getTaxonomy", "gbRecord", function (x) getTaxonomy(x@seqinfo) )

setMethod("getReference", "gbRecord", function (x) getReference(x@seqinfo) )

setMethod("getKeywords", "gbRecord", function (x) getKeywords(x@seqinfo) )

setMethod("getComment", "gbRecord", function (x) getComment(x@seqinfo) )

setMethod("getFeatures", "gbRecord", function (x) x@features)

setMethod("getSequence", "gbRecord",  function (x) getSequence(x@seqinfo))


setMethod("ranges", "gbRecord",
          function (x, join = FALSE, key = TRUE, include = "none", exclude = "") {
            .make_GRanges(getFeatures(x), join = join, include = include,
                          exclude = exclude, key = key)
          })


setMethod("start", "gbRecord",
          function (x, join = FALSE, drop = TRUE) {
            start(getFeatures(x), join = join, drop = drop)
          })


setMethod("end", "gbRecord",
          function (x, join = FALSE, drop = TRUE) {
            end(getFeatures(x), join = join, drop = drop)
          })


setMethod("strand", "gbRecord",
          function (x, join = FALSE) {
            strand(getFeatures(x), join = join)
          })


setMethod("width", "gbRecord",
          function (x, join = FALSE) {
            width(getFeatures(x), join = join)
          })


# listers ----------------------------------------------------------------


setMethod("listQualif", "gbRecord", 
          function (x) {
            lapply(getFeatures(x), listQualif)
          })


# select -----------------------------------------------------------------


setMethod("select", "gbRecord",
          function (x, ..., keys = NULL, cols = NULL) {
            ans <- getFeatures(x)
            ans <- .select(ans, ..., keys = keys)
            ans <- .retrieve(ans, cols = cols)
            ans
          })


# shift ------------------------------------------------------------------


setMethod("shift", "gbRecord",
          function(x, shift, split=FALSE, order=FALSE, updateDb=FALSE)
            .shift_features(x=x, shift=shift, split=split, order=order))



# internal ---------------------------------------------------------------


setMethod('.dbSource', 'gbRecord', function (x) {
  parse_dbsource(getDBSource(x))
})

setMethod(".defline", "gbRecord", function (x) {
  paste0('gi|', getGeneID(x), .dbSource(x), getAccession(x), ' ', getDefinition(x))
})


