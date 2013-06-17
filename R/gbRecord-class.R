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
setClass("gbRecord",
         representation(
           locus = "character", moltype = "character", topology = "character",
           division = "character", update_date = "POSIXlt", create_date = "POSIXlt",
           version = "character", seqid = "character", dblink = "character",
           dbsource = "character", keywords = "character", source = "character",
           organism = "character", taxonomy = "character", references = "character",
           comment = "character", features = "gbFeatureList",
           contig = "gbLocationOrNull", seqinfo = "environment"),
         prototype(
           locus=NA_character_, moltype=NA_character_, topology=NA_character_,
           division=NA_character_, update_date=as.POSIXlt(NA), create_date=as.POSIXlt(NA),
           version=NA_character_, seqid=NA_character_, dblink=NA_character_,
           dbsource=NA_character_, keywords=NA_character_, source=NA_character_,
           organism=NA_character_, taxonomy=NA_character_, references=NA_character_,
           comment=NA_character_, seqinfo=new.env(parent=emptyenv())
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
    if (gb@rettype %ni% c("gb", "gp") || gb@retmode != "text")
      stop("Must use efetch with rettype='gbwithparts','gb', or 'gp' and retmode='text'")
    
    split_gb <- base::unlist(
      strsplit(Rentrez::content(gb, "text"), "\n\n")
    )
    n <- length(split_gb)
    parsed_gb_data <- vector("list", n)
    for (i in seq_len(n))
    {
      gb_data <- base::unlist(strsplit(split_gb[i], "\n"))
      parsed_gb_data[[i]] <- .parseGB( gb_data, with_sequence )
    }
  ##
  ## Parse random textConnections
  ##
  }
  else if (is(gb, "textConnection"))
  {
    con <- gb
    on.exit(close(con))
    parsed_gb_data <- list( .parseGB(readLines(con), with_sequence) )
  ##
  ## Parse GeneBank flat files  
  ##
  }
  else if (tryCatch(all(file.exists(gb)), error = function() FALSE))
  {
    n <- length(gb)
    parsed_gb_data <- vector("list", n)
    for (i in seq_len(n))
    {
      con <- file(gb[i], open="rt")
      parsed_gb_data[[i]] <- .parseGB( readLines(con), with_sequence )
      close(con)
    }
  }
  else
    stop(paste0("'gb' must be the path to a GenBank flat file or an 'efetch' ",
                "object containing GenBank records"))
  
  gbr_list <- list()
  for (gbk in parsed_gb_data)
  {
    gbr <- base::with(gbk, 
                      new("gbRecord",
                          locus = header[["locus"]],
                          moltype = header[["moltype"]],
                          topology = header[["topology"]],
                          division = header[["division"]],
                          update_date = header[["update_date"]],
                          create_date = header[["create_date"]],
                          version = header[["version"]],
                          seqid = header[["seqid"]],
                          dblink = header[["dblink"]],
                          dbsource = header[["dbsource"]],
                          keywords = header[["keywords"]],
                          source = header[["source"]],
                          organism = header[["organism"]],
                          taxonomy = header[["taxonomy"]],
                          references = header[["references"]],
                          comment = header[["comment"]],
                          features = features,
                          contig = contig,
                          seqinfo = seqenv)
    )
    gbr_list <- c(gbr_list, gbr)
  }
  
  if (length(gbr_list) == 1L) 
    gbr_list[[1L]]
  else
    gbRecordList(gbr_list)
}


# show -------------------------------------------------------------------


#' @importFrom XVector toString  subseq
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
                                        if (object@moltype == 'AA') 'aa' else 'bp',
                                        object@moltype, object@topology,
                                        object@division, object@update_date),
                                  offset=13, FORCE=TRUE)),
                sprintf("DEFINITION  %s\n",
                        linebreak(definition(object), offset=13, FORCE=TRUE)),
                sprintf("ACCESSION   %s\n", accession(object)),
                sprintf("VERSION     %s GI:%s\n", object@version, geneid(object)),
                sprintf("DBLINK      Project: %s\n", object@dblink),
                if (object@moltype == "AA") {
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
                        linebreak(object@taxonomy, offset=13, FORCE=TRUE)),
                sprintf("REFERENCE   %s\n", object@references),
                sprintf("COMMENT     %s\n",
                        linebreak(object@comment, offset=13, FORCE=TRUE)),
                if (!is.null(S)) {
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
                if (!is.null(object@contig)) {
                  sprintf("CONTIG      %s\n", linebreak(as(object@contig, "character"),
                                                        offset=13, split=",", FORCE=TRUE))
                })
            }
            
            cat(showme)
          })


# summary ----------------------------------------------------------------


#' @autoImports
setMethod("summary", "gbRecord",
          function (object, n=7, ...) {
            si <- seqinfo(object)
            acc <- seqnames(si)
            len <- unname(seqlengths(si))
            type <- if (object@type == 'AA') 'aa' else 'bp'
            def <- 
              ellipsize(obj=unname(genome(si)),
                        width=getOption("width") - base::nchar(len) - base::nchar(type) - 8)
            cat(sprintf("[[%s]]\n  %i %s: %s\n", acc, len, type, def), sep="")
            summary(object=features(object), n=n, setoff=2)
            
            return(invisible(NULL))
          })


# getters ----------------------------------------------------------------


setMethod("seqinfo", "gbRecord",
          function (x)  {
            tryCatch(get("seqinfo", x@seqinfo),
                     error = function (e) Seqinfo() )
          })


setMethod("seqlengths", "gbRecord",
          function (x) seqlengths(seqinfo(x)))


setMethod("accession", "gbRecord", 
          function (x) seqnames(seqinfo(x)))


setMethod("geneid", "gbRecord", 
          function (x, db = 'gi') {
            db.idx <- which(strsplitN(x@seqid, "|", 1, fixed = TRUE) == db)
            strsplitN(x@seqid, "|", 2, fixed = TRUE)[db.idx]
          })


setMethod("definition", "gbRecord", 
          function (x) genome(seqinfo(x)))


setMethod("features", "gbRecord", 
          function (x) x@features)


#' @autoImports
setMethod("sequence", "gbRecord", 
          function (x) {
            if (exists("sequence", envir=x@seqinfo))
              seq <- base::get("sequence", x@seqinfo)
            else {
              warning("No sequence associated with this record", call.=FALSE)
              seq <- BStringSet()
            }
            seq
          })


setMethod("ranges", "gbRecord",
          function (x, join = FALSE, key = TRUE, include = "none", exclude = "") {
            .make_GRanges(features(x), join = join, include = include,
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
          function (x, join = FALSE) {
            strand(features(x), join = join)
          })


setMethod("width", "gbRecord",
          function (x, join = FALSE) {
            width(features(x), join = join)
          })


# listers ----------------------------------------------------------------


setMethod("listQualif", "gbRecord", 
          function (x) {
            lapply(features(x), listQualif)
          })


# subsetting ----------------------------------------------------------


setMethod("[[", c("gbRecord", "character", "missing"),
          function(x, i) slot(object=x, name=i))


setMethod("$", "gbRecord",
          function(x, name) slot(object=x, name=name))


# select -----------------------------------------------------------------


setMethod("select", "gbRecord",
          function (x, ..., keys = NULL, cols = NULL) {
            ans <- features(x)
            ans <- .select(ans, ..., keys = keys)
            ans <- .retrieve(ans, cols = cols)
            ans
          })


# shift ------------------------------------------------------------------


setMethod("shift", "gbRecord",
          function(x, shift, split=FALSE, order=FALSE, updateDb=FALSE)
            .shift_features(x=x, shift=shift, split=split, order=order))

