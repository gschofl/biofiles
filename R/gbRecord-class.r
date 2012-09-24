
# gbRecord-Class ------------------------------------------------------

#' @include gbFeatureList-class.r
NULL

#' gbRecord class
#' 
#' gbRecord is an S4 class that provides a container for data parsed from
#' a GenBank flat file record. It is implemented as an extension of a
#' \dQuote{\code{\linkS4class{filehashRDS-class}}} with no additional slots.
#'
#' @rdname gbRecord
#' @export
#' @classHierarchy
#' @classMethods
.gbRecord <- setClass("gbRecord", contains="filehashRDS")


setValidity("gbRecord", function (object) {
  # at the moment do nothing but the default checks
  TRUE
})


# constructor ------------------------------------------------------------


#' @keywords internal
setMethod("initialize", "gbRecord",
          function (.Object, dir=character(0), name=character(0), verbose=TRUE) {
            if (missing(dir))
              stop("No database directory provided")
            if (missing(name))
              stop("No database name provided")
            
            .Object@dir <- dir
            .Object@name <- name
            
            if (isValidDb(.Object, verbose=FALSE)) {
              if (verbose) message("Intializing gbRecord")
              
              if (hasNewPath(.Object)) {
                if (verbose) message("Updating db directory location ...")
                updateDirectory(.Object)
              }
            }
            
            .Object
          })


#' @keywords internal
gbRecord <- function (db_dir, header, features, sequence=NULL) {
  
  db <- initGB(db_dir, create=TRUE)
  
  if (is.null(header$accession))
    stop("No accession number available")
  if (is.null(header$definition))
    stop("No sequence definition available")
  if (!is(features, "gbFeatureList"))
    stop("Features must be a 'gbFeatureList' instance")
  
  dbInsert(db, "locus", header[["locus"]])
  dbInsert(db, "length", header[["length"]])
  dbInsert(db, "type", header[["type"]])
  dbInsert(db, "topology", header[["topology"]])
  dbInsert(db, "division", header[["division"]])
  dbInsert(db, "date", header[["date"]])
  dbInsert(db, "definition", header[["definition"]]) # mandatory
  dbInsert(db, "accession", header[["accession"]])   # mandatory
  dbInsert(db, "version", header[["version"]])
  dbInsert(db, "GI", header[["GI"]])
  dbInsert(db, "dblink", header[["dblink"]])
  dbInsert(db, "dbsource", header[["dbsource"]])
  dbInsert(db, "keywords", header[["keywords"]])
  dbInsert(db, "source", header[["source"]])
  dbInsert(db, "organism", header[["organism"]])
  dbInsert(db, "lineage", header[["lineage"]])
  dbInsert(db, "references", header[["references"]])
  dbInsert(db, "comment", header[["comment"]])
  dbInsert(db, "features", features) # mandatory
  dbInsert(db, "sequence", sequence)

  return(invisible(db))
}


# show -------------------------------------------------------------------


setMethod("show", "gbRecord",
          function (object) {
            if(length(object@name) == 0)
              stop("database does not have a name")
            
            Seq <- dbFetch(object, "sequence")
            W <- getOption("width")
            showme <- paste0(
              sprintf("%s database %s with %i features\n", 
                      sQuote(class(object)), sQuote(object@name),
                      length(dbFetch(object, "features"))),
              
              sprintf("LOCUS       %s\n",
                      linebreak(paste(dbFetch(object, "locus"),
                                      dbFetch(object, "length"),
                                      names(dbFetch(object, "length")),
                                      dbFetch(object, "type"),
                                      dbFetch(object, "topology"),
                                      dbFetch(object, "division"),
                                      dbFetch(object, "date")),
                                offset=13, FORCE=TRUE)),
              
              sprintf("DEFINITION  %s\n",
                      linebreak(dbFetch(object, "definition"), offset=13,
                                FORCE=TRUE)),
              
              sprintf("ACCESSION   %s\n", dbFetch(object, "accession")),
              
              sprintf("VERSION     %s GI:%s\n",
                      dbFetch(object, "version"), dbFetch(object, "GI")),
              
              sprintf("DBLINK      Project: %s\n", dbFetch(object, "dblink")),
              
              if (dbFetch(object, "type") == "AA") {
                sprintf("DBSOURCE    %s\n",
                        linebreak(dbFetch(object, "dbsource"), offset=13,
                                  FORCE=TRUE))
              },
              
              sprintf("KEYWORDS    %s\n",
                      linebreak(dbFetch(object, "keywords"), offset=13,
                                FORCE=TRUE)),
              
              sprintf("SOURCE      %s\n",
                      linebreak(dbFetch(object, "source"), offset=13,
                                FORCE=TRUE)),
              
              sprintf("  ORGANISM  %s\n",
                      linebreak(dbFetch(object, "organism"), offset=13,
                                FORCE=TRUE)),
              
              sprintf("            %s\n",
                      linebreak(dbFetch(object, "lineage"), offset=13,
                                FORCE=TRUE)),
              
              sprintf("REFERENCE   %s\n", dbFetch(object, "references")),
              
              sprintf("COMMENT     %s\n",
                      linebreak(dbFetch(object, "comment"), offset=13,
                                FORCE=TRUE)),
              
              if (not.null(Seq)) {
                if (Seq@ranges@width[1L] < W - 14)
                  sprintf("ORIGIN      %s\n", toString(Seq))
                else
                  sprintf("ORIGIN      %s\n             ...\n             %s\n",
                          toString(subseq(Seq, start=1, end=W - 14)),
                          toString(subseq(Seq,
                                          start=length(Seq[[1L]]) -  W + 15,
                                          end=length(Seq[[1L]]))))
              })
            
            cat(showme)
          })


# initGB -----------------------------------------------------------------


setMethod("initGB", "character",
          function(db_dir, create = FALSE, ...) {
            db_dir <- sub("/$", "", db_dir, perl=TRUE)
            if (create) {
              dbCreate(db_dir, "RDS")
            }
            if(!file.exists(db_dir)) {
              stop(sprintf("Database directory %s does not exist",
                           sQuote(db_dir)))
            }
            .gbRecord(dir=normalizePath(db_dir), name=basename(db_dir), ...)
          })


# getters ----------------------------------------------------------------


setMethod("features", "gbRecord", 
          function (x) dbFetch(x, "features"))


setMethod("sequence", "gbRecord", 
          function (x) dbFetch(x, "sequence"))


setMethod("accession", "gbRecord", 
          function (x) dbFetch(x, "accession"))


setMethod("definition", "gbRecord", 
          function (x) dbFetch(x, "definition"))



# subsetting ----------------------------------------------------------


#' @export
setMethod("[[", c("gbRecord", "character", "missing"),
          function(x, i) dbFetch(x, i))


#' @export
setMethod("$", "gbRecord",
          function(x, name) dbFetch(x, name))


# select -----------------------------------------------------------------


setMethod("select", "gbRecord",
          function (x, ..., keys = NULL, cols = NULL) {
            ans <- dbFetch(x, "features")
            ans <- .select(ans, ..., keys = keys)
            ans <- .retrieve(ans, cols = cols)
            ans
          })



# write ------------------------------------------------------------------


setMethod("write", "gbRecord",
          function (x, file = "data") {
            if (file.exists(x@dir)) {
              if (file == "data") {
                file <- file.path(getwd(), paste0(dbFetch(x, "accession"), ".db"))
              }
              dir.create(path=file)
              right <- file.copy(from=Sys.glob(file.path(x@dir, "*")),
                                 to=file, recursive=TRUE)
              
              if (all(right)) {
                cat(paste("Record written to", sQuote(file)))
              }
              
              return(invisible(file))
              
            } else {
              stop(paste("Database file", sQuote(x@dir), "does not exist"))
            }
          })


# shift ------------------------------------------------------------------


setMethod("shift", "gbRecord",
          function(x, shift, split=FALSE, order=FALSE, updateDb=FALSE)
            .shift_features(x=x, shift=shift, split=split, order=order, 
                            updateDb=updateDb))


# revcomp ----------------------------------------------------------------


setMethod("revcomp", "gbRecord",
          function(x, order=FALSE, updateDb=FALSE)
            .revcomp_features(x=x, order=order,
                              updateDb=updateDb))
