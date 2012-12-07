
# gbRecord-Class ------------------------------------------------------

#' @include gbFeatureList-class.r
NULL

#' gbRecord
#' 
#' \dQuote{gbRecord} is an S4 class that provides a container for data
#' parsed from a GenBank record. It is implemented as a
#' \code{\linkS4class{filehashRDS}} database.
#'
#' @rdname gbRecord
#' @export
#' @classHierarchy
#' @classMethods
setClass("gbRecord", contains="filehashRDS")


setValidity("gbRecord", function (object) {
  # at the moment do nothing but the default checks
  TRUE
})


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


# constructor ------------------------------------------------------------


#' \code{gbRecord} instances can be construced by parsing a GenBank flat
#' file or an \code{\linkS4class{efetch}} instance containing one or more
#' GenBank records.
#' If \code{gb} points to a valid \code{gbRecord} database, a \code{gbRecord}
#' object is initialised in the global environment.  
#'
#' @details
#' For a description of the GenBank format see
#' \url{http://www.ncbi.nlm.nih.gov/collab/FT/}
#'
#' @param gb Path to a valid \code{gbRecord} database, a GenBank flat file or
#' an \code{\linkS4class{efetch}} object containing GenBank record(s).
#' @param with_sequence Parse with sequence information if avaliable.
#' @param force Overwrite existing database directories without prompting.
#' @return A (list of) \code{\linkS4class{gbRecord}} object(s).
#' @export
#' @autoImports
gbRecord <- function (gb, with_sequence = TRUE, force = FALSE) {
  
  if (is_gbRecord_db(gb))
    return(init_db(gb, create = FALSE))
  
  if (is(gb, "efetch")) {
    if (gb@rettype %ni% c("gb", "gp") || gb@retmode != "text")
      stop("Must use efetch with rettype='gbwithparts','gb', or 'gp' and retmode='text'")
    
    split_gb <- unlist(strsplit(content(gb, "text"), "\n\n"))
    n <- length(split_gb)
    db_path <- replicate(n, tempfile(fileext=".db"))
    parsed_data <- vector("list", n)
    for (i in seq_len(n)) {
      gb_data <- unlist(strsplit(split_gb[i], "\n"))
      parsed_data[[i]] <- .parseGB(gb_data, db_path[i],
                                   with_sequence=with_sequence,
                                   force=force)
    }
  }
  else if (!isS4(gb) && (is(gb, "connection") || file.exists(gb))) {
    if (is(gb, "textConnection")) {
      con <- gb
      db_path <- tempfile(fileext=".db")
    } else {
      con <- file(gb, open="rt")
      db_path <- paste0(gb, ".db")
    }
    on.exit(close(con))
    parsed_data <- list(.parseGB(gb_data=readLines(con), db_path, with_sequence, force))
    
  } else {
    stop("'gb' must be a valid GenBank flat file or an 'efetch' object containing GenBank records")
  }
  
  gbk_list <- list()
  accn <- character()
  for (gbk in parsed_data) {
    db <- init_db(gbk[["db_dir"]], create = TRUE)
    
    if (is.null(gbk[["header"]][["accession"]]))
      stop("No accession number available")
    if (is.null(gbk[["header"]][["definition"]]))
      stop("No sequence definition available")
    if (!is(gbk[["features"]], "gbFeatureList"))
      stop("Features must be a 'gbFeatureList' instance")
    
    dbInsert(db, "locus", gbk[["header"]][["locus"]])
    dbInsert(db, "length", gbk[["header"]][["length"]])
    dbInsert(db, "type", gbk[["header"]][["type"]])
    dbInsert(db, "topology", gbk[["header"]][["topology"]])
    dbInsert(db, "division", gbk[["header"]][["division"]])
    dbInsert(db, "date", gbk[["header"]][["date"]])
    dbInsert(db, "definition", gbk[["header"]][["definition"]]) # mandatory
    dbInsert(db, "accession", gbk[["header"]][["accession"]])   # mandatory
    dbInsert(db, "version", gbk[["header"]][["version"]])
    dbInsert(db, "GI", gbk[["header"]][["GI"]])
    dbInsert(db, "dblink", gbk[["header"]][["dblink"]])
    dbInsert(db, "dbsource", gbk[["header"]][["dbsource"]])
    dbInsert(db, "keywords", gbk[["header"]][["keywords"]])
    dbInsert(db, "source", gbk[["header"]][["source"]])
    dbInsert(db, "organism", gbk[["header"]][["organism"]])
    dbInsert(db, "lineage", gbk[["header"]][["lineage"]])
    dbInsert(db, "references", gbk[["header"]][["references"]])
    dbInsert(db, "comment", gbk[["header"]][["comment"]])
    dbInsert(db, "features", gbk[["features"]]) # mandatory
    dbInsert(db, "sequence", gbk[["sequence"]])
    
    gbk_list <- c(gbk_list, db)
    accn <- c(accn, gbk[["header"]][["accession"]])
  }
  
  if (length(gbk_list) == 1L) 
    gbk_list <- gbk_list[[1L]]
  else
    names(gbk_list) <- accn
 
  gbk_list
}



#' @keywords internal
#' @autoImports
init_db <- function(db_dir, create = FALSE, ...) {
  db_dir <- sub("/$", "", db_dir, perl=TRUE)
  if (create) {
    dbCreate(db_dir, "RDS")
  }
  if (!file.exists(db_dir)) {
    stop(sprintf("Database directory %s does not exist",
                 sQuote(db_dir)))
  }
  new("gbRecord", dir=normalizePath(db_dir), name=basename(db_dir), ...)
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


setMethod("[[", c("gbRecord", "character", "missing"),
          function(x, i) dbFetch(x, i))


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

