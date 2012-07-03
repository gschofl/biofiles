
# gbRecord-Class ------------------------------------------------------

##' @include gbFeatureList-class.r
NULL

##' gbRecord class
##' 
##' gbRecord is an S4 class that provides a container for data parsed from
##' a GenBank flat file record. It is implemented as an extension of a
##' '\link[filehash]{filehashRDS-class}' with no additional slots.
##'
##' @name gbRecord-class
##' @rdname gbRecord-class
##' @exportClass gbRecord
##' @aliases initGB,initGB-method,gbRecord-method
##' @aliases show,gbRecord-method
##' @aliases getFeatures,gbRecord-method
##' @aliases [[,gbRecord-method
##' @aliases $,gbRecord-method
##' @aliases select,select-method,gbRecord-method
##' @aliases write,gbRecord-method
##' @aliases shiftFeatures,gbFeature-method
.gbRecord <-
  setClass("gbRecord", contains="filehashRDS")


# show-method ---------------------------------------------------------


##' @export
setMethod("show",
          signature(object="gbRecord"),
          function(object) {
            if(length(object@name) == 0)
              stop("database does not have a name")
            cat(sprintf("%s database %s with %i features\n", 
                        sQuote(class(object)), sQuote(object@name),
                        length(object$features)),
                sprintf("LOCUS       %s\n", linebreak(paste(object$locus,
                                                            object$length, names(object$length),
                                                            object$type, object$topology, object$division,
                                                            object$date), offset=13, FORCE=TRUE)),
                sprintf("DEFINITION  %s\n", linebreak(object$definition,
                                                      offset=13, FORCE=TRUE)),
                sprintf("ACCESSION   %s\n", object$accession),
                sprintf("VERSION     %s GI:%s\n", object$version, object$GI),
                sprintf("DBLINK      Project: %s\n", object$dblink),
                if (object$type == "AA") {
                  sprintf("DBSOURCE    %s\n", linebreak(object$dbsource,
                                                        offset=13, FORCE=TRUE))
                },
                sprintf("KEYWORDS    %s\n", linebreak(object$keywords,
                                                      offset=13, FORCE=TRUE)),
                sprintf("SOURCE      %s\n", linebreak(object$source,
                                                      offset=13, FORCE=TRUE)),
                sprintf("  ORGANISM  %s\n", linebreak(object$organism,
                                                      offset=13, FORCE=TRUE)),
                sprintf("            %s\n", linebreak(object$lineage,
                                                      offset=13, FORCE=TRUE)),
                sprintf("REFERENCE   %s\n", object$references),
                sprintf("COMMENT     %s\n", linebreak(object$comment,
                                                      offset=13, FORCE=TRUE)),
                if (!is.null(object$sequence)) {
                  sprintf("ORIGIN      %s\n             ...\n             %s\n",
                          toString(subseq(object$sequence, 
                                          start=1, end=getOption("width")-14)),
                          toString(subseq(object$sequence,
                                          start=length(object$sequence[[1L]]) - 
                                            getOption("width")+15,
                                          end=length(object$sequence[[1L]]))))
                })
          })

# Constructor ---------------------------------------------------------

##' Create and/or initialise a gbRecord filehash database
##'
##' @usage initGB(db_name, create=FALSE, ...)
##'
##' @param db_name name of the database or database object.
##' @param create if \code{TRUE} a new database is created and initialised,
##' else a database object is initialised from an existing database.
##' @param ... other arguments passed to methods
##'
##' @return Returns a \sQuote{\code{gbRecord-class}} object inheriting from 
##' \code{\link[filehash]{filehashRDS-class}} 
##'
##' @export
##' @docType methods
setGeneric("initGB",
           #### initGB-generic ####
           function(db_dir, ...) {
             standardGeneric("initGB")
           })

##' @export
setMethod("initGB",
          #### initGB-method ####
          signature(db_dir="ANY"),
          function(db_dir, create=FALSE, ...) {
            db_dir <- sub("/$", "", db_dir, perl=TRUE)
            if (create) 
              dbCreate(db_dir, "RDS")
            
            if(!file.exists(db_dir))
              stop(sprintf("database directory %s does not exist", sQuote(db_dir)))
            
            .gbRecord(dir=normalizePath(db_dir), name=basename(db_dir), ...)
          })

#' @keywords internal
setMethod("initialize",
          #### initialize-method ####
          signature(.Object="gbRecord"),
          function (.Object, dir=character(0), name=character(0), verbose=TRUE)
          {
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
gbRecord <- function (db_dir, header, features, sequence=NULL)
{
  
  db <- initGB(db_dir, create=TRUE)
  
  if (is.null(header$accession))
    stop("No accession number available")
  if (is.null(header$definition))
    stop("No sequence definition available")
  if (!is(features, "gbFeatureList"))
    stop("Features must be a 'gbFeatureList")
  
  .insertData(db=db, header=header, features=features, sequence=sequence)
  
  db
}

#' @keywords internal
.insertData <- function (db, header, features, sequence)
{
  dbInsert(db, "locus", header$locus)
  dbInsert(db, "length", header$length)
  dbInsert(db, "type", header$type)
  dbInsert(db, "topology", header$topology)
  dbInsert(db, "division", header$division)
  dbInsert(db, "date", header$date)
  dbInsert(db, "definition", header$definition)  ## mandatory
  dbInsert(db, "accession", header$accession)   ## mandatory
  dbInsert(db, "version", header$version)
  dbInsert(db, "GI", header$GI)
  dbInsert(db, "dblink", header$dblink)
  dbInsert(db, "dbsource", header$dbsource)
  dbInsert(db, "keywords", header$keywords)
  dbInsert(db, "source", header$source)
  dbInsert(db, "organism", header$organism)
  dbInsert(db, "lineage", header$lineage)
  dbInsert(db, "references", header$references)
  dbInsert(db, "comment", header$comment)
  dbInsert(db, "features", features)    ## mandatory
  dbInsert(db, "sequence", sequence)
}


# Accessors -----------------------------------------------------------

##' Retrieve features of a GenBank record.
##'
##' @usage getFeatures(x, ...)
##'
##' @param data An instance of \code{\link{gbRecord-class}}.
##'
##' @family accessors
##' @return The \code{\link{gbFeatureList-class}} object contained in a
##' gbRecord database.
##'
##' @docType methods
##' @export
setGeneric("getFeatures", function(x, ...) standardGeneric("getFeatures"))

##' @export
setMethod("getFeatures", "gbRecord", 
          function (x) dbFetch(x, "features")
          )


# Subsetting ----------------------------------------------------------


##' @export
setMethod("[[",
          signature(x="gbRecord", i="character", j="missing"),
          function(x, i) dbFetch(x, i)
          )

##' @export
setMethod("$", "gbRecord",
          function(x, name) dbFetch(x, name)
          )

##' @export
setMethod("select", "gbRecord",
          function (x, subset = "", select = "") {
            ans <- dbFetch(db=x, key="features")
            ans <- .select(x=ans, which=subset)
            ans <- .retrieve(x=ans, which=select)
            ans
          })

##' @export
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



# Shift features ------------------------------------------------------


#' @export
setMethod("shift", "gbRecord",
          function(x, shift, split=FALSE, order=FALSE, update_db=FALSE) {
            .shift_features(x=x, shift=shift, split=split, order=order, 
                            update_db=update_db)
          })

#' @export
setMethod("revcomp", "gbRecord",
          function(x, order=FALSE, update_db=FALSE) {
            .revcomp_features(x=x, order=order,
                              update_db=update_db)
          })
