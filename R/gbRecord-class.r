##' A GenBank record
##'
##' Data parsed from a GenBank flat file record are stored in instances
##' of the gbRecord class. It is implemented as an extension of a
##' '\link[filehash]{filehashRDS-class}' with no additional slots.
##'
##' @name gbRecord-class
##' @rdname gbRecord-class
##' @exportClass gbRecord
##' 
##' @examples
##' getSlots("gbRecord")
setClass("gbRecord", contains="filehashRDS")

##' Create and/or initialise a gbRecord filehash database
##'
##' @usage initGenBank(db_name, create=FALSE, ...)
##'
##' @param db_name name of the database or database object.
##' @param create if \code{TRUE} a new database is created and initialised,
##' else a database object is initialised from an existing database.
##' @param ... other arguments passed to methods
##'
##' @return Returns a '\code{gbRecord}' object inheriting from class 
##' \code{\link[filehash]{filehashRDS-class}} 
##'
##' @export
##' @docType methods
##' @rdname initGenBank-methods
##'
##' @examples
##' ##
setGeneric("initGenBank", function(db_dir, ...) 
  standardGeneric("initGenBank"))

##' @aliases initGenBank,gbRecord,gbRecord-method
##' @rdname initGenBank-methods
setMethod("initGenBank",
          signature(db_dir="ANY"),
          function(db_dir, create=FALSE, ...)
          {
            db_dir <- sub("/$", "", db_dir, perl=TRUE)
            if (create) 
              dbCreate(db_dir, "RDS")
            
            if(!file.exists(db_dir))
              stop(sprintf("database directory %s does not exist", sQuote(db_dir)))

            new("gbRecord", 
                dir=normalizePath(db_dir),
                name=basename(db_dir))
          })


setMethod("initialize",
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
              if (verbose) message("Reintializing gbRecord")
              
              if (hasNewPath(.Object)) {
                if (verbose) message("Updating db directory location ...")
                updateDirectory(.Object)
              }
            }
            
            .Object
          })

### Accessor methods
##' Retrieve features of a GenBank record.
##'
##' @usage features(data)
##'
##' @param data An instance of \code{\link{gbRecord-class}}.
##'
##' @return The \code{gbFeatureList}-class object contained in a gbRecord
##' database.
##'
##' @export
##' @docType methods
##' @rdname gbRecord-class
##'
##' @examples
##' ##
setGeneric("getFeatures", function(data, ...) standardGeneric("getFeatures"))

##' @aliases features,gbRecord,gbRecord-method
##' @rdname gbRecord-class
setMethod("getFeatures",
          signature(data="gbRecord"), 
          function (data) {
            dbFetch(data, "features")
          })


#### Constructor
gbRecord <- function (db_dir, header, features, sequence=NULL) {

  db <- initGenBank(db_dir, create=TRUE)

  if (is.null(header$accession))
    stop("No accession number available")
  if (is.null(header$definition))
    stop("No sequence definition available")
  if (!is(features, "gbFeatureList"))
    stop("Features must be a 'gbFeatureList")
  
  .insertData(db=db, header=header, features=features, sequence=sequence)
  
  db
}

.insertData <- function (db, header, features, sequence) {
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

### show method
##' @export
##' @aliases show,gbRecord,gbRecord-method
##' @rdname gbRecord-class
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


##' Subsetting
##'
##' @export
##' @aliases [[,gbRecord,gbRecord-method
##' @rdname extract-methods
setMethod("[[",
          signature(x="gbRecord", i="character", j="missing"),
          function(x, i) {
            dbFetch(x, i)
          })

##' @export
##' @aliases $,gbRecord,gbRecord-method
##' @rdname extract-methods
setMethod("$",
          signature(x="gbRecord"),
          function(x, name) {
            dbFetch(x, name)
          })

    
### Select method ##########################################################
##' @export
##' @aliases select,gbRecord-method
##' @rdname select-methods
setMethod("select",
          signature(x="gbRecord"),
          function (x, what=c(""), which=list()) {
            x <- dbFetch(db=x, key="features")
            ans <- .select(x, which=which)
            ans <- .retrieve(x=ans, what=what)
            ans
          })
