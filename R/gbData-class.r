##' The main class for GenBank records
##'
##' All data parsed from GenBank flat file records is stored in instances
##' of the gbData class. It is implemented as an extension of a
##' '\link[filehash]{filehashRDS-class}'. Thus, all data is stored as a
##' file-based hash table.
##'
##' @name gbData-class
##' @rdname gbData-class
##' @exportClass gbData
##' 
##' @examples
##' getSlots("gbData")
setClass("gbData", contains="filehashRDS")


##' Create and/or initialise a \code{gbData} filehash database.
##'
##' @usage initGenBank(db_name, create=FALSE, ...)
##'
##' @param db_name name of the database or database object.
##' @param create if \code{TRUE} a new database is created and initialised,
##' else a database object is initialised from an existing database.
##' @param ... other arguments passed to methods
##'
##' @return Returns a '\code{gbData}' object inheriting from class 
##' \code{\link[filehash]{filehashRDS-class}} 
##'
##' @name initGenBank
##' @aliases initGenBank,gbData-method
##' @docType methods
##' @rdname initGenBank-methods
##' @export
##'
##' @examples
##' ##
setGeneric("initGenBank", function(db_name, ...) 
  standardGeneric("initGenBank"))

setMethod("initGenBank", "ANY",
          function(db_name, create=FALSE, ...) {
            db_name <- sub("/$", "", db_name, perl=TRUE)
            if (create) 
              dbCreate(db_name, "RDS")
            
            new("gbData", 
                dir=normalizePath(db_name),
                name=basename(db_name))
          })

## Accesssor methods
setGeneric("features", function(object, ...) standardGeneric("features"))
setGeneric("sequence", function(object, ...) standardGeneric("sequence"))

setMethod("features", "gbData", 
          function (object)  {
            dbFetch(object, "features")
          })

setMethod("sequence", "gbData", 
          function (object) {
            if (length(dbFetch(object, "sequence")) < 2L)
              return(dbFetch(object, "sequence")[[1L]])
            else
              return(dbFetch(object, "sequence"))
          })

#### Constructor
gbData <- function (db_dir, header, features, sequence=NULL) {

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


#### show method
setMethod("show", "gbData",
          function(object) {
            if(length(object@name) == 0)
              stop("database does not have a name")
            cat(sprintf("'%s' database '%s' with %i features\n", 
                        as.character(class(object)), object@name,
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
                                          start=length(object$sequence) - 
                                            getOption("width")+15,
                                          end=length(object$sequence))))
                })
          })

#### Subsetting
setMethod("[[", signature(x = "gbData", i = "character", j = "missing"),
          function(x, i, j) {
            if (i == "features") {
              return(features(x))
            } else if (i == "sequence") {
              return(sequence(x))
            } else {
              dbFetch(x, i) }
          })

setMethod("$", signature(x = "gbData"),
          function(x, name) {
            if (name == "sequence") {
              return(sequence(x)) }
            else if (name == "features") {
              return(features(x)) }
            else {
              dbFetch(x, name)
            }
          })

## Select method
setGeneric("select", function(db, key="ANY", qualifier="ANY", location="ANY")
  standardGeneric("select"))

setMethod("select", signature(db="gbData"), 
          definition=function (db, key, qualifier, location)  {
            .select(x=dbFetch(db, "features"), key, qualifier, location) 
          })
