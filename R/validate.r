GBFIELDS <- paste0("@G@I|accession|comment|date|dblink|dbsource|",
                   "definition|division|features|keywords|length|",
                   "lineage|locus|organism|references|sequence|",
                   "source|topology|type|version")

##' validate a gbRecord database which if referenced by gbFeature or
##' gbFeatureList objects by their '.Dir' slot. 
##' @keywords internal
hasValidDb <- function (x, verbose=TRUE, fields=GBFIELDS) {
  
  if (!.hasSlot(x, ".Dir")) {
    if (verbose) message("Object has no '.Dir' slot")
    return( FALSE )
  }
  
  if (!file.exists(x@.Dir)) {
    if (verbose) 
      message(gettextf("directory %s does not exist", sQuote(x@.Dir)))
    return( FALSE )
  }
  
  if (sum(idx <- grepl(fields, list.files(x@.Dir))) != 20L) {
    if (verbose)
      message(paste("Field(s) are missing from the database:",
                    paste(strsplit(fields, split="\\|")[[1]][!idx], collapse=",")))
    return( FALSE )
  }
  
  TRUE
}

##' validate a gbRecord database (i.e. check if the db directory contains
##' all fields).
##' @keywords internal
isValidDb <- function (x, verbose=TRUE, fields=GBFIELDS) {
  
  if (sum(idx <- grepl(fields, list.files(x@dir))) != 20L) {
    if (verbose)
      message(paste("Field(s) are missing from the database:",
                    paste(strsplit(fields, split="\\|")[[1]][!idx], collapse=",")))
    return( FALSE )
  }
  
  TRUE
}

## check if the db directory has been moved from the location
## where it was instantiated
hasNewPath <- function (x) {
  ## check if the db directory has moved
  !identical(x@dir, x$features@.Dir) 
}

## if yes update the gbFeatureList@.Dir and gbFeature@.Dir
## slots to  point to the current directory.
updateDirectory <- function (db) {
  newPath <- db@dir
  data <- dbFetch(db, "features")
  data <- 
    `slot<-`(data, name=".Dir", check=FALSE, value=newPath)
  data <- 
    `slot<-`(data, name=".Data", check=FALSE,
             value=lapply(data@.Data, `slot<-`, name=".Dir",
                          check=FALSE, value=newPath))
  dbDelete(db, "features")
  dbInsert(db, "features", data)
}



