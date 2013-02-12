#' validate a gbRecord database which if referenced by gbFeature or
#' gbFeatureList objects by their '.Dir' slot. 
#' @keywords internal
hasValidDb <- function (object, verbose=TRUE) {
  
  if (!.hasSlot(object, "db") && !is(object@db, "gbRecord")) {
    if (verbose) {
      message("Object has no 'db' slot")
    }
    return(FALSE)
  }
  
  if (!file.exists(object@db@dir)) {
    if (verbose)
      message(sprintf("Directory %s does not exist.", sQuote(object@db@name)))
    return(FALSE)
  }
  
  if (any(idx <- is.na(match(.GBFIELDS, dir(object@db@dir))))) {
    if (verbose)
      message(sprintf("Field(s) %s are missing from database.",
                      sQuote(paste(.GBFIELDS[-idx], collapse=","))))
    return(FALSE)
  }
  
  TRUE
}


#' validate a gbRecord database (i.e. check if the db directory contains
#' all fields).
#' @keywords internal
isValidDb <- function (object, verbose=TRUE) {
  idx <- is.na(charmatch(.GBFIELDS, dir(object@dir)))
  if (any(idx)) {
    if (verbose) {
      message(sprintf("%s are missing from database.",
                      sQuote(paste(.GBFIELDS[idx], collapse=", "))))
    }
    return(FALSE)
  } 
  return(TRUE)
}

# check if the db directory has been moved from the location
# where it was instantiated
hasNewPath <- function (x) {
  !identical(x@dir, x$features@.Info@db@dir) 
}

# if yes update the gbFeatureList@.Dir and gbFeature@.Dir
# slots to  point to the current directory.
#' @autoImports
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


is.gbRecordDb <- function (gb) {
  if (tryCatch(file.exists(gb), error = function(e) FALSE) &&
      file.info(gb)[["isdir"]]) {
    all(not.na(charmatch(.GBFIELDS, dir(gb))))
  } else {
    FALSE
  }
}

