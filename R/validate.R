#' @keywords internal
hasValidDb <- function (object, verbose=TRUE) {
  
  if (is(object, "gbFeature") || is(object, "gbFeatureList")) {
    db <- slot(slot(seqinfo(object), "db"), "dir")
  } else if (!isS4(object)) {
    attr <- attributes(object)
    if (!all(c("accession", "definition", "dir") %in% names(attr))) {
      return(FALSE)
    }
    db <- attr$dir
  }

  if (!tryCatch(file.exists(db), error=function(e) TRUE)) {
    if (verbose)
      message(sprintf("Directory %s does not exist.", sQuote(db)))
    return(FALSE)
  }
  
  isValidDb(db, verbose = verbose)
}


#' @keywords internal
isValidDb <- function (db, verbose=TRUE) {
  idx <- is.na(charmatch(.GBFIELDS, dir(db)))
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
hasNewPath <- function (object) {
  !identical(object@dir, seqinfo(features(object))@db@dir)
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

