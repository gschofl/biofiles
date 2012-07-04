
# Generics ------------------------------------------------------------

setGeneric("start", function (x, ...) standardGeneric("start"))

setGeneric("start<-", function (x, ...) standardGeneric("start<-"))

setGeneric("end", function (x, ...) standardGeneric("end"))

setGeneric("end<-", function (x, ...) standardGeneric("end<-"))

setGeneric("width", function (x, ...) standardGeneric("width"))

setGeneric("strand", function (x, ...) standardGeneric("strand"))

setGeneric("strand<-", function (x, ...) standardGeneric("strand<-"))

setGeneric("partial", function (x, ...)  standardGeneric("partial"))


## getter-generics ========================================================== 


##' get genomic location of a GenBank feature
##'
##' @usage getLocation(x, attributes=TRUE, join=FALSE)
##'
##' @param x A \code{\link{gbFeature}} or \code{\link{gbFeatureList}} object
##' @param attributes set the \code{accession}, \code{definition},
##' \code{database} attributes.
##' @param join combine compound locations
##'
##' @return A data frame
##'
##' @docType methods
##' @export
setGeneric("getLocation",
           function (x, attributes = FALSE, join = FALSE, ...) 
             standardGeneric("getLocation")
           )


##' get index of a GenBank feature.
##'
##' @usage getIndex(x)
##'
##' @param x A gbFeature or gbFeatureList object
##'
##' @return A numeric vector of feature indeces
##'
##' @docType methods
##' @export
setGeneric("getIndex",
           function (x, attributes = FALSE, ...) 
             standardGeneric("getIndex")
           )


##' @docType methods
##' @export
setGeneric("getKey",
           function(x, attributes = FALSE, ...) 
             standardGeneric("getKey")
           )


##' @docType methods
##' @export
setGeneric("getQualifier",
           function(x, which = "", attributes = FALSE, ...)
             standardGeneric("getQualifier")
           )


##' @docType methods
##' @export
setGeneric("dbXref",
           function(x, db = NULL, ...) 
             standardGeneric("dbXref")
           )


##' @docType methods
##' @export
setGeneric("getSequence",
           function(x, ...) 
             standardGeneric("getSequence")
           )


##' @docType methods
##' @export
setGeneric("hasKey",
           function(x, key, ...)
             standardGeneric("hasKey")
           )


##' @docType methods
##' @export
setGeneric("hasQualifier",
           function(x, qualifier, ...)
             standardGeneric("hasQualifier")
           )


##' Retrieve features of a GenBank record.
##'
##' @usage getFeatures(x, ...)
##'
##' @param data An instance of \code{\link{gbRecord-class}}.
##'
##' @return The \code{\link{gbFeatureList-class}} object contained in a
##' gbRecord database.
##'
##' @docType methods
##' @export
setGeneric("getFeatures",
           function (x, ...)
             standardGeneric("getFeatures")
           )


## shift-generic ============================================================ 
##
##' shift location of features in a GenBank record
##'
##' @usage shift(x, shift=0L, split=FALSE, order=FALSE, update_db=FALSE)
##'
##' @param x A gbLocation, gbFeature, gbFeatureList, or gbRecord object
##' (gbFeatureLists must include a 'source' field).
##' @param shift Number of basepairs (or aa residues) to shift.
##' @param split (For gbFeatureList and gbRecord objects) Should a feature
##' that spans across the end of the sequence be split.
##' @param order (For gbFeatureList and gbRecord objects) Should the
##' resulting gbFeatureList be reordered.
##' @param update_db Should filehash database be updated with new feature
##' locations.
##'
##' @return A gbLocation, gbFeature, or gbFeatureList object
##'
##' @docType methods
##' @export
setGeneric("shift", function (x, shift = 0L, ...) standardGeneric("shift"))


## select-generic =========================================================== 
##
##' Select elements from a GenBank Record
##' 
##' @usage select(x, subset = c(""), select =c(""))
##' 
##' @param x A \sQuote{\code{gbRecord}} or \sQuote{\code{gbFeatureList}}
##' object
##' @param subset Which elements to select from indicated by index, key,
##' location, or qualifier value.
##' @param select Which information to be retrieved from the subsetted
##' elements.
##' 
##' @export
##' @docType methods
setGeneric("select",
           function(x, subset = "", select = "")
             standardGeneric("select")
           )


## initGB-generic =========================================================== 
##
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
           function(db_dir, ...)
             standardGeneric("initGB")
           )
