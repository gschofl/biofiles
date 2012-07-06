
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


##' Retrieve feature locations from a GenBank record
##'
##' @usage getLocation(x, attributes = TRUE, join = FALSE)
##'
##' @param x A \code{\link{gbFeature}} or \code{\link{gbFeatureList}} object
##' @param attributes Set the \code{accession}, \code{definition},
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

##' Retrieve feature indices from a GenBank record
##'
##' @usage getIndex(x, attributes = FALSE)
##'
##' @param x A \code{gbFeature} or \code{gbFeatureList} object.
##' @param attributes Set the \code{accession}, \code{definition},
##' \code{database} attributes.
##'
##' @return A numeric vector of feature indeces
##'
##' @docType methods
##' @export
setGeneric("getIndex",
           function (x, attributes = FALSE, ...) 
             standardGeneric("getIndex")
           )

##' Retrieve feature keys from a GenBank record
##'
##' @usage getKey(x, attributes = FALSE)
##' 
##' @param x A \code{gbFeature} or \code{gbFeatureList} object.
##' @param attributes Set the \code{accession}, \code{definition},
##' \code{database} attributes.
##'
##' @docType methods
##' @export
setGeneric("getKey",
           function(x, attributes = FALSE, ...) 
             standardGeneric("getKey")
           )

##' Retrieve feature keys from a GenBank record
##'
##' @usage getKey(x, attributes = FALSE)
##' 
##' @param x A \code{gbFeature} or \code{gbFeatureList} object.
##' @param which Character vector giving the name(s) of the qualifiers
##' to retrieve.
##' @param attributes Set the \code{accession}, \code{definition},
##' \code{database} attributes.
##'
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
##' @return A \code{\link{gbLocation-class}}, \code{\link{gbFeature-class}},
##' or\code{\link{gbFeatureList-class}} object
##'
##' @docType methods
##' @export
setGeneric("shift", function (x, shift = 0L, ...) standardGeneric("shift"))


## select-generic =========================================================== 
##
##' Select elements from a GenBank Record
##' 
##' This function extracts features or information contained in features
##' from \code{gbRecord}s or \code{gbFeatureList}s.
##' 
##' \sQuote{subset} specifications are given in a character string separated
##' by semicolons.
##' 
##' \describe{
##'    \item{by index}{
##'    For example: \code{"idx=1,3,4"}, \code{"idx=1:4"},
##'    \code{"index=1,3,5:8"}, \code{"idx=1;idx=3;idx=4"}
##'    }
##'    \item{by location}{
##'    For example: \code{"loc=10000:20000"}, \code{"loc=:20000"}, 
##'    \code{"location=:20000,40000:60000"}
##'    }
##'    \item{by key}{
##'    For example: \code{"key=CDS"}, \code{"key=CDS,gene"}
##'    }
##'    \item{by qualifier values}{
##'    For example: \code{"product=replication"}, \code{"pseudo"}
##'    }
##'    \item{by any combination of the above}{
##'    For example: \code{"location=:10000;key=CDS,gene"}
##'    }
##' }
##' 
##' @usage select(x, subset = c(""), select =c(""))
##' 
##' @param x A \sQuote{\code{gbRecord}} or \sQuote{\code{gbFeatureList}}
##' object
##' @param subset A character string specifying the elements to select
##' indicated by index, key, location, or qualifier value. See Details.
##' @param select A character string specifying the information to be
##' returned from the subsetted elements. Supported values are \dQuote{idx}
##' or \dQuote{index}, \dQuote{location} or \dQuote{range}, \dQuote{start},
##' \dQuote{end}, \dQuote{strand}, \dQuote{key}, or any qualifier tag (e.g.,
##' \dQuote{locus_tag}, \dQuote{product}). The default will return a
##' \code{\link{gbFeature-class}} or \code{\link{gbFeatureList-class}} object
##' 
##' @return Depending on the value of \code{select}.
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
