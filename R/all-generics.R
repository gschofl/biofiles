#' @include utils.r
#' @include validate.r
NULL


# The biofiles API -------------------------------------------------------

##    Basic getters/setters in
##    gbLocation-class, gbRange-class, gbFeature-class, gbFeatureList-class
##      start, end, strand, width
##      start<-, end<-, strand<-
##
##    Getters/setters in gbLocation-class
##      range, partial, accession
##
##    Getters/setters in gbFeature-class, gbFeatureList-class
##      index, key, range, location, qualif, dbxref, sequence,
##      accession, definition
##      key<-, qualif<-
##
##    Getters/setters in gbRecord-class
##      features, sequence, accession, definition
##
##    Testing methods in gbFeature-class, gbFeatureList-class
##      hasKey, hasQualif
##
##    Show methods all classes
##      show
##
##    Subsetting Methods
##      [, [[, $
##
##    The "start" and "end" generics are defined in the stats package.
##    The "shift", "start<-", and "end<-" generics are defined in the
##    IRanges package.
##    The "range" generic is defined in the base package
##
##    We need to override the "width" and "shift" generics from IRanges because
##    they don't provide a ... argument.
##    We also redefine the "sequence" generic from Base.

# getter/setter generics -------------------------------------------------

### The "start" and "end" generics are defined in the stats package.
#' @rdname start
#' @export
#' @genericMethods
setGeneric("start")


#' @rdname start
#' @export
#' @genericMethods
setGeneric("start<-")


#' @rdname end
#' @export
#' @genericMethods
setGeneric("end")


#' @rdname end
#' @export
#' @genericMethods
setGeneric("end<-")


### The "strand" generic is defined in the BiocGenerics package.
#' @rdname strand
#' @export
#' @genericMethods
setGeneric("strand")


#' @rdname strand
#' @export
#' @genericMethods
setGeneric("strand<-")

### The "width" generic is defined in the IRanges package. We need
### to override it because they don't provide a dotdotdot interface.
#' @rdname width
#' @export
#' @genericMethods
setGeneric("width", signature="x", function (x, ...) {
  standardGeneric("width")
})


#' @rdname partial
#' @export
#' @genericMethods
setGeneric("partial", signature="x", function (x, ...) {
  standardGeneric("partial")
})


#' @rdname accession
#' @export
#' @genericMethods
setGeneric("accession", signature="x", function (x, ...) {
  standardGeneric("accession")
})


#' @rdname definition
#' @export
#' @genericMethods
setGeneric("definition", signature="x", function (x, ...) {
  standardGeneric("definition")
})


#' @rdname location
#' @export
#' @genericMethods
setGeneric("location", signature="x",
           function (x, attributes = FALSE, join = FALSE, ...) {
             standardGeneric("location")
           })


#' @rdname annotation
#' @export
#' @genericMethods
setGeneric("annotation")


#' Return feature indices from a GenBank record
#'
#' @param x A \code{\linkS4class{gbFeature}} or
#' \code{\linkS4class{gbFeatureList}} instance.
#' @param attributes Set the \code{accession}, \code{definition},
#' \code{database} attributes.
#' @param ... Additional arguments passed to methods.
#' @return A numeric vector of feature indeces.
#' @rdname index
#' @export
#' @genericMethods
setGeneric("index", signature="x",
           function (x, attributes = FALSE, ...) {
             standardGeneric("index")
           })


#' Get/set feature keys from a GenBank Record
#'
#' @param x A \code{\linkS4class{gbFeature}} or
#' \code{\linkS4class{gbFeatureList}} instance.
#' @param attributes Set the \code{accession}, \code{definition},
#' \code{database} attributes.
#' @param ... Additional arguments passed to methods.
#' @rdname key
#' @export
#' @genericMethods
setGeneric("key", signature="x",
           function(x, attributes = FALSE, ...) {
             standardGeneric("key")
             })


#' @rdname key
#' @export
#' @genericMethods
setGeneric("key<-", signature="x",
           function(x, value, ...) {
             standardGeneric("key<-")
           })


#' Get/set feature qualifiers from a GenBank record
#' 
#' @param x A \code{\linkS4class{gbFeature}} or
#' \code{\linkS4class{gbFeatureList}} instance.
#' @param which (Optional) A character vector giving the name(s) of the
#' qualifiers to retrieve.
#' @param attributes Set the \code{accession}, \code{definition},
#' \code{database} attributes.
#' @param ... Additional arguments passed to methods.
#' @rdname qualif
#' @export
#' @genericMethods
setGeneric("qualif", signature=c("x", "which"),
           function(x, which, attributes = FALSE, ...) {
             standardGeneric("qualif")
           })


#' @rdname qualif
#' @export
#' @genericMethods
setGeneric("qualif<-", signature=c("x", "which"),
           function(x, which, value, ...) {
             standardGeneric("qualif<-")
           })


#' Get/set db_xref from a GenBank record
#' 
#' @param x A \code{\linkS4class{gbFeature}} or
#' \code{\linkS4class{gbFeatureList}} instance.
#' @param db (Optional) A character vector giving the database names of the
#' desired db_xrefs.
#' @param ... Additional arguments passed to methods.
#' @rdname dbxref
#' @export
#' @genericMethods
setGeneric("dbxref", signature="x",
           function(x, db = NULL, ...) {
             standardGeneric("dbxref")
           })


#' Get sequences of a GenBank records or features.
#'
#' @param x A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeature}},
#'  or \code{\linkS4class{gbFeatureList}} instance.
#' @param ... Additional arguments passed to methods.
#' @return An \code{\linkS4class{XStringSet}} object.
#' @rdname sequence
#' @export
#' @genericMethods
setGeneric("sequence", signature="x", function(x, ...) {
  standardGeneric("sequence")
})


#' Retrieve features of a GenBank record.
#'
#' @param x A \code{\linkS4class{gbRecord}} instance.
#' @param ... Additional arguments passed to methods.
#' @return The \code{\linkS4class{gbFeatureList}} object contained in a
#' gbRecord database.
#' @rdname features
#' @export
#' @genericMethods
setGeneric("features", signature="x", function (x, ...) {
  standardGeneric("features")
})


# test-generics ----------------------------------------------------------


#' @rdname hasKey
#' @export
#' @genericMethods
setGeneric("hasKey", signature=c("x","key"), function(x, key, ...) {
  standardGeneric("hasKey")
})


#' @rdname hasQualif
#' @export
#' @genericMethods
setGeneric("hasQualif", signature=c("x","qualifier"),
           function(x, qualifier, ...) {
             standardGeneric("hasQualif")
           })


# shift ------------------------------------------------------------------



#' Shift location of features in a GenBank record
#'
#' @usage shift(x, shift=0L, split=FALSE, order=FALSE, updateDb=FALSE)
#'
#' @param x A gbLocation, gbFeature, gbFeatureList, or gbRecord instance
#' (gbFeatureLists must include a 'source' field).
#' @param shift Number of basepairs (or aa residues) to shift.
#' @param split (For gbFeatureList and gbRecord objects) Should a feature
#' that spans across the end of the sequence be split.
#' @param order (For gbFeatureList and gbRecord objects) Should the
#' resulting gbFeatureList be reordered.
#' @param updateDb Should filehash database be updated with new feature
#' locations.
#'
#' @return A \code{\link{gbLocation-class}}, \code{\link{gbFeature-class}},
#' or\code{\link{gbFeatureList-class}} object
#' @rdname shift
#' @export
#' @genericMethods
setGeneric("shift", signature="x",
           function(x, shift=0L, use.names=TRUE, ...) {
             standardGeneric("shift")
           })


# revcomp ----------------------------------------------------------------


#' reverse complement features in a GenBank record
#'
#' @param x A \code{\linkS4class{gbRecord}} or
#' \code{\linkS4class{gbFeatureList}} or instance.
#' (\code{gbFeatureList} instances must include a \sQuote{source} field).
#' @param order Should the resulting gbFeatureList be reordered.
#' @param update_db Should the sequence and the new feature locations be
#' updated in the underlying filehash database?
#' @param ... Additional arguments passed to methods.
#' 
#' @return A \code{\linkS4class{gbFeatureList}} instance
#' @rdname revcomp
#' @export
#' @genericMethods
setGeneric("revcomp", signature="x",
           function (x, order=FALSE, updateDb=FALSE, ...){
             standardGeneric("revcomp") 
           })


# view -------------------------------------------------------------------


#' View all features in a gbFeatureList
#' 
#' @param x A \code{\linkS4class{gbFeatureList}} instance.
#' @param n How many features to show (Default: all).
#' @param ... Additional arguments passed to methods.
#' @return NULL
#' @rdname view
#' @export
#' @genericMethods
setGeneric("view", signature="x", function (x, n, ...){
  standardGeneric("view")
})


# select -----------------------------------------------------------------


#' Select features from a GenBank Record
#' 
#' This function extracts features or information contained in features
#' from \code{gbRecord}s or \code{gbFeatureList}s.
#' 
#' @details
#' Queries can be are provided as named values using predefined
#' keywords and \dQuote{qualifier=value} pairs:
#' 
#' Possible keywords are:
#' 
#' \describe{
#'   \item{index/idx}{
#'     For example: \code{idx=c(3,4,5,6)}, \code{idx=100:150},
#'     \code{index=c(1,12:20)}
#'   }
#'   \item{range}{
#'     For example: \code{range="10000..25000"},
#'     \code{range="..10000,20000..25000"},
#'     \code{range="30000.."}
#'   }
#'   \item{key}{
#'     For example: \code{key="CDS"}, \code{key=c("CDS", "gene")},
#'     \code{key="CDS|gene"}
#'   }
#'   \item{arbitrary qualifiers}{
#'     For example: \code{product="ribosomal"}, \code{locus_tag=c("CPSIT_0123",
#'     "CPSIT_0124", "CPSIT_0125")}, \code{pseudo=TRUE}
#'   }
#' }
#' 
#' Alternatively (or additionally) queries can be described in a single
#' character string of \dQuote{tag=value} pairs passed to the \code{keys}
#' argument.
#' 
#' Different \dQuote{tag=value} pairs are separated by semicolons.
#' 
#' \describe{
#'   \item{index}{
#'     For example: \code{"idx=1,3,4"}, \code{"idx=1:4"},
#'     \code{"index=1,3,5:8"}, \code{"idx=1;idx=3;idx=4"}
#'   }
#'   \item{by range}{
#'     For example: \code{"range=10000..20000"}, \code{"range=..20000"}, 
#'     \code{"range=..20000,40000..60000"}
#'   }
#'   \item{by key}{
#'       For example: \code{"key=CDS"}, \code{"key=CDS,gene"}
#'   }
#'   \item{by qualifier values}{
#'     For example: \code{"product=replication"}, \code{"pseudo"}
#'   }
#'   \item{by any combination of the above}{
#'     For example: \code{"range=..10000;key=CDS,gene"}
#'   }
#' }
#' 
#' @param x A \sQuote{\code{gbRecord}} or \sQuote{\code{gbFeatureList}}
#' instance
#' @param ... Named values that specify the features to select. These are
#' merged with the values of \code{keys} to create the actual query. See
#' Details.
#' @param keys A character string composed of \sQuote{\emph{key=value}}-pairs
#' separated by semicolons that specify the elements to select. See Details.
#' @param cols A character string of \sQuote{\emph{keys}} that
#' indicate the data to be retrieved from the selected features.
#' Supported \sQuote{\emph{keys}} are \dQuote{idx} or \dQuote{index},
#' \dQuote{location} or \dQuote{range}, \dQuote{start}, \dQuote{end},
#' \dQuote{strand}, \dQuote{key}, or any qualifier tag (e.g.,
#' \dQuote{locus_tag}, \dQuote{product}). If no \sQuote{\emph{keys}} is given
#' a  \code{\linkS4class{gbFeature}} or \code{\linkS4class{gbFeatureList}}
#' object is returned.
#' 
#' @return Depending on the value of \code{select}.
#' @rdname select
#' @export
#' @genericMethods
setGeneric("select", signature="x",
           function(x, ..., keys = NULL, cols = NULL) {
             standardGeneric("select")
           })

