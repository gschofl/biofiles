#' @include utils.R
NULL


# The biofiles API -------------------------------------------------------

##    Basic getters/setters in
##    gbLocation-class, gbFeature-class, gbFeatureList-class
##      start, end, strand, width
##      start<-, end<-, strand<-
##
##    Getters/setters in gbLocation-class
##      ranges, fuzzy, accession
##
##    Getters/setters in gbFeature-class, gbFeatureList-class
##      index, key, location, ranges, sequence, seqinfo
##      qualif, dbxref, locusTag, product, proteinID, note, translation
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
##    Genetic Write methods
##      write.GenBank, write.FeatureTable
##
##    The "start" and "end" generics are defined in the stats package.
##    
##    The "shift", "start<-", "end<-", and "ranges" generics are defined in the
##    IRanges package.
##
##    We need to override the "width" and "shift" generics from IRanges because
##    they don't provide a ... argument.
##    We also redefine the "sequence" generic from Base.

# getter/setter generics -------------------------------------------------

### "start" and "end" are defined as S3 generics in the stats package.
### Here we explicitely set them S4.
#' Get or set the start of genomic features
#' 
#' @usage start(x, join=FALSE, ...)
#' @param x A \code{gbFeature} or \code{gbFeatureList} object.
#' @param join Join compound genomic locations into a single range.
#' @param ... Further arguments passed to methods.
#' @return An integer vector
#' @rdname start
#' @export
#' @genericMethods
#' @importFrom IRanges start
setGeneric("start")


#' @rdname start
#' @export
#' @genericMethods
#' @importFrom IRanges "start<-"
setGeneric("start<-")


#' Get or set the end of genomic features
#' 
#' @usage start(x, join=FALSE, ...)
#' @param x A \code{gbFeature} or \code{gbFeatureList} object.
#' @param join Join compound genomic locations into a single range.
#' @param ... Further arguments passed to methods.
#' @return An integer vector
#' @rdname end
#' @export
#' @genericMethods
#' @importFrom IRanges end
setGeneric("end")


#' @rdname end
#' @export
#' @genericMethods
#' @importFrom IRanges "end<-"
setGeneric("end<-")


### The "strand" generic is defined in the BiocGenerics package.
#' Get or set the strand of genomic features
#'
#' @usage strand(x, join=FALSE, ...)
#' @param x A \code{gbFeature} or \code{gbFeatureList}object.
#' @param join Join compound genomic locations into a single range.
#' @param ... Further arguments passed to methods.
#' @return An integer vector of 1 (plus strand), -1 (minus strand), or
#' \code{NA}
#' @rdname strand
#' @export
#' @genericMethods
#' @importFrom BiocGenerics strand
setGeneric("strand")


#' @rdname strand
#' @export
#' @genericMethods
#' @importFrom BiocGenerics "strand<-"
setGeneric("strand<-")


### The "width" generic is defined in the IRanges package. We need
### to override it because they don't provide a dotdotdot interface.
#' Get the width of genomic features
#'
#' @usage width(x, ...)
#' @param x A \code{gbFeature} or \code{gbFeatureList}object.
#' @param ... Further arguments passed to methods.
#' @return An integer vector.
#' @rdname width
#' @export
#' @genericMethods
#' @importFrom IRanges width
setGeneric("width", signature="x", function (x, ...) {
  standardGeneric("width")
})


### The "ranges" generic is defined in the IRanges package.
#' Get or set the range of genomic features
#' 
#' @usage ranges(x, join=FALSE, include="none", exclude="")
#' @param x A \code{gbFeature} or \code{gbFeatureList} object.
#' @param join Join compound genomic locations into a single range.
#' @param key Include feature keys with ranges.
#' @param include Include qualifiers as metadata columns. Can be "none",
#' "all", or a character vector of qualifier tags.
#' @param exclude Exclude specific qualifiers.
#' @param ... Further arguments passed to methods.
#' @return A \code{\linkS4class{GRanges}} object.
#' @rdname ranges
#' @export
#' @importFrom IRanges ranges
setGeneric("ranges")


#' @rdname ranges
#' @export
#' @importFrom IRanges "ranges<-"
setGeneric("ranges<-")


#' @rdname fuzzy
#' @export
#' @genericMethods
setGeneric("fuzzy", signature="x", function (x, ...) {
  standardGeneric("fuzzy")
})


### The "seqinfo" generic is defined in the BiocGenerics package.
#' Get sequence information about genomic features
#' 
#' @usage seqinfo(x)
#' @rdname seqinfo
#' @export
#' @genericMethods
#' @importFrom GenomicRanges seqinfo
setGeneric("seqinfo")


#' @usage accession(x)
#' @rdname seqinfo
#' @export
#' @genericMethods
setGeneric("accession", signature="x", function (x, ...) {
  standardGeneric("accession")
})


#' @usage definition(x)
#' @rdname seqinfo
#' @export
#' @genericMethods
setGeneric("definition", signature="x", function (x, ...) {
  standardGeneric("definition")
})


#' @usage seqlengths(x)
#' @rdname seqinfo
#' @export
#' @genericMethods
#' @importFrom GenomicRanges seqlengths
setGeneric("seqlengths")


### The "annotation" generic is defined in the BiocGenerics package.
#' @rdname annotation
#' @export
#' @genericMethods
#' @importFrom BiocGenerics annotation
setGeneric("annotation")


#' @rdname summary
#' @export
#' @genericMethods
setGeneric("summary")


#' Get indices of GenBank features
#'
#' @param x A \code{\linkS4class{gbFeature}} or
#' \code{\linkS4class{gbFeatureList}} instance.
#' @param ... Additional arguments passed to methods.
#' @return A numeric vector of feature indeces.
#' @rdname index
#' @export
#' @genericMethods
setGeneric("index", signature="x",
           function (x, ...) {
             standardGeneric("index")
           })


#' Get genomic locations of GenBank features
#'
#' @param x A \code{\linkS4class{gbFeature}} or
#' \code{\linkS4class{gbFeatureList}} instance.
#' @param ... Additional arguments passed to methods.
#' @return A list of \code{\linkS4class{gbLocation}} objects
#' @rdname location
#' @export
#' @genericMethods
setGeneric("location", signature="x",
           function (x, ...) {
             standardGeneric("location")
           })


#' Get/set keys of GenBank features
#'
#' @param x A \code{\linkS4class{gbFeature}} or
#' \code{\linkS4class{gbFeatureList}} instance.
#' @param ... Additional arguments passed to methods.
#' @rdname key
#' @export
#' @genericMethods
setGeneric("key", signature="x",
           function(x, ...) {
             standardGeneric("key")
             })


#' @rdname key
#' @export
#' @genericMethods
setGeneric("key<-", signature="x",
           function(x, value, ...) {
             standardGeneric("key<-")
           })


#' Get/set qualifiers of GenBank features
#' 
#' @param x A \code{\linkS4class{gbFeature}} or
#' \code{\linkS4class{gbFeatureList}} instance.
#' @param which (Optional) A character vector giving the name(s) of the
#' qualifiers to retrieve.
#' @param ... Additional arguments passed to methods.
#' @rdname qualif
#' @export
#' @genericMethods
setGeneric("qualif", signature=c("x", "which"),
           function(x, which, ...) {
             standardGeneric("qualif")
           })


#' @rdname qualif
#' @export
#' @genericMethods
setGeneric("qualif<-", signature=c("x", "which"),
           function(x, which, value, ...) {
             standardGeneric("qualif<-")
           })


#' Get the \code{db_xref}s of GenBank features
#' 
#' @param x A \code{\linkS4class{gbFeature}} or
#' \code{\linkS4class{gbFeatureList}} instance.
#' @param db (Optional) A character vector giving the database names of the
#' desired \code{db_xref}s.
#' @param ... Additional arguments passed to methods.
#' @rdname dbxref
#' @export
#' @genericMethods
setGeneric("dbxref", signature="x",
           function(x, db = NULL, ...) {
             standardGeneric("dbxref")
           })


#' Get sequences of GenBank features
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


# write generics ---------------------------------------------------------


#' Write GenBank records or features to file in GenBank format
#'
#' @details
#' For a description of the GenBank format see
#' \url{http://www.ncbi.nlm.nih.gov/collab/FT/}
#' 
#' @param x A \code{\linkS4class{gbRecord}} instance.
#' @param file A connection or a character string naming the file to write to.
#' @param header if \code{FALSE} exclude the Genbank header.
#' @param append if \code{TRUE} the data is appended to the connection.
#' @export
#' @genericMethods
setGeneric("write.GenBank", signature="x",
           function(x, file, header = TRUE, append = FALSE, ...) {
             standardGeneric("write.GenBank")
           })


#' Write GenBank records or features to file in Feature Table format
#'
#' Feature Tables are simple five-column tab-delimited tables specifying the
#' location and type of each feature. They can be used as input for tbl2asn
#' or Sequin to generate annotation.
#' 
#' @param x A \code{\linkS4class{gbRecord}} instance.
#' @param file A connection or a character string naming the file to write to.
#' @param tablename Optional table name to appear in the first line
#' of the feature table.
#' @param dbname Data base name associated with the CDS qualifier protein_id.
#' @param sequence if \code{TRUE}, additionally autput fasta file
#' @param append if \code{TRUE} the data is appended to the connection.
#' @export
#' @genericMethods
setGeneric("write.FeatureTable", signature="x",
           function(x, file, tablename = "", dbname = "",
                    sequence = FALSE, append = FALSE, ...) {
             standardGeneric("write.FeatureTable")
           })


# list-generics ----------------------------------------------------------


#' @rdname listQualif
#' @export
#' @genericMethods
setGeneric("listQualif", signature="x", function(x, ...) {
  standardGeneric("listQualif")
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
#' @usage shift(x, shift=0L, split=FALSE, order=FALSE)
#'
#' @param x A \code{\linkS4class{gbFeatureList}} or
#' \code{\linkS4class{gbRecord}} instance (gbFeatureLists must 
#' be complete and include a 'source' field).
#' @param shift Number of basepairs (or aa residues) to shift.
#' @param split Split features that after the shift extends across the end of
#' the sequence.
#' @param order Reorder features after the shift.
#'
#' @return A \code{\linkS4class{gbFeatureList}} object.
#' @rdname shift
#' @export
#' @genericMethods
setGeneric("shift", signature="x",
           function(x, shift=0L, use.names=TRUE, ...) {
             standardGeneric("shift")
           })


# view -------------------------------------------------------------------


#' View all features in a \code{gbFeatureList}
#' 
#' @param x A \code{\linkS4class{gbFeatureList}} instance.
#' @param n How many features to show (Default: all).
#' @param ... Additional arguments passed to methods.
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

