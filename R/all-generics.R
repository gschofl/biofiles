#' @include utils.R
#' @useDynLib biofiles
#' @importFrom GenomicRanges seqlengths seqinfo
NULL


# The biofiles API -------------------------------------------------------

##    Basic getters/setters in
##    gbLocation, gbFeature, gbFeatureList, gbRecord
##      start, end, strand, width
##      start<-, end<-, strand<-
##      ranges
##
##    Getters/setters in gbLocation-class
##      fuzzy, accession
##
##    Getters/setters in gbFeature-class, gbFeatureList-class
##      index, key, location, ranges, sequence, seqinfo
##      qualif, dbxref, locusTag, product, proteinID, note, translation
##      key<-, qualif<-
##
##    Getters/setters in gbRecord-class
##      locus, reference accession, definition
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

## "start" is defined as an S3 generics in the stats package.
#' Get or set the start of genomic features
#' 
#' @usage start(x, join = FALSE, ...)
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param join Join compound genomic locations into a single range.
#' @param ... Further arguments passed to methods.
#' @return An integer vector or a list of integer vectors.
#' @seealso
#'   \code{\link{end}}, \code{\link{strand}}, \code{\link{width}}, \code{\link{ranges}}
#' @rdname start-methods
#' @export
#' @docType methods
#' @importFrom stats start
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' start(x)
#' 
#' cds <- x["CDS"]
#' start(cds)
#' 
setGeneric("start")

#' @rdname start-methods
#' @export
#' @docType methods
#' @importFrom IRanges "start<-"
setGeneric("start<-")

#' Get or set the end of genomic features
#' 
#' @usage end(x, join = FALSE, ...)
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param join Join compound genomic locations into a single range.
#' @param ... Further arguments passed to methods.
#' @return An integer vector or a list of integer vectors.
#' @seealso
#'   \code{\link{start}}, \code{\link{strand}}, \code{\link{width}}, \code{\link{ranges}}
#' @rdname end-methods
#' @export
#' @docType methods
#' @importFrom stats end
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' end(x)
#' 
#' cds <- x["CDS"]
#' end(cds)
#' 
setGeneric("end")

#' @rdname end-methods
#' @export
#' @docType methods
#' @importFrom IRanges "end<-"
setGeneric("end<-")

### The "strand" generic is defined in the BiocGenerics package.
#' Get or set the strand of genomic features
#'
#' @usage strand(x, join = FALSE, ...)
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param join Join compound genomic locations into a single range.
#' @param ... Further arguments passed to methods.
#' @return An integer vector (or a list thereof) of 1 (plus strand), -1 (minus strand), or
#' \code{NA}
#' @seealso
#'   \code{\link{start}}, \code{\link{end}}, \code{\link{width}}, \code{\link{ranges}}
#' @rdname strand-methods
#' @export
#' @docType methods
#' @importFrom BiocGenerics strand
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' strand(x)
#' 
setGeneric("strand")

#' @rdname strand-methods
#' @export
#' @docType methods
#' @importFrom BiocGenerics "strand<-"
setGeneric("strand<-")

### The "width" generic is defined in the IRanges package. We need
### to define joint_width internally instead of width(x, join=TRUE)
#' Get the width of genomic features
#'
#' @usage width(x)
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @return An integer vector or a list of integer vectors.
#' @seealso
#'   \code{\link{start}}, \code{\link{end}}, \code{\link{strand}}, \code{\link{ranges}}
#' @rdname width-methods
#' @export
#' @docType methods
#' @importFrom IRanges width
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' width(x)
#'
setGeneric("width")

#' @rdname width-methods
#' @export
#' @docType methods
setGeneric("joint_width", signature="x", function(x) standardGeneric("joint_width"))

### The "ranges" generic is defined in the IRanges package.
#' Extract features as \code{"\linkS4class{GRanges}"} objects.
#' 
#' @usage ranges(x, join = FALSE, key = TRUE, include = "none", exclude = "")
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param join Join compound genomic locations into a single range.
#' @param key Include feature keys with ranges.
#' @param include Include qualifiers as metadata columns. Can be "none",
#' "all", or a character vector of qualifier tags.
#' @param exclude Exclude specific qualifiers.
#' @param ... Further arguments passed to methods.
#' @return A \code{\linkS4class{GRanges}} or \code{\linkS4class{GRangesList}} object.
#' @seealso
#'   \code{\link{start}}, \code{\link{end}}, \code{\link{width}}, \code{\link{strand}},
#'   \code{\link{location}}, \code{\link{key}}, \code{\link{qualif}}
#' @rdname ranges-methods
#' @export
#' @importFrom IRanges ranges
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' 
#' ## default "GRanges" object
#' ranges(x)
#' 
#' ## subset CDSs and include "product", "note", "protein_id"
#' ranges(x["CDS"], include=c("product", "note", "protein_id"))
#'
#' ## subset CDSs and exclude "translation"
#' ranges(x["CDS"], include="all", exclude="translation")
#'
setGeneric("ranges")


#' @rdname ranges-methods
#' @export
#' @importFrom IRanges "ranges<-"
setGeneric("ranges<-")

#' Has a feature unclear start/end positions?
#' 
#' @usage fuzzy(x)
#' @rdname fuzzy-methods
#' @export
#' @docType methods
setGeneric("fuzzy", signature="x", function (x, ...) {
  standardGeneric("fuzzy")
})

#' Access the various fields of a GenBank record.
#' 
#' @usage getLocus(x)
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @rdname accessor-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' 
#' getLocus(x)
#' getLength(x)
#' getGeneID(x)
#' getReference(x)
#' 
setGeneric('getLocus', function (x, ...) standardGeneric('getLocus'))

#' @usage getLength(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric("getLength", function (x, ...) standardGeneric('getLength'))

#' @usage getMoltype(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric('getMoltype', function (x, ...) standardGeneric('getMoltype'))

#' @usage getTopology(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric('getTopology', function (x, ...) standardGeneric('getTopology'))

#' @usage getDivision(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric('getDivision', function (x, ...) standardGeneric('getDivision'))

#' @usage getDate(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric('getDate', function (x) standardGeneric('getDate'))

#' @usage getDefinition(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric("getDefinition", function (x, ...) standardGeneric("getDefinition"))

#' @usage getAccession(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric("getAccession", function (x, ...) standardGeneric("getAccession"))

#' @usage getVersion(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric("getVersion", function (x, ...) standardGeneric("getVersion"))

#' @usage getGeneID(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric("getGeneID", function (x, ...) standardGeneric("getGeneID"))

#' @usage getDBLink(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric('getDBLink', function (x) standardGeneric('getDBLink'))

#' @usage getDBSource(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric('getDBSource', function (x) standardGeneric('getDBSource'))

#' @usage getSource(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric('getSource', function (x) standardGeneric('getSource'))

#' @usage getOrganism(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric('getOrganism', function (x) standardGeneric('getOrganism'))

#' @usage getTaxonomy(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric('getTaxonomy', function (x) standardGeneric('getTaxonomy'))

#' @usage getKeywords(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric('getKeywords', function (x) standardGeneric('getKeywords'))

#' @usage getReference(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric('getReference', function (x) standardGeneric('getReference'))

#' @usage getComment(x)
#' @rdname accessor-methods
#' @export
#' @docType methods
setGeneric('getComment', function (x) standardGeneric('getComment'))


# header, features, sequence ------------------------------------------------


#' Get the feature table from a GenBank record.
#'
#' @param x A \code{\linkS4class{gbRecord}} instance.
#' @param ... Additional arguments passed to methods.
#' @return The \code{\linkS4class{gbFeatureList}} of a Genbank record.
#' @rdname getFeatures-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' getFeatures(x)
#' 
setGeneric("getFeatures", function(x, ...) standardGeneric("getFeatures"))


#' @rdname getFeatures-methods
#' @export
#' @docType methods
#' @examples
#' ft(x)
#' 
setGeneric("ft", function(x, ...) standardGeneric("ft"))


#' Get the sequence from a GenBank record.
#' 
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} instance.
#' @param ... Additional arguments passed to methods.
#' @return An \code{\linkS4class{XStringSet}} object, containing either the
#' complete sequence(s) of the record(s), or of the selected feature(s)
#' @rdname getSequence-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' 
#' ## extract the full-length sequence of the record.
#' getSequence(x)
#' 
#' ## extract coding sequences only
#' getSequence(x["CDS"])
#' 
setGeneric("getSequence", function(x, ...) standardGeneric("getSequence"))


#' Extract the header from a \code{"\linkS4class{gbRecord}"} object.
#' 
#' @param x A \code{"\linkS4class{gbRecord}"}, \code{"\linkS4class{gbFeature}"},
#'  or \code{"\linkS4class{gbFeatureList}"} instance.
#' @param ... Additional arguments passed to methods.
#' @return A \code{"\linkS4class{gbHeader}"} instance
#' @rdname getHeader-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' getHeader(x)
#' 
setGeneric("getHeader", function(x, ...) standardGeneric("getHeader"))

#' @rdname getHeader-methods
#' @export
#' @docType methods
#' @examples
#' header(x)
#' 
setGeneric("header", function(x, ...) standardGeneric("header"))


#' Summarise a GenBank record. 
#'
#' @usage summary(object, n = 8, ...)
#' @param object An object of class\code{\linkS4class{gbFeature}},
#' \code{\linkS4class{gbFeatureList}}, \code{\linkS4class{gbRecord}}, or
#' \code{\linkS4class{gbRecordList}}.
#' @param n For \code{list}-like objects, how many elements should
#' be summarized in head and tail.
#' @param ... Arguments to be passed to or from other methods.
#' @rdname summary-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' summary(x)
#' 
setGeneric("summary")


#' Access the indices of GenBank features
#'
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param ... Additional arguments passed to methods.
#' @return A numeric vector of feature indeces.
#' @rdname index-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' index(x)
#' 
setGeneric("index", signature="x", function (x, ...) {
  standardGeneric("index")
})


#' Access genomic locations of GenBank features
#'
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}}, or
#' \code{\linkS4class{gbRecord}} object.
#' @param ... Additional arguments passed to methods.
#' @return A list of \code{\linkS4class{gbLocation}} objects
#' @rdname location-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' location(x)
#'
setGeneric("location", signature="x", function (x, ...) {
  standardGeneric("location")
})


#' Get/set keys of GenBank features
#'
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param ... Additional arguments passed to methods.
#' @rdname key-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' key(x)
#'
setGeneric("key", signature="x", function(x, ...) {
  standardGeneric("key")
})


#' @rdname key-methods
#' @export
#' @docType methods
setGeneric("key<-", signature="x", function(x, value, ...) {
  standardGeneric("key<-")
})


#' Get/set qualifiers of GenBank features
#' 
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}}, or
#' \code{\linkS4class{gbRecord}} object.
#' @param which (Optional) A character vector giving the name(s) of the
#' qualifiers to retrieve.
#' @param fixed If \code{TRUE}, \code{which} is matched against qualifiers as is,
#' if \code{FALSE} it is treated as a regular expression.
#' @param use.names If \code{TRUE}, return a \code{data.frame} using \code{which}
#' as column names, if \code{FALSE} return, if possible, a character vector or
#' a list.
#' @param ... Additional arguments passed to methods.
#' @rdname qualif-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' qualif(x[[1]], 'db_xref')
#' 
#' ## use shortcuts to common qualifiers
#' proteinID(x["CDS"])
#'
setGeneric("qualif", signature=c("x", "which"), function(x, which, ...) {
  standardGeneric("qualif")
})


#' @rdname qualif-methods
#' @export
#' @docType methods
setGeneric("qualif<-", signature=c("x", "which"), function(x, which, value, ...) {
  standardGeneric("qualif<-")
})


#' Access the \code{db_xref}s of GenBank features
#' 
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}}, or
#' \code{\linkS4class{gbRecord}} object.
#' @param db (Optional) A character vector giving the database names of the
#' desired \code{db_xref}s.
#' @param ... Additional arguments passed to methods.
#' @return A named character vector (or list of named character vectors)
#' of db_xrefs.
#' @rdname dbxref-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' 
#' ## all db_xrefs associated with CDSs
#' dbxref(x["CDS"])
#' 
#' ## retrieve the TaxId from the "source" field.
#' dbxref(x[[1]], "taxon")
#'
setGeneric("dbxref", signature="x", function(x, db = NULL, ...) {
  standardGeneric("dbxref")
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
#' @param sequence if \code{FALSE} exclude the sequence.
#' @param append if \code{TRUE} the data is appended to the connection.
#' @export
#' @docType methods
setGeneric("write.GenBank", signature="x",
           function(x, file, header = TRUE, sequence=TRUE, append = FALSE, ...) {
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
#' @docType methods
setGeneric("write.FeatureTable", signature="x",
           function(x, file, tablename = "", dbname = "",
                    sequence = FALSE, append = FALSE, ...) {
             standardGeneric("write.FeatureTable")
           })


# list-generics ----------------------------------------------------------


#' Access the names of Genbank qualifiers
#'
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}},
#' or \code{\linkS4class{gbRecord}} object.
#' @param ... Additional arguments to be passed to or from methods.
#' @return A character vector (or list of character vectors) of 
#'    qualifiers.
#' @seealso
#'    \code{\link{listUniqueQualifs}}, \code{\link{hasQualif}}
#' @rdname listQualif-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' listQualif(x["source"])
#'
setGeneric("listQualif", signature="x", function(x, ...) {
  standardGeneric("listQualif")
})


# test-generics ----------------------------------------------------------

#' Has a feature a specific key?
#'  
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param ... Additional arguments to be passed to or from methods.
#' @return A logical vector or a list of logical vectors.
#' @seealso
#'    \code{\link{key}}
#' @rdname hasKey-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' hasKey(x, 'CDS')
#'
setGeneric("hasKey", signature=c("x","key"), function(x, key, ...) {
  standardGeneric("hasKey")
})

#' Has a feature a specific qualifier?
#'
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureList}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param ... Additional arguments to be passed to or from methods.
#' @return A logical vector or a list of logical vectors.
#' @seealso
#'  \code{\link{listQualif}}, to extract a list of available qualifiers for
#'  each feature; \code{\link{listUniqueQualif}}, for a vector of all
#'  unique qualifiers present in an object. 
#' @rdname hasQualif-methods
#' @export
#' @docType methods
#' @examples
#' gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
#' x <- gbRecord(gbk_file)
#' hasQualif(x, 'CDS')
#'
setGeneric("hasQualif", signature=c("x","qualifier"), function(x, qualifier, ...) {
  standardGeneric("hasQualif")
})


# shift and revcomp -------------------------------------------------------------


#' Shift location of features in a GenBank record
#'
#' @usage shift(x, shift=0L, split=FALSE, order=FALSE)
#' @param x A \code{\linkS4class{gbFeatureList}} or
#' \code{\linkS4class{gbRecord}} instance (gbFeatureLists must 
#' be complete and include a 'source' field).
#' @param shift Number of basepairs (or aa residues) to shift.
#' @param split Split features that after the shift extends across the end of
#' the sequence.
#' @param order Reorder features after the shift.
#'
#' @return A \code{\linkS4class{gbFeatureList}} object.
#' @rdname shift-methods
#' @export
#' @docType methods
setGeneric("shift", signature="x", function(x, shift=0L, use.names=TRUE, ...) {
  standardGeneric("shift")
})

#' Reverse-complement features in a GenBank record
#' 
#' @usage revcomp(x, order=TRUE)
#' @param x A \code{\linkS4class{gbFeatureList}} or
#' \code{\linkS4class{gbRecord}} object (gbFeatureLists must 
#' be complete and include a 'source' field).
#' @param order Reorder features after reverse-complementing them.
#' @rdname revcomp-methods
#' @export
#' @docType methods
setGeneric("revcomp", signature="x", function(x, ...) standardGeneric("revcomp"))


# view -------------------------------------------------------------------


#' View all features in a \code{gbFeatureList}
#' 
#' @param x A \code{\linkS4class{gbFeatureList}} instance.
#' @param n How many features to show (Default: all).
#' @param ... Additional arguments passed to methods.
#' @rdname view-methods
#' @keywords internal
#' @export
#' @docType methods
setGeneric("view", signature="x", function(x, n, ...) {
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
#' @rdname select-methods
#' @export
#' @docType methods
setGeneric("select", signature="x", function(x, ..., keys = NULL, cols = NULL) {
  standardGeneric("select")
})


# internal ---------------------------------------------------------------


#' @keywords internal
setGeneric('.seqinfo', function(x) standardGeneric('.seqinfo'))

#' @keywords internal
setGeneric('.header', function(x) standardGeneric('.header'))

#' @keywords internal
setGeneric('.sequence', function(x) standardGeneric('.sequence'))

#' @keywords internal
setGeneric('.locus', function(x) standardGeneric('.locus'))

#' @keywords internal
setGeneric('.features', function(x) standardGeneric('.features'))

#' @keywords internal
setGeneric('.contig', function(x) standardGeneric('.contig'))

#' @keywords internal
setGeneric('.dbSource', function(x) standardGeneric('.dbSource'))

#' @keywords internal
setGeneric('.defline', function(x) standardGeneric('.defline'))

