#' @include utils.R
#' @useDynLib biofiles
NULL

# The biofiles API -------------------------------------------------------

##    Basic getters/setters in gbLocation-class, gbFeature-class, gbFeatureTable-class,
##    gbRecord-class:
##      start, end, span, strand, joint_range, start<-, end<-, strand<-, ranges
##
##    Getters/setters in gbLocation-class:
##      fuzzy, accession
##
##    Getters/setters in gbFeature-class, gbFeatureTable-class:
##      index, key, location, ranges, sequence, seqinfo
##      qualif, dbxref, locusTag, product, proteinID, note, translation
##      key<-, qualif<-
##
##    Getters/setters in gbRecord-class:
##      locus, reference accession, definition
##
##    Testing methods in gbFeature-class, gbFeatureTable-class:
##      hasKey, hasQualif
##
##    Show methods all classes:
##      show
##
##    Subsetting Methods:
##      [, [[, $
##
##    Genetic Write methods:
##      write.GenBank, write.FeatureTable
##    
##    The "shift" and "ranges" generics are defined in the IRanges package.

# getter/setter generics -------------------------------------------------

#' Get or set the start position of genomic features
#' 
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param join Join compound genomic locations into a single range.
#' @param ... Further arguments passed to methods.
#' @return An integer vector or a list of integer vectors.
#' @seealso
#'   \code{\link{end}}, \code{\link{strand}}, \code{\link{span}}, \code{\link{ranges}}
#' @rdname start-methods
#' @importFrom stats start
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' 
#' start(x)
#' cds <- x["CDS"]
#' start(cds)
setGeneric("start")

#' @param check if \code{FALSE}, don't perform validity checks.
#' @param value The start information to set on \code{x}.
#' @importFrom IRanges "start<-"
#' @rdname start-methods
setGeneric("start<-")

#' Get or set the end position of genomic features
#' 
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param join Join compound genomic locations into a single range.
#' @param ... Further arguments passed to methods.
#' @return An integer vector or a list of integer vectors.
#' @seealso
#'   \code{\link{start}}, \code{\link{strand}}, \code{\link{span}}, \code{\link{ranges}}
#' @rdname end-methods
#' @importFrom stats end
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' 
#' end(x)
#' cds <- x["CDS"]
#' end(cds)
setGeneric("end")

#' @param check if \code{FALSE}, don't perform validity checks.
#' @param value The end information to set on \code{x}.
#' @importFrom IRanges "end<-"
#' @rdname end-methods
setGeneric("end<-")

### The "strand" generic is defined in the BiocGenerics package.
#' Get or set the strand of genomic features
#'
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param join Join compound genomic locations into a single range.
#' @param ... Further arguments passed to methods.
#' @return An integer vector (or a list thereof) of 1 (plus strand), -1 (minus strand), or
#' \code{NA}
#' @seealso
#'   \code{\link{start}}, \code{\link{end}}, \code{\link{span}}, \code{\link{ranges}}
#' @rdname strand-methods
#' @importFrom BiocGenerics strand
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' strand(x)
setGeneric("strand")

#' @rdname strand-methods
#' @importFrom BiocGenerics strand<-
setGeneric("strand<-")

#' Get the span of genomic features.
#'
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param join Join compound genomic locations into a single range.
#' @param ... Further arguments passed to methods.
#' @return An integer vector or a list of integer vectors.
#' @seealso
#'   \code{\link{start}}, \code{\link{end}}, \code{\link{strand}}, \code{\link{ranges}}
#' @rdname span-methods
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' span(x)
setGeneric("span", signature = "x", function (x, ...) {
  standardGeneric("span")
})

#' @rdname span-methods
#' @export
setGeneric("joint_range", signature = "x", function(x) standardGeneric("joint_range"))

### The "ranges" generic is defined in the IRanges package.
#' Extract features as \code{"\linkS4class{GRanges}"} objects.
#' 
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param join Join compound genomic locations into a single range.
#' @param key Include feature keys with ranges.
#' @param include Include qualifiers as metadata columns. Can be "none",
#' "all", or a character vector of qualifier tags.
#' @param exclude Exclude specific qualifiers.
#' @param ... Further arguments passed to methods.
#' @return A \code{\linkS4class{GRanges}} or \code{\linkS4class{GRangesList}} object.
#' @seealso
#'   \code{\link{start}}, \code{\link{end}}, \code{\link{span}}, \code{\link{strand}},
#'   \code{\link{location}}, \code{\link{key}}, \code{\link{qualif}}
#' @rdname ranges-methods
#' @importFrom IRanges ranges
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' 
#' ## default "GRanges" object.
#' ranges(x)
#' 
#' ## subset CDSs and include "product", "note", "protein_id" as metadata.
#' ranges(x["CDS"], include = c("product", "note", "protein_id"))
#'
#' ## subset CDSs and exclude "translation"
#' ranges(x["CDS"], include = "all", exclude = "translation")
setGeneric("ranges")

#' Has a feature fuzzy locations?
#' 
#' With a GenBank location like \emph{complement(<123..150)} we don't know the exact
#' start position of the feature. Use \code{fuzzy} to test for fuzzy locations.
#' 
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param ... Further arguments passed to methods.
#' @return A logical matrix.
#' @rdname fuzzy-methods
#' @export
#' @examples
#' l <- as.gbLocation("complement(<123..150)")
#' fuzzy(l)
#' 
#' ## note that start() or end() return exact positions even if they are fuzzy.
#' start(l)
setGeneric("fuzzy", signature = "x", function (x, ...) {
  standardGeneric("fuzzy")
})

#' Access the various fields of a GenBank record.
#' 
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param ... Further arguments passed to methods.
#' @rdname accessor-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' 
#' getLocus(x)
#' getLength(x)
#' getGeneID(x)
#' getReference(x)
#' getDate(x)
setGeneric('getLocus', function (x, ...) standardGeneric('getLocus'))

#' @rdname accessor-methods
#' @export
setGeneric("getLength", function (x, ...) standardGeneric('getLength'))

#' @rdname accessor-methods
#' @export
setGeneric('getMoltype', function (x, ...) standardGeneric('getMoltype'))

#' @rdname accessor-methods
#' @export
setGeneric('getTopology', function (x, ...) standardGeneric('getTopology'))

#' @rdname accessor-methods
#' @export
setGeneric('getDivision', function (x, ...) standardGeneric('getDivision'))

#' @rdname accessor-methods
#' @export
setGeneric('getDate', function (x) standardGeneric('getDate'))

#' @rdname accessor-methods
#' @export
setGeneric("getDefinition", function (x, ...) standardGeneric("getDefinition"))

#' @rdname accessor-methods
#' @export
setGeneric("getAccession", function (x, ...) standardGeneric("getAccession"))

#' @rdname accessor-methods
#' @export
setGeneric("getVersion", function (x, ...) standardGeneric("getVersion"))

#' @rdname accessor-methods
#' @export
setGeneric("getGeneID", function (x, ...) standardGeneric("getGeneID"))

#' @rdname accessor-methods
#' @export
setGeneric('getDBLink', function (x) standardGeneric('getDBLink'))

#' @rdname accessor-methods
#' @export
setGeneric('getDBSource', function (x) standardGeneric('getDBSource'))

#' @rdname accessor-methods
#' @export
setGeneric('getSource', function (x) standardGeneric('getSource'))

#' @rdname accessor-methods
#' @export
setGeneric('getOrganism', function (x) standardGeneric('getOrganism'))

#' @rdname accessor-methods
#' @export
setGeneric('getTaxonomy', function (x) standardGeneric('getTaxonomy'))

#' @rdname accessor-methods
#' @export
setGeneric('getKeywords', function (x) standardGeneric('getKeywords'))

#' @rdname accessor-methods
#' @export
setGeneric('getReference', function (x) standardGeneric('getReference'))

#' @rdname accessor-methods
#' @export
setGeneric('getComment', function (x) standardGeneric('getComment'))


# header, features, sequence ------------------------------------------------


#' Get the feature table from a GenBank record.
#'
#' @param x A \code{\linkS4class{gbRecord}} instance.
#' @param ... Additional arguments passed to methods.
#' @return The \code{\linkS4class{gbFeatureTable}} of a Genbank record.
#' @rdname getFeatures-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' getFeatures(x)
setGeneric("getFeatures", function(x, ...) standardGeneric("getFeatures"))

#' @rdname getFeatures-methods
#' @export
#' @examples
#' ft(x)
setGeneric("ft", function(x, ...) standardGeneric("ft"))

#' Get the sequence from a GenBank record.
#' 
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} instance.
#' @param ... Additional arguments passed to methods.
#' @return An \code{\linkS4class{XStringSet}} object, containing either the
#' complete sequence(s) of the record(s), or of the selected feature(s)
#' @rdname getSequence-methods
#' @export
#' @examples
#' \dontrun{
#' gbk_file <- system.file("extdata", "S_cerevisiae_mito.gbk", package = "biofiles")
#' x <- gbRecord(gbk_file)
#' 
#' ## extract the full-length sequence of the record.
#' getSequence(x)
#' 
#' ## extract coding sequences only
#' getSequence(x["CDS"])
#' }
setGeneric("getSequence", function(x, ...) standardGeneric("getSequence"))

#' Extract the header from a \code{"\linkS4class{gbRecord}"} object.
#' 
#' @param x A \code{"\linkS4class{gbRecord}"}, \code{"\linkS4class{gbFeature}"},
#'  or \code{"\linkS4class{gbFeatureTable}"} instance.
#' @param ... Additional arguments passed to methods.
#' @return A \code{"\linkS4class{gbHeader}"} instance
#' @rdname getHeader-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' getHeader(x)
setGeneric("getHeader", function(x, ...) standardGeneric("getHeader"))

#' @rdname getHeader-methods
#' @export
#' @examples
#' header(x)
setGeneric("header", function(x, ...) standardGeneric("header"))

#' Summarise a GenBank record. 
#'
#' @param object An object of class\code{\linkS4class{gbFeature}},
#' \code{\linkS4class{gbFeatureTable}}, \code{\linkS4class{gbRecord}}, or
#' \code{\linkS4class{gbRecordList}}.
#' @param ... Arguments to be passed to methods.
#' @rdname summary-methods
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' summary(x)
setGeneric("summary")

#' Access the indices of GenBank features
#'
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param ... Additional arguments passed to methods.
#' @return A numeric vector of feature indeces.
#' @rdname index-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' index(x)
setGeneric("index", signature = "x", function (x, ...) {
  standardGeneric("index")
})

#' Access genomic locations of GenBank features
#'
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}}, or
#' \code{\linkS4class{gbRecord}} object.
#' @param join Join compound genomic locations to a single range.
#' @param ... Additional arguments passed to methods.
#' @return A list of \code{\linkS4class{gbLocation}} objects
#' @rdname location-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' location(x)
setGeneric("location", signature = "x", function (x, ...) {
  standardGeneric("location")
})

#' Get/set keys of GenBank features
#'
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param ... Additional arguments passed to methods.
#' @rdname key-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' key(x)
setGeneric("key", signature = "x", function(x, ...) {
  standardGeneric("key")
})

#' @param check if \code{FALSE}, don't perform validity checks.
#' @param value The key information to set on \code{x}.
#' @rdname key-methods
#' @export
setGeneric("key<-", signature = "x", function(x, check = TRUE, value) {
  standardGeneric("key<-")
})

#' Get/set qualifiers of GenBank features
#' 
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}}, or
#' \code{\linkS4class{gbRecord}} object.
#' @param which A character vector giving the name(s) of the
#' qualifiers to retrieve or set.
#' @param fixed If \code{TRUE}, \code{which} is matched against qualifiers as is,
#' if \code{FALSE} it is treated as a regular expression.
#' @param use.names If \code{TRUE}, return a \code{data.frame} using \code{which}
#' as column names, if \code{FALSE} return, if possible, a character vector or
#' a list.
#' @param ... Additional arguments passed to methods.
#' @return A \code{data.frame}.
#' @rdname qualif-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' qualif(x[[1]], 'db_xref')
#' 
#' ## use shortcuts to common qualifiers
#' proteinID(x["CDS"])
#' locusTag(x["CDS"])
setGeneric("qualif", signature = "x", function(x, which = "", ...) {
  standardGeneric("qualif")
})

#' @param check if \code{FALSE}, don't perform validity checks.
#' @param value The qualifier information to set on \code{x}.
#' @rdname qualif-methods
#' @export
setGeneric("qualif<-", signature = "x", function(x, which, check = TRUE, value) {
  standardGeneric("qualif<-")
})

#' Access the \code{db_xref}s of GenBank features
#' 
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}}, or
#' \code{\linkS4class{gbRecord}} object.
#' @param db (Optional) A character vector giving the database names of the
#' desired \code{db_xref}s.
#' @param ... Additional arguments passed to methods.
#' @return A named character vector (or list of named character vectors)
#' of db_xrefs.
#' @rdname dbxref-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' 
#' ## all db_xrefs associated with CDSs
#' dbxref(x["CDS"])
#' 
#' ## retrieve the TaxId from the "source" field.
#' dbxref(x[[1]], "taxon")
setGeneric("dbxref", signature = "x", function(x, db = NULL, ...) {
  standardGeneric("dbxref")
})


# write/save ----------------------------------------------------------------


#' Save and load \code{gbRecord} objects.
#' 
#' Serialise and unserialise \code{\linkS4class{gbRecord}}s using
#' \code{\link{saveRDS}} and \code{\link{readRDS}}
#' 
#' @param x A \code{\linkS4class{gbRecord}} or \code{\linkS4class{gbRecordList}} instance.
#' @param file A character string naming the file to write to or read from.
#' If \code{NULL}, the accession number will be used to construct a file name.
#' @param dir Target directory. (Default: current working directory)
#' @param ... Arguments passed to \code{\link{saveRDS}}.
#' @rdname saveRecord-methods
#' @export
#' @examples
#' \dontrun{
#' aca <- genomeRecordFromNCBI("Bacteria/Acaryochloris_marina", verbose = TRUE)
#' aca
#' saveRecord(aca)
#' rm(aca)
#' aca <- loadRecord("./NC_009925_NC_009926_NC_009927_NC_009928_NC_009929_NC_0099__.rds")
#' aca
#' }
setGeneric("saveRecord", function(x, file = NULL, dir = ".", ...) standardGeneric("saveRecord"))


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
#' @param ... Additional arguments passed to methods.
#' @rdname write.GenBank-methods
#' @export
#' @examples
#' \dontrun{
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' write.GenBank(x, file = "data/marine_metagenome.gbk")
#' 
#' ## write selected features to file.
#' write.GenBank(x["CDS"], file = "data/marine_metagenome_cds.gbk", header = FALSE, sequence = FALSE)
#' }
setGeneric("write.GenBank", function(x, file, append = FALSE, ...) {
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
#' @param sequence if \code{TRUE}, additionally output a fasta file.
#' @param append if \code{TRUE} the data is appended to the connection.
#' @param ... Additional arguments passed to methods.
#' @rdname write.FeatureTable-methods
#' @export
#' @examples
#' \dontrun{
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' write.FeatureTable(x, file = "data/marine_metagenome.tbl")
#' }
setGeneric("write.FeatureTable", signature = "x",
           function(x, file, tablename = "", dbname = "",
                    sequence = FALSE, append = FALSE, ...) {
             standardGeneric("write.FeatureTable")
           })


# list-generics ----------------------------------------------------------


#' List the names of Genbank qualifiers.
#'
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecord}} object.
#' @param ... Additional arguments to be passed to or from methods.
#' @return A character vector (or list of character vectors) of 
#'    qualifier names.
#' @seealso
#'    \code{\link{uniqueQualifs}}, \code{\link{hasQualif}}
#' @rdname qualifList-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' qualifList(x["source"])
setGeneric("qualifList", function(x, ...) standardGeneric("qualifList"))


#' Tabulate Genbank qualifiers
#' 
#' Extract a frequency table (or list of tables in the case of
#' \code{gbRecordList}s) of qualifier names.
#'
#' @param x A \code{\linkS4class{gbFeatureTable}}, \code{\linkS4class{gbRecord}},
#' or \code{\linkS4class{gbRecordList}} object.
#' @param ... Additional arguments to be passed to or from methods.
#' @return A \code{\link{table}} (or list of \code{table}s) of 
#'    qualifiers names.
#' @seealso
#'    \code{\link{uniqueQualifs}}, \code{\link{hasQualif}}
#' @rdname qualifTable-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' qualifTable(x)
setGeneric("qualifTable", function(x, ...) standardGeneric("qualifTable"))


#' Tabulate Genbank features
#' 
#' Extract a frequency table (or list of tables in the case of
#' \code{gbRecordList}s) of feature keys.
#'
#' @inheritParams qualifTable
#' @return A \code{\link{table}} (or list of \code{table}s) of 
#'    feature keys.
#' @seealso
#'    \code{\link{qualifTable}}
#' @rdname featureTable-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' featureTable(x)
setGeneric("featureTable", function(x, ...) standardGeneric("featureTable"))



# test-generics ----------------------------------------------------------

#' Has a feature a specific key?
#'  
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param key A feature key. 
#' @param ... Additional arguments to be passed to or from methods.
#' @return A logical vector or a list of logical vectors.
#' @seealso
#'    \code{\link{key}}
#' @rdname hasKey-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' hasKey(x, 'CDS')
setGeneric("hasKey", signature = c("x","key"), function(x, key, ...) {
  standardGeneric("hasKey")
})

#' Has a feature a specific qualifier?
#'
#' @param x A \code{\linkS4class{gbFeature}}, \code{\linkS4class{gbFeatureTable}},
#' \code{\linkS4class{gbRecord}}, or \code{\linkS4class{gbRecordList}} object.
#' @param qualifier A character string. Name of a qualifier.
#' @param ... Additional arguments to be passed to or from methods.
#' @return A logical vector or a list of logical vectors.
#' @seealso
#'  \code{\link{qualifList}}, to extract a list of available qualifiers for
#'  each feature; \code{\link{uniqueQualifs}}, for a vector of all
#'  unique qualifiers present in an object. 
#' @rdname hasQualif-methods
#' @export
#' @examples
#' load(system.file("extdata", "marine_metagenome.rda", package = "biofiles"))
#' hasQualif(x, 'CDS')
setGeneric("hasQualif", signature = c("x", "qualifier"), function(x, qualifier, ...) {
  standardGeneric("hasQualif")
})


# shift and revcomp -------------------------------------------------------------


#' Shift the location of features in a GenBank record
#'
#' @note \code{shift} does not currently handle compound locations In a shifted
#' feature table compound locations get merged.
#' @usage shift(x, shift = 0L, split = FALSE, order = TRUE, ...)
#' @param x A \code{\linkS4class{gbFeatureTable}} or
#' \code{\linkS4class{gbRecord}} instance (gbFeatureTables must 
#' be complete and include a 'source' field).
#' @param shift Number of basepairs (or aa residues) to shift.
#' @param split Split features that after the shift extends across the end of
#' the sequence.
#' @param order Reorder features after the shift.
#' @param ... Additional arguments passed to methods.
#'
#' @return A \code{\linkS4class{gbFeatureTable}} object.
#' @rdname shift-methods
#' @export
#' @examples
#' load(system.file("extdata", "S_cerevisiae_mito.rda", package = "biofiles"))
#' 
#' ## shift the S. cerevisiae mitochondrion such that cytochrome b is the first CDS
#' cytb <- start(filter(x, product = "^cytochrome b$")[[1]])[1]
#' x2 <- shift(x, shift = -cytb + 1, split = TRUE)
setGeneric("shift", signature = "x", function(x, shift = 0L, split = FALSE, order = TRUE, ...) {
  standardGeneric("shift")
})

#' Reverse-complement features in a GenBank record
#' 
#' @usage revcomp(x, order = TRUE, ...)
#' @param x A \code{\linkS4class{gbFeatureTable}} or
#' \code{\linkS4class{gbRecord}} object (gbFeatureTables must 
#' be complete and include a 'source' field).
#' @param order Reorder features after reverse-complementing them.
#' @param ... Additional arguments passed to methods.
#' @rdname revcomp-methods
#' @export
#' @examples
#' load(system.file("extdata", "S_cerevisiae_mito.rda", package = "biofiles"))
#' xr <- revcomp(x)
setGeneric("revcomp", signature = "x", function(x, order = TRUE, ...) standardGeneric("revcomp"))


# view -------------------------------------------------------------------


#' View all features in a \code{gbFeatureTable}
#' 
#' @param x A \code{\linkS4class{gbFeatureTable}} instance.
#' @param n How many features to show (Default: all).
#' @param ... Additional arguments passed to methods.
#' @rdname view-methods
#' @keywords internal
#' @export
setGeneric("view", signature = "x", function(x, n, ...) {
  standardGeneric("view")
})


# filter --------------------------------------------------------------------


#' Return a subset of features or annotations from a GenBank Record
#' 
#' \code{filter} returns a subset of features from \code{\linkS4class{gbRecordList}},
#' \code{\linkS4class{gbRecord}} or \code{\linkS4class{gbFeatureTable}} objects,
#' based on filters provided as \emph{key}, \emph{range}, or \emph{qualifier} values.
#'
#' @details
#' Filters are provided as named values using keywords and/or
#' \dQuote{\emph{qualifier = value}} pairs:
#' 
#' Permissible keywords are:
#' 
#' \describe{
#'   \item{index/idx}{
#'     For example: \code{idx = c(3,4,5,6)}, \code{idx = 100:150},
#'     \code{index = c(1,12:20)}
#'   }
#'   \item{range}{
#'     For example: \code{range = "10000..25000"},
#'     \code{range = "..10000,20000..25000"},
#'     \code{range = "30000.."}
#'   }
#'   \item{key}{
#'     For example: \code{key = "CDS"}, \code{key = c("CDS", "gene")}
#'   }
#'   \item{arbitrary qualifiers}{
#'     For example: \code{product = "ribosomal"}, \code{locus_tag = c("CPSIT_0123",
#'     "CPSIT_0124", "CPSIT_0125")}, \code{pseudo = TRUE}
#'   }
#' }
#' 
#' 
#' @param x A \sQuote{\code{gbRecord}} or \sQuote{\code{gbFeatureTable}}
#' instance.
#' @param ... For \code{filter}: named values that specify the features to select. These are
#' merged with the values of \code{keys} to create the actual query. See
#' Details; for \code{select}: see \code{.cols}.
#' @param .cols A character vector of \sQuote{\emph{keys}} that specify annotaions
#' to be returned as a \code{data.frame} from the filtered features. If \code{NULL},
#' a \sQuote{\code{gbFeatureTable}} is returned.
#' Supported \sQuote{\emph{keys}} are \dQuote{index} or \dQuote{idx}, \dQuote{start},
#' tag (e.g., \dQuote{locus_tag}, \dQuote{product}, \dQuote{db_xref}). Specific
#' \code{db_xref}s can by queried using, e.g. \dQuote{db_xref.GI} or
#' \dQuote{db_xref.GeneID}.
#' 
#' @return Depending on the value of \code{.col} a \code{gbRecordList},
#' \code{gbRecord}, or\code{gbFeatureTable} or a \code{data.frame}.
#' @rdname manip-methods
#' @export
#' @examples
#' load(system.file("extdata", "S_cerevisiae_mito.rda", package = "biofiles"))
#' 
#' ## filter all hydrophobic tRNAs from the yeast mitochondrion
#' hydrophobic <- c("Val", "Ile", "Leu", "Met", "Phe", "Trp", "Cys")
#' trna <- filter(x, key = "tRNA", product = hydrophobic)
#' 
#' ## select start, end, orientation, product, and GeneID
#' df <- select(trna, "start", "end", "strand", "product", "db_xref.GeneID")
#' df
#' 
#' ## combine the above steps into one
#' cols <- c("start", "end", "strand", "product", "db_xref.GeneID")
#' filter(x, key = "tRNA", product = hydrophobic, .cols = cols)
#' 
#' ## filter all CDS from position 60,000 bp onward
#' filter(x, key = "CDS", range = "60000..")
setGeneric("filter", signature = "x", function(x, ...) standardGeneric("filter"))


# select -----------------------------------------------------------------


#' \code{select} returns a specified subset of annotations from GenBank Features
#' as a \code{data.frame}.
#' @rdname manip-methods
#' @export
setGeneric("select", signature = "x", function(x, ...) standardGeneric("select"))


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

