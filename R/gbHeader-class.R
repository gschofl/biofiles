setOldClass("environment")

# The items in a LOCUS record are: 
# locus name (1) sequence length (2), bp or aa (3), if 'bp' molecule type (4),
# topology (5, optional), GenBank division (6), and modification date (7)
#
# Legal values for the division code include:
#   
# PRI - primate sequences
# ROD - rodent sequences
# MAM - other mammalian sequences
# VRT - other vertebrate sequences
# INV - invertebrate sequences
# PLN - plant, fungal, and algal sequences
# BCT - bacterial sequences
# VRL - viral sequences
# PHG - bacteriophage sequences
# SYN - synthetic sequences
# UNA - unannotated sequences
# EST - EST sequences (Expressed Sequence Tags) 
# PAT - patent sequences
# STS - STS sequences (Sequence Tagged Sites) 
# GSS - GSS sequences (Genome Survey Sequences) 
# HTG - HTGS sequences (High Throughput Genomic sequences) 
# HTC - HTC sequences (High Throughput cDNA sequences) 
# ENV - Environmental sampling sequences
# CON - Constructed sequences
# TSA - Transcriptome Shotgun Assembly sequences
#
#' Generator object for the \code{\linkS4class{gbLocus}} class
#'
#' The generator object for the \code{\linkS4class{gbLocus}} reference class.
#'
#' @param ... List of arguments (see NOTE)
#' @section Methods:
#' \describe{
#' \item{\code{#new(lnm, len, mtp, div, top, mdt, cdt)}:}{
#'    Create a new \code{\linkS4class{gbLocus}} object}
#' }
#' 
#' @note Arguments to the \code{#new} method must be named arguments:
#' \itemize{
#' \item{lnm}{ Locus name; stored in the \code{lnm} field. }
#' \item{len}{ Sequence lenght; stored in the \code{len} field. }
#' \item{mtp}{ Molecule type; stored in the \code{mtp} field). }
#' \item{div}{ Genbank division; stored in the \code{div} field. }
#' \item{top}{ Topology; stored in the \code{top} field. }
#' \item{mdt}{ Modification date; stored in the \code{mdt} field. }
#' \item{cdt}{ Create date; stored in the \code{cdt} field. }
#' } 
#'  
#' @seealso
#'    \code{\linkS4class{gbLocus}}
#' @rdname gbLocus
#' @keywords classes internal
#' @export
.gbLocus <- setRefClass(
  'gbLocus',
  fields = list(
    lnm = 'character',
    len = 'integer',
    mtp = 'character',
    top = 'character',
    div = 'character',
    cdt = 'POSIXlt',
    mdt = 'POSIXlt'
  ),
  methods = list(
    initialize = function(lnm, len, mtp, top, div, cdt, mdt) {
      if (!nargs()) return()
      lnm <<- lnm
      len <<- as.integer(len)
      mtp <<- mtp
      top <<- top
      div <<- div
      cdt <<- as.POSIXlt(cdt, format="%d-%b-%Y") 
      mdt <<- as.POSIXlt(mdt, format="%d-%b-%Y")
    },
    to_string = function() {
      type <- if (mtp == "AA") "aa" else "bp"
      smtp <- ifelse(substring(mtp, 3, 3)=='-', mtp, paste0(blanks(3), mtp %|AA|% blanks(2)))
      sprintf("%-12s%-17s %+10s %s %-10s %-8s %s %s",
              "LOCUS", lnm, len, type, smtp, top %|NA|% '', div,
              toupper(format(mdt, "%d-%b-%Y")))
    },
    show = function() {
      methods::show(ellipsize(gsub('\\s+', ' ', to_string())))
    }
  )
)


#' Class \code{"gbLocus"}
#'
#' A container for GenBank LOCUS records 
#' @name gbLocus-class
#' @section Fields:
#' \describe{
#' \item{\code{lnm}:}{ Locus name. Usually the accession number.}
#' \item{\code{len}:}{ Sequence length; In bp or aa, depending on \code{mtp}. }
#' \item{\code{mtp}:}{ Molecule type; NA, DNA, RNA, tRNA (transfer RNA), rRNA (ribosomal RNA), 
#'  mRNA (messenger RNA), uRNA (small nuclear RNA), or AA (protein sequence). RNAs
#'  can be prefixes ss- (single-stranded), ds- (double-stranded), or ms- (mixed-stranded)}
#' \item{\code{div}:}{ Genbank division. }
#' \item{\code{top}:}{ Topology; linear, circular, or missing (\code{NA}). }
#' \item{\code{mdt}:}{ Modification date. }
#' \item{\code{cdt}:}{ Create date. }
#' }
#' @section Extends: All reference classes extend and inherit methods from
#'    \code{"\linkS4class{envRefClass}"}.
#' @seealso
#'    \code{\link{.gbLocus}}
#' @keywords classes
#' @examples
#'
#' showClass("gbLocus")
#'
NULL


#' @keywords internal
gbLocus <- function(locus_line) {
  tokens <- usplit(locus_line, split="\\s+")[-1]
  # GenBank format: 'bp', GenPept: 'aa'
  gb <- if (tokens[3] == 'bp') TRUE else FALSE 
  date_idx <- length(tokens)
  divi_idx <- date_idx - 1
  topo_idx <- date_idx - 2
  if (gb && date_idx < 7 || !gb && date_idx < 6) {
    # topology is missing
    topo_idx <- NULL
  }
  .gbLocus(
    lnm = tokens[1],
    len = tokens[2],
    mtp = if (gb) tokens[4] else 'AA',
    top = tokens[topo_idx] %||% NA_character_,
    div = tokens[divi_idx] %||% NA_character_,
    cdt = tokens[date_idx],
    mdt = tokens[date_idx]
  )
}


#' gbHeader-class
#' 
#' Metadata about a sequence record.
#' 
#' @slot locus \code{gbLocus}; locus name, sequence length, molecule type,
#'      topology (optional), division code, modification date.
#' @slot definition \code{character}; Description of the sequence.
#' @slot accession \code{character}; The primary accession number. A unique, 
#'      unchanging identifier assigned to each GenBank sequence record.
#' @slot version \code{character}; The primary accession number and a numeric
#'      version number
#' @slot seqid \code{character}; Gene Identifier ("GI"); An integer key assigned
#'      to the sequence by NCBI.
#' @slot dblink \code{character}; Cross-references to resources that
#'      support the existence a sequence record, such as the Project Database
#'      and the NCBI Trace Assembly Archive.
#' @slot dbsource
#' @slot keywords \code{character}; Short description of gene products and other
#'      information about an entry.
#' @slot source \code{character}; Common name of the organism.
#' @slot organism \code{character}; Formal scientific name of the organism.
#' @slot taxonomy \code{character}; Taxonomic classification levels.
#' @slot references \code{character}; Citations for all articles containing data
#'      reported in the entry.
#' @slot comment \code{character}; Remarks.
#' 
#' @name gbHeader-class
#' @rdname gbHeader-class
#' @exportClass gbHeader
# .gbHeader <- setRefClass(
#   'gbHeader',
#   fields = list(
#       locus = 'gbLocus',
#       definition = 'character',
#       accession
#     )
#   )
setClass("gbHeader",
         slots = list(
           locus = "gbLocus", definition = "character", accession = "character",
           version = "character", seqid = "character", dblink = "character",
           dbsource = "character", keywords = "character", source = "character",
           organism = "character", taxonomy = "character", references = "character",
           comment = "character"
         ),
         prototype = list(
           definition = NA_character_, accession = NA_character_,
           version = NA_character_, seqid = NA_character_, dblink = NA_character_,
           dbsource = NA_character_, keywords = NA_character_, source = NA_character_,
           organism = NA_character_, taxonomy = NA_character_, references = NA_character_,
           comment = NA_character_
         ) 
)


#' seqinfo-class
#' 
#' @name seqinfo-class
#' @rdname seqinfo-class
#' @exportClass seqinfo
setClass("seqinfo", contains="environment")


setMethod('initialize', 'seqinfo', function (.Object) {
  .Object@.xData <- new.env(parent=emptyenv())
  .Object@.xData$header <- new("gbHeader")
  .Object@.xData$sequence <- new("DNAStringSet")
  .Object
})


setMethod('.header', 'seqinfo', function (x) get('header', x))


setMethod('.sequence', 'seqinfo', function (x) {
  tryCatch(get("sequence", x), error = function (e) {
    warning("No sequence associated with this object", call.=FALSE)
    new('BStringSet')
  })
})

setMethod("getLocus", "seqinfo", function (x) .header(x)@locus$lnm)


setMethod("getLength", "seqinfo", function (x) .header(x)@locus$len)


setMethod("getMoltype", "seqinfo", function (x) .header(x)@locus$mtp)


setMethod("getTopology", "seqinfo", function (x) .header(x)@locus$top)


setMethod("getDivision", "seqinfo", function (x) .header(x)@locus$div)


setMethod("getDate", "seqinfo", function (x) {
  c(create_date = .header(x)@locus$cdt,
    update_date = .header(x)@locus$mdt)
})


setMethod("getDefinition", "seqinfo", function (x) .header(x)@definition)


setMethod("getAccession", "seqinfo", function (x) .header(x)@accession)


setMethod("getVersion", "seqinfo", function (x) .header(x)@version)


setMethod("getGeneID", "seqinfo", 
          function (x, db = 'gi') {
            seqid <- .header(x)@seqid
            if (is.na(seqid)) {
              seqid
            }
            else {
              db.idx <- which(strsplitN(seqid, "|", 1, fixed = TRUE) == db)
              strsplitN(seqid, "|", 2, fixed = TRUE)[db.idx]
            }
          })


setMethod("getDBLink", "seqinfo", function (x) .header(x)@dblink)


setMethod("getDBSource", "seqinfo", function (x) .header(x)@dbsource)


setMethod("getSource", "seqinfo", function (x) .header(x)@source)


setMethod("getOrganism", "seqinfo", function (x) .header(x)@organism)


setMethod("getTaxonomy", "seqinfo", function (x) .header(x)@taxonomy)


setMethod("getReference", "seqinfo", function (x) .header(x)@references)


setMethod("getKeywords", "seqinfo", function (x) .header(x)@keywords)


setMethod("getComment", "seqinfo", function (x) .header(x)@comment)


setMethod("getSequence", "seqinfo", function (x) .sequence(x))


setMethod("show", "seqinfo", function (object) {
  if (is.na(getAccession(object))) {
    acc <- len <- def <- ""
  }
  else {
    acc <- getAccession(object)
    len <- paste0(getLength(object), " ", getMoltype(object))
    def <- getDefinition(object)
    acc <- pad(acc, nchar(acc) + 2, "right")
    len <- pad(len, nchar(len) + 2, "right")
    def <- ellipsize(def, width=getOption("width") - 
                     nchar(acc) - nchar(len) - 3)
  }
  cat(sprintf("%s%s%s", acc, len, def))
})






