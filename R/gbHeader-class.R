#### Layout of the entry fields in a GenBank Record ####
#
# LOCUS - Name for the entry. Mandatory keyword.
# 
# DEFINITION - A concise description of the sequence. Mandatory keyword.
# 
# ACCESSION  - The primary accession number. Mandatory keyword.
# 
# VERSION - The primary accession number and a numeric version number,
# followed by an integer key (GI) assigned to the sequence. Mandatory keyword.
# 
# DBLINK - Cross-references to resources that support the existence a sequence
# record (Project Database, NCBI Trace Assembly Archive). Optional keyword.
#
# DBSOURCE - GenPept only.
# 
# KEYWORDS - Other information about an entry. Mandatory keyword.
# 
# SEGMENT	- The order in which this entry appears in a series of discontinuous
# sequences from the same molecule. Optional keyword.
# 
# SOURCE	- Common name of the organism or the name most frequently used
# in the literature. Mandatory keyword. One Subkeyword.
#
#   ORGANISM - Formal scientific name of the organism.
#            - Taxonomic classification levels.
#
# REFERENCE - Citations for all articles containing the data reported
# in this entry. Mandatory keyword. Six subkeywords. 
# 
#   AUTHORS	- Authors of the citation. Optional subkeyword.
# 
#   CONSRTM	- The collective names of consortia. Optional subkeyword.
# 
#   TITLE	- Full title of citation. Optional subkeyword.
# 
#   JOURNAL	- The journal name, volume, year, and page numbers of the citation.
#   Mandatory subkeyword.
# 
#   PUBMED - The PubMed unique identifier for a citation. Optional subkeyword.
# 
#   REMARK	- The relevance of a citation. Optional subkeyword.
#
# COMMENT  - Cross-references to other sequence entries, comparisons to
# other collections, notes of changes in LOCUS names, and other remarks.
# Optional keyword.
# 
# FEATURES - Feature Table. Optional keyword.
# 
# CONTIG - Information about how individual sequence records can be combined to
# form larger-scale biological objects, such as chromosomes or complete genomes.
# A special join() statement on the CONTIG line provides the accession numbers
# and basepair ranges of the underlying records which comprise the object.
#
# ORIGIN - Start of sequence data. Mandatory keyword.
# 
# // 	- Entry termination. Mandatory at the end of an entry.
#


#### Layout of a LOCUS record ####
#
#
# Keyword LOCUS, locus name (1) sequence length (2), bp or aa (3),
# if 'bp' molecule type (4), topology (5, optional), GenBank division (6),
# and modification date (7)
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
    is_empty = function() {
      sum(length(lnm), length(len), length(mtp), length(top),
          length(div), length(cdt), length(mdt)) != 7
    },
    to_string = function() {
      type <- if (mtp == "AA") "aa" else "bp"
      smtp <- ifelse(substring(mtp, 3, 3)=='-', mtp, paste0(blanks(3), mtp %|AA|% blanks(2)))
      sprintf("%-12s%-17s %+10s %s %-10s %-8s %s %s",
              "LOCUS", lnm, len, type, smtp, top %|NA|% '', div,
              toupper(format(mdt, "%d-%b-%Y")))
    },
    show = function() {
      if (is_empty()) {
        showme <- sprintf("An empty %s instance.\n", sQuote(class(.self)))
      } else {
        showme <- paste0(
          sprintf("A %s instance:\n", sQuote(class(.self))),
          ellipsize(to_string())
          )
      }
      cat(showme, "\n")
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


gbReferences <- function(ref_lines) {
  return("Not implemented yet")
  ref_idx <- grep("^REFERENCE", ref_lines)
  idx <- Map(seq, ref_idx, c(ref_idx[-1] - 1, length(ref_lines)))
  ref <- ref_lines[idx[[1]]]
  auth_line <- grep("^  AUTHORS", ref, value=TRUE)
  .parseAuthors <- function(auth_line) {}
}


#' Generator object for the \code{\linkS4class{gbHeader}} class
#'
#' The generator object for the \code{\linkS4class{gbHeader}} reference class.
#'
#' @param ... List of arguments; must be named arguments
#' corresponding to the fields of a \code{\linkS4class{gbHeader}} object
#' @section Methods:
#' \describe{
#' \item{\code{#new(...)}:}{
#'    Create a new \code{\linkS4class{gbHeader}} object}
#' }
#' 
#' @seealso
#'    \code{\linkS4class{gbHeader}}
#' @rdname gbHeader
#' @keywords classes internal
#' @export
.gbHeader <- setRefClass(
  'gbHeader',
  fields = list(
    locus = 'gbLocus',
    definition = 'character',
    accession = 'character',
    version = 'character',
    seqid = 'character', ## NCBI GI identifier
    dblink = 'character',
    dbsource = 'character', ## GenPept only
    keywords = 'character',
    source = 'character',
    organism = 'character',
    taxonomy = 'character',
    references = 'character',
    comment = 'character'
  ),
  methods = list(
    is_empty = function() {
      locus$is_empty() && sum(length(definition), length(accession)) != 2
    },
    to_string = function(write_to_file = FALSE) {
      'Generate a character string representation of a GenBank file header'
      if (write_to_file) {
        loc <- locus$to_string()
        w <- 79
        f <- FALSE
      } else {
        loc <- ellipsize(locus$to_string())
        w <- getOption("width") - 2
        f <- TRUE
      }
      o <- i <- 12
      paste0(
        loc,
        sprintf("\n%-12s%s\n", "DEFINITION", linebreak(definition, width=w, indent=0, offset=o, FORCE=f)),
        sprintf("%-12s%s\n", "ACCESSION", accession),
        sprintf("%-12s%-12s%s%s\n", "VERSION", version, "GI:", strsplitN(seqid, "|", 2, fixed = TRUE)),
        if (!all(is.na(dblink))) {
          sprintf("%-12s%s%s\n", "DBLINK", "Project: ", dblink)
        } else {
          ""
        },
        sprintf("%-12s%s\n", "KEYWORDS", linebreak(keywords, width=w, offset=o, FORCE=f)),
        sprintf("%-12s%s\n", "SOURCE", linebreak(source, width=w, offset=o, FORCE=f)),
        sprintf("%-12s%s\n", "  ORGANISM",
                paste0(organism, '\n', linebreak(taxonomy, width=w, indent=i, offset=o, FORCE=f))),
        sprintf("%-12s%-3s(bases %s to %s)\n%-12s%s\n%-12s%s\n%-12s%s\n",
                "REFERENCE", 1, 1, locus$len,
                "  AUTHORS", "authors",
                "  TITLE", "title",
                "  JOURNAL", "journal"),
        if (!is.na(comment)) {
          sprintf("%-12s%s\n", "COMMENTS", linebreak(comment, width=w, offset=o, FORCE=f))
        } else {
          ""
        })
    },
    show = function() {
      'Method for automatically printing a Genbank file header.'
      if (is_empty()) {
        showme <- sprintf("An empty %s instance.\n", sQuote(class(.self)))
      } else {
        showme <- paste0(
          sprintf("A %s instance:\n", sQuote(class(.self))),
          to_string(write_to_file = FALSE)
        )
      }
      cat(showme, "\n")
    },
    write = function(file = "", append = FALSE, sep = "\n") {
      'Write a GenBank file header to file.'
      cat(to_string(write_to_file = TRUE), file = file, append = append, sep = sep)
    }
  )
)


#' Class \code{"gbHeader"}
#'
#' A container for GenBank file headers. 
#' @name gbHeader-class
#' @section Fields:
#' \describe{
#' \item{\code{locus}:}{ A \code{\linkS4class{gbLocus}} object. }
#' \item{\code{definition}:}{ \code{character}; Description of the sequence. }
#' \item{\code{accession}:}{ \code{character}; The primary accession number.
#'      A unique, unchanging identifier assigned to each GenBank sequence record. }
#' \item{\code{version}:}{ \code{character}; The primary accession number and a
#'      numeric version number }
#' \item{\code{seqid}:}{ \code{character}; Gene Identifier ("GI"); An integer key
#'      assigned to the sequence by NCBI. }
#' \item{\code{dblink}:}{ \code{character}; Cross-references to resources that
#'      support the existence a sequence record, such as the Project Database
#'      and the NCBI Trace Assembly Archive. }
#' \item{\code{dbsource}:}{  \code{character}; GenPept files only }
#' \item{\code{keywords}:}{ code{character}; Short description of gene products
#'      and other information about an entry. }
#' \item{\code{source}:}{ \code{character}; Common name of the organism. }
#' \item{\code{organism}:}{ \code{character}; Formal scientific name of the
#'      organism. }
#' \item{\code{taxonomy}:}{ \code{character}; Taxonomic classification levels. }
#' \item{\code{references}:}{ Citations for all articles containing data
#'      reported in the entry. Yet to be implemented. }
#' \item{\code{comment}:}{ \code{character}; Remarks. }
#' }
#' @section Extends: All reference classes extend and inherit methods from
#'    \code{"\linkS4class{envRefClass}"}.
#' @seealso
#'    \code{\link{.gbHeader}}
#' @keywords classes
#' @examples
#'
#' showClass("gbHeader")
#'
NULL


#' @keywords internal
gbHeader <- function(gb_header) {
  # generate a vector with the positions of the main GenBank keywords.
  # keywords are all capitals, beginning in column 1 of a record.
  gbk_idx <- grep("^[A-Z//]+", gb_header)
  gbk_kwd <- strsplitN(gb_header[gbk_idx], split=" +", 1)
  ##
  ## LOCUS (Mandatory)
  locus <- gbLocus(gb_header[gbk_idx[gbk_kwd == "LOCUS"]])
  ##
  ## DEFINITION (Mandatory)
  def_idx <- which(gbk_kwd == "DEFINITION")
  def_line <- gb_header[seq.int(gbk_idx[def_idx], gbk_idx[def_idx + 1] - 1)]
  definition <- paste0(trim(sub("DEFINITION  ", "", def_line)), collapse=" ")
  ##
  ## ACCESSION (Mandatory)
  acc_line <- gb_header[gbk_idx[gbk_kwd == "ACCESSION"]]
  accession <- strsplitN(acc_line, split="\\s+", 2L)
  ##
  ## VERSION and GI (Mandatory)
  ver_line <- gb_header[gbk_idx[gbk_kwd == "VERSION"]]
  version  <- usplit(ver_line, split="\\s+")[2L]
  seqid    <- paste0('gi|', usplit(ver_line, split="GI:", fixed=TRUE)[2L])
  ##
  ## DBLINK (Optional)
  if (length(db_line <- gb_header[gbk_idx[gbk_kwd == "DBLINK"]]) > 0L) {
    dblink <- usplit(db_line, split="Project: ", fixed=TRUE)[2L]
  } else {
    dblink <- NA_character_
  }
  ##
  ## DBSOURCE (GenPept only; sometimes more than one line)
  if (length(dbsrc_idx <- which(gbk_kwd == "DBSOURCE")) > 0L) {
    dbs_lines <- gb_header[seq.int(gbk_idx[dbsrc_idx], gbk_idx[dbsrc_idx + 1] - 1)]
    dbsource <- paste(gsub("^ +", "", sub("DBSOURCE", "", dbs_lines)), collapse="\n")
  } else {
    dbsource <- NA_character_
  }
  ##
  ## KEYWORDS (Mandatory)
  key_line <- gb_header[gbk_idx[gbk_kwd == "KEYWORDS"]]
  keywords <- sub("KEYWORDS    ", "", key_line)
  ##
  ## SOURCE with ORGANISM and the complete lineage (Mandatory)
  src_idx <- which(gbk_kwd == 'SOURCE')
  source_lines <- gb_header[seq.int(gbk_idx[src_idx], gbk_idx[src_idx + 1] - 1)]                  
  source <- sub("SOURCE      ", "", source_lines[1L])
  organism <- sub("  ORGANISM  ", "", source_lines[2L])
  taxonomy <- paste(gsub("^ +", "", source_lines[-c(1L, 2L)]), collapse=" ")
  ##
  ## REFERENCES (Mandatory?)
  if (length(ref_idx <- which(gbk_kwd == "REFERENCE")) > 0L) {
    ref_lines <-
      gb_header[
        seq.int(
          gbk_idx[ref_idx[1]],
          (gbk_idx[ref_idx[length(ref_idx)] + 1] - 1) %|NA|% length(gb_header)
        )]
    references <- gbReferences(ref_lines)
  } else {
    references <- "Not available"
  }
  ##
  ## COMMENT (Optional)
  if (length(gbk_idx[gbk_kwd == "COMMENT"]) > 0L) {
    com_lines <- gb_header[seq.int(gbk_idx[gbk_kwd == "COMMENT"], length(gb_header))]
    comment <- paste(gsub("^ +", "", sub("COMMENT", "", com_lines)), collapse="\n")
  } else {
    comment <- NA_character_
  }
  .gbHeader(
    locus=locus,
    definition=definition,
    accession=accession,
    version=version,
    seqid=seqid,
    dblink=dblink,
    dbsource=dbsource,
    keywords=keywords,
    source=source,
    organism=organism,
    taxonomy=taxonomy,
    references=references,
    comment=comment
  )
} 


setOldClass("environment")

#' seqinfo-class
#' 
#' @name seqinfo-class
#' @rdname seqinfo-class
#' @exportClass seqinfo
setClass("seqinfo", contains="environment")


setMethod('initialize', 'seqinfo', function(.Object) {
  .Object@.xData <- new.env(parent=emptyenv())
  .Object
})

## Internal Getters

setMethod('.header', 'seqinfo', function(x) {
  tryCatch(get("header", x), error = function(e) {
    warning("No header associated with this object", call.=FALSE)
    .gbHeader$new()
  })
})

setMethod('.sequence', 'seqinfo', function(x) {
  tryCatch(get("sequence", x), error = function(e) {
    warning("No sequence associated with this object", call.=FALSE)
    new('BStringSet')
  })
})

setMethod('.locus', 'seqinfo', function(x) .header(x)$locus)

## Getters

setMethod("getLocus", "seqinfo", function(x) .locus(x)$lnm)


setMethod("getLength", "seqinfo", function(x) .locus(x)$len)


setMethod("getMoltype", "seqinfo", function(x) .locus(x)$mtp)


setMethod("getTopology", "seqinfo", function(x) .locus(x)$top)


setMethod("getDivision", "seqinfo", function(x) .locus(x)$div)


setMethod("getDate", "seqinfo", function(x) {
  c(create_date = .locus(x)$cdt,
    update_date = .locus(x)$mdt)
})


setMethod("getDefinition", "seqinfo", function(x) .header(x)$definition)


setMethod("getAccession", "seqinfo", function(x) .header(x)$accession)


setMethod("getVersion", "seqinfo", function(x) .header(x)$version)


setMethod("getGeneID", "seqinfo", 
          function(x, db = 'gi') {
            seqid <- .header(x)$seqid
            if (is.na(seqid)) {
              seqid
            } else {
              db.idx <- which(strsplitN(seqid, "|", 1, fixed = TRUE) == db)
              strsplitN(seqid, "|", 2, fixed = TRUE)[db.idx]
            }
          })


setMethod("getDBLink", "seqinfo", function(x) .header(x)$dblink)


setMethod("getDBSource", "seqinfo", function(x) .header(x)$dbsource)


setMethod("getSource", "seqinfo", function(x) .header(x)$source)


setMethod("getOrganism", "seqinfo", function(x) .header(x)$organism)


setMethod("getTaxonomy", "seqinfo", function(x) .header(x)$taxonomy)


setMethod("getReference", "seqinfo", function(x) .header(x)$references)


setMethod("getKeywords", "seqinfo", function(x) .header(x)$keywords)


setMethod("getComment", "seqinfo", function(x) .header(x)$comment)


setMethod("show", "seqinfo", function(object) {
  if (length(getAccession(object)) == 0L) {
    acc <- len <- def <- ""
  } else {
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






