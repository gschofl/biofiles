#' @include parser-general.R
NULL

.gbk_mandatory <- c("LOCUS", "DEFINITION", "ACCESSION", "VERSION", "FEATURES", "//")

#' Parser for GenBank/GenPept records.
#' 
#' @param x A character vector
#' @return A \code{\linkS4class{gbRecord}} instance.
#' @keywords internal
gbk_record <- function(rec) {
  # get a vector with the positions of the main GenBank fields
  rec_idx <- grep("^[A-Z//]+", rec)
  rec_kwd <- strsplitN(rec[rec_idx], " +", 1L)
  gbk_contig <- gbk_sequence <-  NULL
  
  # Check the presence of mandatory fields
  if (any(is.na(charmatch(.gbk_mandatory, rec_kwd)))) {
    stop("mandatory fields are missing from the GenBank file")
  }
  
  ## get positions of features, origin, contig and end_of_record
  ftb_idx <- rec_idx[rec_kwd == "FEATURES"]
  seq_idx <- rec_idx[rec_kwd == "ORIGIN"]
  ctg_idx <- rec_idx[rec_kwd == "CONTIG"]
  end_idx <- rec_idx[rec_kwd == "//"]
  ftb_end_idx <- rec_idx[which(rec_kwd == "FEATURES") + 1] - 1
  
  ## HEADER
  x <- rec[seq.int(ftb_idx - 1)]
  seqenv <- seqinfo(gbk_header(x), NULL)
  
  ## SEQUENCE
  if (length(seq_idx) > 0L) {
    # if "//" is right after "ORIGIN" there is no sequence
    # and gb_sequence stays set to NULL
    if (end_idx - seq_idx > 1L) {
      gbk_sequence <- rec[seq.int(seq_idx + 1, end_idx - 1)]
    }
    ## CONTIG
  } else if (length(ctg_idx) > 0L) {
    contig_line <- strsplitN(collapse(rec[seq.int(ctg_idx, end_idx-1)], ''),
                             'CONTIG', 2L, fixed = TRUE)
    gb_contig <- gbLocation(contig_line)
  }
  seqenv$sequence <- 
    parse_sequence(seq = gbk_sequence, acc = getAccession(seqenv),
                   seqtype = getMoltype(seqenv), src = "gbk")
  
  ## FEATURES
  ft <- rec[seq.int(ftb_idx + 1, ftb_end_idx)]
  ft <- parse_features(x = ft, seqinfo = seqenv)
  new_gbRecord(seqinfo = seqenv, features = ft, contig = gbk_contig) 
}

#' @keywords internal
gbk_header <- function(x) {
  # generate a vector with the positions of the main GenBank keywords.
  # keywords are all capitals, beginning in column 1 of a record.
  gbk_idx <- grep("^[A-Z//]+", x)
  gbk_kwd <- strsplitN(x[gbk_idx], split = " +", 1)

  ## LOCUS (Mandatory)
  locus <- gbk_locus(x[gbk_idx[gbk_kwd == "LOCUS"]])

  ## DEFINITION (Mandatory)
  def_idx <- which(gbk_kwd == "DEFINITION")
  def_line <- x[seq.int(gbk_idx[def_idx], gbk_idx[def_idx + 1] - 1)]
  definition <- collapse(sub("DEFINITION  ", "", def_line), ' ')

  ## ACCESSION (Mandatory)
  acc_line <- x[gbk_idx[gbk_kwd == "ACCESSION"]]
  accession <- strsplitN(acc_line, split = "\\s+", 2L)

  ## VERSION and GI (Mandatory)
  ver_line <- x[gbk_idx[gbk_kwd == "VERSION"]]
  version  <- usplit(ver_line, split = "\\s+")[2L]
  seqid    <- paste0('gi|', usplit(ver_line, split = "GI:", fixed = TRUE)[2L])

  ## DBLINK (Optional)
  if (length(db_line <- x[gbk_idx[gbk_kwd == "DBLINK"]]) > 0L) {
    dblink <- usplit(db_line, split = "Project: ", fixed = TRUE)[2L]
  } else {
    dblink <- NA_character_
  }

  ## DBSOURCE (GenPept only; sometimes more than one line)
  if (length(dbsrc_idx <- which(gbk_kwd == "DBSOURCE")) > 0L) {
    dbs_lines <- x[seq.int(gbk_idx[dbsrc_idx], gbk_idx[dbsrc_idx + 1] - 1)]
    dbsource <- collapse(gsub("^ +", "", sub("DBSOURCE", "", dbs_lines)), "\n")
  } else {
    dbsource <- NA_character_
  }

  ## KEYWORDS (Mandatory)
  key_line <- x[gbk_idx[gbk_kwd == "KEYWORDS"]]
  keywords <- sub("KEYWORDS    ", "", key_line)

  ## SOURCE with ORGANISM and the complete lineage (Mandatory)
  src_idx <- which(gbk_kwd == 'SOURCE')
  source_lines <- x[seq.int(gbk_idx[src_idx], gbk_idx[src_idx + 1] - 1)]                  
  source <- sub("SOURCE      ", "", source_lines[1L])
  organism <- sub("  ORGANISM  ", "", source_lines[2L])
  taxonomy <- collapse(gsub("^ +", "", source_lines[-c(1L, 2L)]), ' ')

  ## REFERENCES (Mandatory?)
  if (length(ref_idx <- which(gbk_kwd == "REFERENCE")) > 0L) {
    ref_lines <-
      x[
        seq.int(
          gbk_idx[ref_idx[1]],
          (gbk_idx[ref_idx[length(ref_idx)] + 1] - 1) %|na|% length(x)
        )]
    references <- gbk_reference_list(ref_lines)
  } else {
    references <- .gbReferenceList()
  }

  ## COMMENT (Optional)
  if (length(gbk_idx[gbk_kwd == "COMMENT"]) > 0L) {
    com_lines <- x[seq.int(min(gbk_idx[gbk_kwd == "COMMENT"]), length(x))]
    comment <- collapse(gsub("^ +", "", sub("COMMENT", "", com_lines)), "\n")
  } else {
    comment <- NA_character_
  }
  
  .gbHeader(
    locus = locus,
    definition = definition,
    accession = accession,
    version = version,
    seqid = seqid,
    dblink = dblink,
    dbsource = dbsource,
    keywords = keywords,
    source = source,
    organism = organism,
    taxonomy = taxonomy,
    references = references,
    comment = comment
  )
}

#' @keywords internal
gbk_locus <- function(locus_line) {
  tokens <- usplit(locus_line, split = "\\s+")[-1]
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

#' @keywords internal
gbk_reference <- function(ref) {
  ## split by subkeywords
  ref_idx <- grep("^ {0,3}[A-Z]+", ref)
  ref_list <- ixsplit(ref, ref_idx, include_i = TRUE, collapse_x = TRUE)
  ## 
  kwd   <- vapply(ref_list, strsplitN, '\\s+', 1L, FUN.VALUE = "")
  field <- vapply(ref_list, strsplitN, '^[A-Z]+\\s+(?!\\S)\\s', 2L, perl = TRUE, FUN.VALUE = "")
  ##
  ref <- set_reference()
  ref$refline(field[kwd == "REFERENCE"])
  ref$authors(field[kwd == "AUTHORS"])
  ref$consrtm(field[kwd == "CONSRTM"])
  ref$title(field[kwd == "TITLE"])
  ref$journal(field[kwd == "JOURNAL"])
  ref$pubmed(field[kwd == "PUBMED"])
  ref$remark(field[kwd == "REMARK"])
  ref$yield()
}

#' @keywords internal
gbk_reference_list <- function(ref_lines) {
  ## split references
  ref_idx <- grep("REFERENCE", ref_lines, fixed = TRUE, ignore.case = FALSE)
  ref_list <- ixsplit(ref_lines, ref_idx)
  .gbReferenceList(ref = lapply(ref_list, gbk_reference))
}

