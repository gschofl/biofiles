#' @include parser-general.R
NULL

.embl_mandatory <- c("ID", "AC", "DT", "DE", "KW", "OS", "RN", "FH")

#' Parser for Embl Records.
#' 
#' @param x A character vector
#' @return A \code{\linkS4class{gbRecord}} instance.
#' @keywords internal
embl_record <- function(rec) {
  # get a vector with the positions of the main Embl fields by grepping the
  # spacers (XX) and adding one to the index
  rec_idx <- c(1, grep("^XX", rec) + 1L)
  rec_kwd <- strsplitN(rec[rec_idx], "\\s+", 1L)
  embl_contig <- embl_sequence <-  NULL
  
  # Check the presence of mandatory fields
  if (any(is.na(charmatch(.embl_mandatory, rec_kwd)))) {
    stop("mandatory fields are missing from the Embl file")
  }
  
  ## get positions of features, sequence, contig and end_of_record
  ftb_idx <- rec_idx[rec_kwd == "FH"] + 2L
  seq_idx <- rec_idx[rec_kwd == "SQ"] + 1L
  ctg_idx <- rec_idx[rec_kwd == "CO"]
  end_idx <- grep("^//", rec)
  ftb_end_idx <- rec_idx[which(rec_kwd == "FH") + 1] - 1
  
  ## HEADER
  x <- rec[seq.int(ftb_idx - 3)]
  seqenv <- seqinfo(embl_header(x), NULL)
  
  ## SEQUENCE
  if (length(seq_idx) > 0L) {
    # if "//" is right after "ORIGIN" there is no sequence
    # and gb_sequence stays set to NULL
    if (end_idx - ori_idx > 1L) {
      embl_sequence <- rec[seq.int(seq_idx + 1, end_idx - 1)]
    }
    ## CONTIG
  } else if (length(ctg_idx) > 0L) {
    contig_line <- strsplitN(collapse(rec[seq.int(ctg_idx, end_idx - 1)], ''), '^CO\\s+', 2L)
    embl_contig <- gbLocation(contig_line)
  }
  seqenv$sequence <-
    parse_sequence(seq = embl_sequence, acc = getAccession(seqenv)[1],
                   seqtype = getMoltype(seqenv), src = "embl")
  
  ## FEATURES
  ft <- rec[seq.int(ftb_idx, ftb_end_idx - 1)]
  ft <- parse_features(x = ft, seqinfo = seqenv)
  new_gbRecord(seqinfo = seqenv, features = ft, contig = embl_contig) 
}

#' @keywords internal
embl_header <- function(x) {
  ## remove all the XX lines
  x <- compactXX(x)
  # generate a vector of the main Embl keywords.
  kwd <- strsplitN(x, split = "\\s+", 1)
  
  ## ID [identification] and DT [date ] (Mandatory)
  locus_line <- embl_line(x, kwd, 'ID')
  date_line <- embl_line(x, kwd, 'DT')
  locus <- embl_locus(locus_line, date_line)
  
  ## AC [accession number] (Mandatory)
  acc_line <- embl_line(x, kwd, 'AC', ' ')
  accession <- usplit(acc_line, split = ";\\s*")
  
  ## DE [description] (Mandatory)
  definition <- embl_line(x, kwd, 'DE', ' ')
  
  ## SV [sequence version] (Optional)
  version <- embl_line(x, kwd, 'SV') %||% NA_character_
  seqid    <- NA_character_

  ## PR [project identifier] (Optional)
  dblink <- if (length(db_line <- embl_line(x, kwd, "PR")) > 0L) {
    usplit(usplit(db_line, split = "Project: ?")[2L], ";\\s*")
  } else {
    NA_character_
  }
  
  ## DR [Database cross-reference] (Optional)
  dbsource <- embl_line(x, kwd, "DR", ' ')

  ## KW [keyword] (Mandatory)
  keywords <- embl_line(x, kwd, "KW", ' ')
  
  ## OS [organism species]; OC [organism classification] (Mandatory)             
  source <- embl_line(x, kwd, "OS")
  organism <- sub(" \\(.+\\)$", "", source)
  taxonomy <- embl_line(x, kwd, "OC", ' ')
  
  ## REFERENCES (Mandatory)
  references <- if (length(ref_idx <- which(kwd == "RN")) > 0L) {
    ref_lines <- x[grepl('^R', x)]
    embl_reference_list(ref_lines)
  } else {
    .gbReferenceList()
  }
  
  ## CC [comments or notes] (Optional)
  comment <- embl_line(x, kwd, 'CC', ' ') %||% NA_character_

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
embl_locus <- function(locus_line, date_line) {
  ## The tokens represent:
  ## 1. Primary accession number
  ## 2. Sequence version number (Currently we don't use SV)
  ## 3. Topology: 'circular' or 'linear'
  ## 4. Molecule type
  ## 5. Data class
  ## 6. Taxonomic division (not always present)
  ## 7. Sequence length
  tokens <- usplit(locus_line, split = ";\\s+")
  ntok <- length(tokens)
  seqlen <- usplit(tokens[ntok], " ")
  dates <- usplit(date_line, "\\s.+")
  .gbLocus(
    lnm = tokens[1],
    len = as.integer(seqlen[1]),
    mtp = if (seqlen[2] != 'BP.') 'AA' else tokens[4],
    top = tokens[3],
    div = if (ntok > 6) paste0(tokens[5], '/', tokens[6]) else tokens[5], ## combine data class and taxonomic division
    cdt = dates[1],
    mdt = dates[2]
  )
}

#' @keywords internal
embl_reference <- function(x) {
  kwd <- strsplitN(x, split = "\\s+", 1)
  ## Reference cross-reference
  rx <- strsplit(embl_line(x, kwd, "RX"), ";\\s+") 
  crossref <- vapply(rx, function(x) sub(".$", "", x[2]), FUN.VALUE = "")
  names(crossref) <- vapply(rx, "[", 1, FUN.VALUE = "")
  crossref <- c(crossref["PUBMED"], crossref[names(crossref) != "PUBMED"])
  ref <- set_reference()
  ref$refline( trim(embl_line(x, kwd, "RN"), "\\[|\\]") ) ## Reference Number
  ref$authors( embl_line(x, kwd, "RA", ' ') )             ## Reference Author
  ref$consrtm( embl_line(x, kwd, "RG", ' ') )             ## Reference Group
  ref$title( embl_line(x, kwd, "RT", ' ') )               ## Reference Title
  ref$journal( embl_line(x, kwd, "RL", ' ') )             ## Reference Location
  ref$pubmed( crossref )                                  ## Reference cross-reference
  ref$remark( embl_line(x, kwd, "RC", ' ') )              ## Reference Comment
  ref$yield()
}

#' @keywords internal
embl_reference_list <- function(ref_lines) {
  ## split references
  ref_idx <- grep("^RN", ref_lines)
  ref_list <- ixsplit(ref_lines, ref_idx)
  .gbReferenceList(ref = lapply(ref_list, embl_reference))
}

