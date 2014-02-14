#' Quickly list all qualifier names
#' 
#' @usage uniqueQualifs(x)
#' @param x A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeatureTable}},
#' or, \code{\linkS4class{gbFeature}} instance
#' @return A character vector of qualifier names
#' @export
uniqueQualifs <- Compose("unique", "unlist", "qualifList")


#' @usage locusTag(x)
#' @rdname qualif-methods
#' @export
locusTag <- Partial("qualif", which = "locus_tag", use.names = FALSE)


#' @usage product(x)
#' @rdname qualif-methods
#' @export
product <- Partial("qualif", which = "product", use.names = FALSE)


#' @usage note(x)
#' @rdname qualif-methods
#' @export
note <- Partial("qualif", which = "note", use.names = FALSE)


#' @usage proteinID(x)
#' @rdname qualif-methods
#' @export
proteinID <- Partial("qualif", which = "protein_id", use.names = FALSE)


#' @usage geneID(x)
#' @rdname qualif-methods
#' @export
geneID <- Partial("qualif", which = "gene", use.names = FALSE)


.translation <- Partial("qualif", which = "translation", use.names = FALSE)
#' @usage translation(x)
#' @rdname qualif-methods
#' @importFrom Biostrings AAStringSet
#' @export
translation <- function(x) AAStringSet(.translation(x))


#' Retrieve the sequence of a contig
#' 
#' ## EXPERIMENTAL ##
#' 
#' @param x gbRecord
#' @param merge
#' @importFrom Biostrings unlist DNAStringSet
#' @importFrom IRanges metadata "metadata<-"
#' @importFrom reutils efetch
#' @export
getContigSeq <- function(x, merge = TRUE) { 
  db <- switch(getMoltype(x), AA = "protein", "nuccore")
  contig <- .contig(x)
  s <- start(contig)
  e <- end(contig)
  w <- width(contig)
  str <- strand(contig)
  acc <- getAccession(contig)
  dna <- DNAStringSet()
  for (i in seq_along(acc)) {
    if (!nzchar(acc[i])) {
      dna <- c(dna, DNAStringSet(dup('N', w[i])))
      dna@ranges@NAMES[i] <- paste0('Gap:', w[i])
    } else {
      f <- efetch(acc[i], db, "fasta", "xml", seqstart = s[i], seqstop = e[i])
      if (str[i]  == -1) {
        dna <- c(dna, reverseComplement(DNAStringSet(f$xmlValue("//TSeq_sequence"))))
      } else {
        dna <- c(dna, DNAStringSet(f$xmlValue("//TSeq_sequence")))
      }
      dna@ranges@NAMES[i] <- paste0('Acc:', f$xmlValue("//TSeq_accver"),
                                    ';GI:', f$xmlValue("//TSeq_gi"),
                                    ';SID:', f$xmlValue("//TSeq_sid"),
                                    ';TaxId:', f$xmlValue("//TSeq_taxid"),
                                    ';defline:', f$xmlValue("//TSeq_defline"))
    }
  }
  
  if (merge) {
    is <- c(1, cumsum(width(dna)) + 1)
    is <- is[-length(is)]
    ie <- cumsum(width(dna))
    r <- IRanges(start = is, end = ie, names = names(dna))
    res <- DNAStringSet(unlist(dna))
    metadata(res) <- list(ranges = r)
    return(res)
  }
  
  dna
}


gbReader <- function(verbose = FALSE) {
  txt <- character()
  update <- function(str) {
    con <- textConnection(str)
    on.exit(close(con))
    txt <<- c(txt, readLines(con))
  }
  value <- function() {
    txt
  }
  record <- function() {
    if (verbose) {
      cat("Parsing:", strsplitN(txt[2], 'DEFINITION  ', 2), sep = "\n")
    }
    parse_gb_record(txt)
  }
  list(
    update = update,
    value = value,
    record = record
  )
}
