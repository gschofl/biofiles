#' Quickly list all qualifier names.
#' 
#' @usage uniqueQualifs(...)
#' @param ... A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeatureTable}},
#' or, \code{\linkS4class{gbFeature}} instance.
#' @return A character vector of qualifier names.
#' @export
uniqueQualifs <- Compose("unique", "unlist", "qualifList")

#' Return the \emph{locus_tag} qualifiers from GenBank features.
#' 
#' @usage locusTag(...)
#' @param ... A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeatureTable}},
#' or, \code{\linkS4class{gbFeature}} instance.
#' @return A character vector of \emph{locus_tag}s.
#' @export
locusTag <- Partial("qualif", which = "locus_tag", use.names = FALSE)

#' Return the \emph{product} qualifiers from GenBank features.
#' 
#' @usage product(...)
#' @param ... A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeatureTable}},
#' or, \code{\linkS4class{gbFeature}} instance.
#' @return A character vector of \emph{product}s.
#' @export
product <- Partial("qualif", which = "product", use.names = FALSE)

#' Return the \emph{note} qualifiers from GenBank features.
#' 
#' @usage note(...)
#' @param ... A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeatureTable}},
#' or, \code{\linkS4class{gbFeature}} instance.
#' @return A character vector of \emph{note}s.
#' @export
note <- Partial("qualif", which = "note", use.names = FALSE)

#' Return the \emph{protein_id} qualifiers from GenBank features.
#' 
#' @usage proteinID(...)
#' @param ... A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeatureTable}},
#' or, \code{\linkS4class{gbFeature}} instance.
#' @return A character vector of \emph{protein_id}s.
#' @export
proteinID <- Partial("qualif", which = "protein_id", use.names = FALSE)

#' Return the \emph{gene} qualifiers from GenBank features.
#' 
#' @usage geneID(...)
#' @param ... A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeatureTable}},
#' or, \code{\linkS4class{gbFeature}} instance.
#' @return A character vector of \emph{gene}s.
#' @export
geneID <- Partial("qualif", which = "gene", use.names = FALSE)


.translation <- Partial("qualif", which = "translation", use.names = FALSE)
#' Return the \emph{translation}s from GenBank features.
#' 
#' @usage translation(x)
#' @param x A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeatureTable}},
#' or, \code{\linkS4class{gbFeature}} instance.
#' @return An \code{\linkS4class{AAStringSet}}.
#' @importFrom Biostrings AAStringSet
#' @export
translation <- function(x) AAStringSet(.translation(x))


#' Retrieve the sequence of a contig
#' 
#' ## EXPERIMENTAL ##
#' 
#' For GenBank records that contain a CONTIG field, try to get the contig
#' sequences and (optionally) merge them into a single sequence.
#' 
#' @param x A \code{\linkS4class{gbRecord}} object.
#' @param merge Merge the retrieved contig sequences.
#' @return A \code{\linkS4class{DNAStringSet}} instance.
#' @importFrom Biostrings unlist DNAStringSet
#' @importFrom S4Vectors "metadata<-" width
#' @importFrom reutils efetch
#' @export
getContigSeq <- function(x, merge = TRUE) { 
  db <- switch(getMoltype(x), AA = "protein", "nuccore")
  contig <- .contig(x)
  s <- start(contig)
  e <- end(contig)
  w <- span(contig)
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
    parse_gbk_record(txt)
  }
  list(
    update = update,
    value = value,
    record = record
  )
}
