#' Quickly list all qualifier names
#' 
#' @usage listUniqueQualifs(x)
#' @param x A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeatureList}},
#' or, \code{\linkS4class{gbFeature}} instance
#' @return A character vector of qualifier names
#' @export
listUniqueQualifs <- Compose("unique", "unlist", "listQualif")


#' @usage locusTag(x)
#' @rdname qualif-methods
#' @export
locusTag <- Partial("qualif", which="locus_tag", use.names=FALSE)


#' @usage product(x)
#' @rdname qualif-methods
#' @export
product <- Partial("qualif", which="product", use.names=FALSE)


#' @usage note(x)
#' @rdname qualif-methods
#' @export
note <- Partial("qualif", which="note", use.names=FALSE)


#' @usage proteinID(x)
#' @rdname qualif-methods
#' @export
proteinID <- Partial("qualif", which="protein_id", use.names=FALSE)


#' @usage geneID(x)
#' @rdname qualif-methods
#' @export
geneID <- Partial("qualif", which="gene", use.names=FALSE)


.translation <- Partial("qualif", which="translation", use.names=FALSE)
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
  db <- switch(getMoltype(x), AA="protein", "nuccore")
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
      f <- efetch(acc[i], db, "fasta", "xml", seqstart=s[i], seqstop=e[i])
      if (str[i] == -1) {
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
    r <- IRanges(start=is, end=ie, names=names(dna))
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
      cat("Parsing:", strsplitN(txt[2], 'DEFINITION  ', 2), sep="\n")
    }
    parse_gb_record(txt)
  }
  list(
    update = update,
    value = value,
    record = record
  )
}


#' @importFrom RCurl getCurlHandle basicTextGatherer curlOptions
#' @importFrom RCurl curlSetOpt chunkToLineReader curlPerform
#' @importFrom iterators iter
#' @importFrom foreach foreach "%do%"
fetchGbRecords <- function(urls, ..., curl=getCurlHandle()) {
  content <- list()
  foreach(url = iter(urls)) %do% {
    content[[url]] <- gbReader(...)
    opts <- curlOptions(url=url, writefunc=chunkToLineReader(content[[url]]$update)$read)
    curl <- curlSetOpt(.opts=opts, curl=curl, ...)
    tryCatch(curlPerform(curl=curl), error = function(e) {
      warning(e$message, immediate.=TRUE)
    })
  }
  content
}


fetchHeader <- function(urls, ..., curl=getCurlHandle()) {
  content <- list()
  foreach(url = iter(urls)) %do% {
    content[[url]] <- basicTextGatherer()
    opts <- curlOptions(url=url, headerfunc=content[[url]]$update, 
                        writefunc=content[[url]]$update, nobody=TRUE,
                        header=TRUE)
    curl <- curlSetOpt(.opts=opts, curl=curl, ...)
    tryCatch(curlPerform(curl=curl), error = function(e) {
      warning(e$message, immediate.=TRUE)
    })
  }
  content
}


contentLength <- function(urls, unit="Mb") {
  unit <- match.arg(unit, c("b", "kb", "Mb", "Gb"))
  headers <- fetchHeader(urls)
  res <- foreach(url = iter(urls), .combine="c") %do% {
    h <- usplit(headers[[url]]$value(), '\n|\r\n')
    len <- as.numeric(strsplitN(unique(grep("Content-Length", h, value=TRUE)), ':\\s+', 2))
    len <- switch(unit, b=len, kb=(len/1024), Mb=(len/1024^2), Gb=(len/1024^3))
    setNames(len, basename(url))
  }
  attr(res, "unit") <- unit
  res
}


chooseFromList <- function(list, prompt = "", default = NULL) {
  while (TRUE) {
    which <- tolower(usplit(readline(prompt), '\\s+')) %||% default
    if (is.null(which)) {
      list <- NULL
      break
    } else if (!collapse(which, '') %in% c("a","all","q","quit")) {
      if (!any(is.na(idx <- suppressWarnings(as.numeric(which))))
          && max(idx) <= length(list)
          && min(idx) > 0) {
        list <- list[idx]
        break
      } else {
        message("Please enter 'all', 'quit', or numbers between 1 and ", length(list))
      }
    } else {
      if (which[1] == 'q' || which[1] == 'quit') {
        list <- NULL
      } 
      break
    }
  }
  list
}


#' Retrieve genomes in GenBank format directly form NCBI
#' 
#' ## EXPERIMENTAL ##
#' 
#' @param which Path to ftp directory on NCBI's
#' \href{ftp://ftp.ncbi.nih.gov/genomes/}{Genomes FTP site}. Examples would be
#' \code{Bacteria/Acetinobacter*} or \code{Viruses/Bat_coronavirus*}.
#' If there are multiple matching directories the user will be prompted to choose
#' one. If there are multiple matching \code{gbk} files the user will also be
#' prompted to choose one or more.
#' @param ignore.case Ignore case when matching.
#' @param ... Arguments passed on to \code{\link[RCurl]{curlOptions}}.
#' @return A \code{\linkS4class{gbRecord}} or \code{\linkS4class{gbRecordList}}
#' object.
#' @export
#' @examples
#' \dontrun{
#' gbk <- genomeRecordFromNCBI("Bacteria/Chlamydia_muridarum", verbose = TRUE)
#' 
#' }
#' 
genomeRecordFromNCBI <- function(which, ignore.case=TRUE, .parse=TRUE, ...) {
  if (missing(which)) {
    stop("\"which\" is missing with no default", call.=TRUE)
  }
  MAX <- function(x, n) if (length(x) > n) c(x[1:n], "...") else x
  base_url <- "ftp://ftp.ncbi.nih.gov/genomes/" 
  which <- usplit(which, "/", fixed=TRUE)
  
  g <- basicTextGatherer()
  curl <- getCurlHandle(header=FALSE, ftplistonly=TRUE, ftp.use.epsv=FALSE)

  while(length(which) > 0) {
    base_url <- paste0(base_url, which[1], '/')
    tryCatch(curlPerform(url=base_url, writefunction=g$update, curl=curl), error=function(e) {
      stop(e$message, call.=TRUE)
    })
    which <- which[-1]
    x <- usplit(g$value(), '\n')
    target <- x[grep(which[1] %|na|% '\\.gbk$', x, ignore.case=ignore.case)]
    
    if (length(target) > 1 && !all(grepl('\\.gbk$', target))) {
      cat(sprintf("[%s] %s", seq_along(target), target), sep="\n")
      msg <- paste0("Choose target directory (", collapse(MAX(seq_along(target), 4)), "|q(uit)]: ")
      target <- chooseFromList(target, msg, 1)
      if (length(target) > 1) {
        warning("Multiple target directories. Will use the first only")
        target <- target[1]
      }
      which[1] <- target
    } else if (length(target) == 1 && !all(grepl('\\.gbk$', target))) {
      which[1] <- target
    } else if (all(grepl('\\.gbk$', target))) {
      urls <- paste0(base_url, target)
      size <- contentLength(urls)
      cat(sprintf("[%s] %s (%s %s)", seq_along(target), target, round(size, 2),
                  attr(size, "unit")), sep="\n")
      msg <- paste0("Choose gbk file(s) to download: [a(ll)|",
                    collapse(MAX(seq_along(target), 4)), "|q(uit)]: ")
      target <- chooseFromList(target, msg)
      urls <- paste0(base_url, target)
      break
    } else {
      stop("Something went horribly wrong")
    }
    g$reset()
  }
  curl <- curlSetOpt(ftplistonly=FALSE, curl=curl)
  res <- fetchGbRecords(urls, ..., curl=curl)
  gbkList <- list()
  for (gb in names(res)) {
    nm <- paste0(strsplitN(dirname(gb), "/", 1, 'end', fixed=TRUE), '/', basename(gb))
    gbkList[[nm]] <- if (.parse) {
      res[[gb]]$record() 
    } else {
      res[[gb]]$value()
    }
  }
  if (length(gbkList) == 1) {
    gbkList[[1]]
  } else {
    if (.parse) {
      gbRecordList(gbkList)
    } else {
      gbkList
    }
  }
}


