#' @include gbLocation-class.R
NULL

#' gbInfo
#' 
#' \dQuote{gbInfo} is an S4 class that provides basic infomation
#' about a genomic sequence. It extends the \linkS4class{Seqinfo}.
#' 
#' @slot db A \code{\linkS4class{filehashRDS}} instance.
#' @slot seqnames Identifier of the sequence. Usually the Accession number.
#' @slot seqlengths Length of the sequence.
#' @slot is_circular Is the sequence circular.
#' @slot genome The definition of the sequence data.       
#' 
#' @rdname gbInfo
#' @export
#' @classHierarchy
#' @classMethods
setClass(Class="gbInfo", representation(db = "filehashRDS"),
         contains="Seqinfo")


gbInfo <- function (object) {
  if (is(object, "gbRecord") && isValidDb(object)) {
    new("gbInfo", db = object, 
        seqnames = dbFetch(object, 'accession'),
        seqlengths = as.integer(dbFetch(object, 'length')),
        is_circular = dbFetch(object, 'topology') == "circular",
        genome = dbFetch(object, 'definition'))
  } else {
    stop("Need a valid 'gbRecord' instance to initialize 'gbInfo'")
  }
}


setValidity2("gbInfo", function (object) {
  if (!is(object@db, "gbRecord")) {
    return(FALSE)
  }
  if (!isValidDb(object@db)) {
    return(FALSE)
  }
  if (dbFetch(object@db, 'accession') != object@seqnames) {
    return(FALSE)
  }
  if (dbFetch(object@db, 'length') != object@seqlengths) {
    return(FALSE)
  }
  if ((dbFetch(object@db, 'topology') == "circular") != object@is_circular) {
    return(FALSE)
  }
  if (dbFetch(object@db, 'definition') != object@genome) {
    return(FALSE)
  }
  
  return(TRUE)
})


createGRanges <- function (x, include_qual = "none", exclude_qual = "") {
  if (!hasValidDb(x)) {
    stop("Is no valid gbRecord")
  }  
  db <- init_db(x@.Dir, verbose=FALSE)
  seqid <- dbFetch(db, "accession")
  seqinfo <- Seqinfo(seqnames = seqid,
                     seqlengths = unname(dbFetch(db, "length")),
                     isCircular = dbFetch(db, "topology") == "circular",
                     genome = dbFetch(db, "definition"))
  
  qual <- DataFrame(key = key(x))
  if (include_qual != "none") {
    which <- if (include_qual == "all") {
      setdiff(listUniqueQualifs(x), exclude_qual)
    } else {
      which <- setdiff(include_qual, exclude_qual)
    }
    qual <- cbind(qual, DataFrame(qualif(x, which, attributes=FALSE)))
  }

  GRanges(seqnames=Rle(seqid),
          ranges=IRanges(start = start(x), end = end(x), names = index(x)),
          strand=Rle(strand(x)), qual, seqinfo = seqinfo)
}


setMethod("sequence", "GRanges",
          function (x, seq, ...) {
            
            if (is(seq, "gbRecord")) {
              seq <- sequence(seq)
            }
            
            if (is(seq, "DNAStringSet")) {
              if (length(seq) > 1) {
                warning("'seq' contains multiple sequences. Only the first will be used")
              }
              seq <- seq[[1]]
            }
            
            if (!is(seq, "DNAString")) {
              stop("'seq' must be a 'DNAString' object")
            }
            
            start <- biofiles::start(x)
            end <- biofiles::end(x)
            strand <- biofiles::strand(x)
            seqs <- as(Views(seq, start, end), "DNAStringSet")
            o <- c(seq_along(seqs)[strand == -1],
                   seq_along(seqs)[strand == 1])
            append(reverseComplement(seqs[strand == -1]),
                   seqs[strand == 1])[order(o)]
          })



