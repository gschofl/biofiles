#' @include gbLocation-class.R
NULL

#' gbInfo
#' 
#' \dQuote{gbInfo} is an S4 class that provides basic infomation
#' about a genomic sequence.
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


#' @export
gbInfo <- function (object) {
  if (is(object, "gbRecord") && isValidDb(object)) {
    new("gbInfo", db = object, 
        seqnames = dbFetch(object, 'accession'),
        seqlengths = as.integer(dbFetch(object, 'length')),
        is_circular = dbFetch(object, 'topology') == "circular",
        genome = dbFetch(object, 'definition'))
  } else if (!isS4(object) && hasValidDb(object)) {
    db_dir <- attr(object, "dir")
    db <- new("gbRecord", dir=normalizePath(db_dir),
              name=basename(db_dir), verbose=FALSE)
    new("gbInfo", db = db, 
        seqnames = dbFetch(db, 'accession'),
        seqlengths = as.integer(dbFetch(db, 'length')),
        is_circular = dbFetch(db, 'topology') == "circular",
        genome = dbFetch(db, 'definition'))
  } else {
    stop("Need a valid 'gbRecord' instance to initialize 'gbInfo'")
  }
}


setValidity2("gbInfo", function (object) {
  if (!is(object@db, "gbRecord")) {
    return(FALSE)
  }
  if (!isValidDb(object@db@dir)) {
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


.make_GRanges <- function (x, join = FALSE, include = "none", exclude = "", key = TRUE) {
  if (!hasValidDb(x)) {
    stop("No valid gbFeature")
  }
  if (key) {
    qual <- DataFrame(key = key(x))
  } else {
    qual <- DataFrame()
  }
  if (any(include != "none")) {
    if (all(include == "all")) {
      which <- setdiff(listUniqueQualifs(x), exclude)
    } else {
      which <- setdiff(include, exclude)
    }
    qual <- c(qual, DataFrame(.simplify(.qualAccess(x, which, fixed=TRUE),
                                        unlist=FALSE)))
  }
  
  start <- start(x, join = join)
  end <- end(x, join = join)
  strand <- strand(x, join = join)
  lt <- .qualAccess(x, 'locus_tag', fixed=TRUE)
  gene <- .qualAccess(x, 'gene', fixed=TRUE)
  names <- ifelse(is.na(lt), gene, lt)
    
  if (is.list(start)) {
    i <- update_indeces(start)
    start <- unlist(start)
    end <- unlist(end)
    strand <- unlist(strand)
    names <- names[i]
    if (length(qual) > 0)
      qual <- qual[i, , drop=FALSE] 
  }

  GRanges(seqnames=Rle(accession(x)), ranges=IRanges(start, end, names=names),
          strand=strand, qual, seqinfo = seqinfo(x))
}


update_indeces <- function(x) {
  j_len <- vapply(x, length, numeric(1))
  j_idx <- which(j_len > 1)
  j_len <- j_len[j_idx]
  sort(c(seq_along(x), rep(j_idx, j_len - 1)))
}
