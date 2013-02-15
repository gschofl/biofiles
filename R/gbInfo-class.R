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


.make_GRanges <- function (x, join = TRUE, with_qual = "none", without_qual = "") {
  if (!hasValidDb(x)) {
    stop("Is no valid gbFeature")
  }  
  
  qual <- DataFrame(key = key(x))
  
  if (any(with_qual != "none")) {
    if (all(with_qual == "all")) {
      which <- setdiff(listUniqueQualifs(x), without_qual)
    } else {
      which <- setdiff(with_qual, without_qual)
    }
    qual <- cbind(qual, DataFrame(as.list(qualif(x, which))))
  }

  GRanges(seqnames=Rle(accession(x)),
          ranges=IRanges(start(x, join=join), end(x, join=join), names=index(x)),
          strand=strand(x, join = join), qual, seqinfo = seqinfo(x))
}

