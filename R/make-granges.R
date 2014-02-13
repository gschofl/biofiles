#' @importFrom IRanges DataFrame Rle
#' @importFrom GenomicRanges GRanges GRangesList Seqinfo
NULL

.GRanges <- function(x, join = FALSE, include = "none", exclude = "", key = TRUE) {
  if (length(x) == 0) {
    return(GRanges())
  }
  if (key) {
    qual <- DataFrame(key = key(x))
  } else {
    qual <- DataFrame()
  }
  if (any(include != "none")) {
    if (all(include == "all")) {
      which <- setdiff(uniqueQualifs(x), exclude)
    } else {
      which <- setdiff(include, exclude)
    }
    qual <- c(qual, DataFrame(.simplify(.qual_access(x, which, fixed=TRUE),
                                        unlist=FALSE)))
  }
  
  start <- start(x, join = join)
  end <- end(x, join = join)
  strand <- strand(x, join = join)
  lt <- .qual_access(x, 'locus_tag', fixed = TRUE)
  gene <- .qual_access(x, 'gene', fixed = TRUE)
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

  seqinfo <- Seqinfo(seqnames=getAccession(x),
                     seqlengths=getLength(x),
                     isCircular=getTopology(x)=='circular',
                     genome=getDefinition(x))
  
  GRanges(seqnames=Rle(getAccession(x)), ranges=IRanges(start, end, names=names),
          strand=strand, qual, seqinfo = seqinfo)
}


.IRanges <- function(x) {
  if (length(x) == 0) {
    return(IRanges())
  }
  jr <- joint_range(x)
  IRanges(start = jr[, 1], end = jr[, 2], names = index(x))
}


update_indeces <- function(x) {
  j_len <- vapply(x, length, numeric(1))
  j_idx <- which(j_len > 1)
  j_len <- j_len[j_idx]
  sort(c(seq_along(x), rep(j_idx, j_len - 1)))
}
