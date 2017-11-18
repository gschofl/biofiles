.GRanges <- function(x, join = FALSE, include = "none", exclude = "", key = TRUE) {
  if (length(x) == 0) {
    return(GenomicRanges::GRanges())
  }
  if (key) {
    qual <- S4Vectors::DataFrame(key = key(x))
  } else {
    qual <- S4Vectors::DataFrame()
  }
  if (any(include != "none")) {
    if (all(include == "all")) {
      which <- setdiff(uniqueQualifs(x), exclude)
    } else {
      which <- setdiff(include, exclude)
    }
    qual <- c(qual, S4Vectors::DataFrame(.simplify(.qual_access(x, which, fixed = TRUE),
                                        unlist = FALSE)))
  }
  
  my_start  <- start(x, join = join)
  my_end    <- end(x, join = join)
  my_strand <- strand(x, join = join)
  lt <- .qual_access(x, which = 'locus_tag', fixed = TRUE)
  gene <- .qual_access(x, 'gene', fixed = TRUE)
  names <- ifelse(is.na(lt), gene, lt)
 
  if (is.list(my_start)) {
    i <- update_indices(my_start)
    my_start <- unlist(my_start)
    my_end <- unlist(my_end)
    my_strand <- unlist(my_strand)
    names <- names[i]
    if (length(qual) > 0)
      qual <- qual[i, , drop = FALSE] 
  }

  seqinfo <- GenomeInfoDb::Seqinfo(
    seqnames = getAccession(x),
    seqlengths = getLength(x),
    isCircular = getTopology(x) == 'circular',
    genome = getDefinition(x)
  )
  
  GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(getAccession(x)),
    ranges   = IRanges::IRanges(my_start, my_end, names = names),
    strand   = S4Vectors::Rle(my_strand),
    qual,
    seqinfo  = seqinfo
  )
}


.IRanges <- function(x) {
  if (length(x) == 0) {
    return(IRanges::IRanges())
  }
  jr <- matrix(joint_range(x), ncol = 2)
  IRanges::IRanges(start = jr[, 1], end = jr[, 2], names = index(x))
}


update_indices <- function(x) {
  j_len <- vapply(x, length, numeric(1))
  j_idx <- which(j_len > 1)
  j_len <- j_len[j_idx]
  sort(c(seq_along(x), rep(j_idx, j_len - 1)))
}
