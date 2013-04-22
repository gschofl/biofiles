#' @autoImports
.make_GRanges <- function (x, join = FALSE, include = "none", exclude = "", key = TRUE) {
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
