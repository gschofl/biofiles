recycle <- function (val, len) {
  lv <- length(val)
  if (len > lv) {
    val <- c(rep(val, len%/%lv), val[seq_len(len%%lv)])
  }
  val
}


#' @autoImports
showInfo <- function(object) {
  if (length(object) == 0) {
    seqname <- seqlen <- genome <- ""
  } else {
    sn <- seqnames(object)
    sl <- paste(seqlengths(object), names(seqlengths(object)))
    g <- genome(object)
    seqname <- pad(sn, nchar(sn) + 2, "right")
    seqlen <- pad(sl, nchar(sl) + 2, "right")
    genome <- ellipsize(g, width=getOption("width") - 
                          nchar(seqname) - nchar(seqlen) - 3)
  }
  cat(sprintf("%s%s%s", seqname, seqlen, genome))
}


#' @autoImports
ellipsize <- function(obj, width = getOption("width"), ellipsis = "...") {
  str <- encodeString(obj)
  ifelse(nchar(str) > width - 1,
         paste0(substring(str, 1, width - nchar(ellipsis) - 1), ellipsis),
         str)
}


#' @autoImports
is_compound <- function (x) {
  if (is(x, "gbFeatureList")) {
    return(vapply(x, function (f) !is.na(f@location@compound), logical(1)))
  } else if (is(x, "gbFeature")) {
    return(!is.na(x@location@compound))
  } else if (is(x, "gbLocation")) {
    return(!is.na(x@compound))
  }
}


#' @autoImports
getCompounds <- function (x) {
  x <- x[which(is_compound(x))]
  if (length(x) == 0) return(NA_real_) 
  cL <- vapply(x, function (f) nrow(f@location@range), integer(1))
  cL
}


#' @autoImports
expandIds <- function (x) {
  cmp_pos <- Position(is_compound, x)
  if (is.na(cmp_pos)) {
    return(list(ids=index(x, FALSE), keys=key(x, FALSE)))
  }
  xHead <- x[seq_len(cmp_pos - 1)]
  xCmp <- x[cmp_pos]
  L <- length(x)
  xTail <- if (cmp_pos + 1 <= L) {
    x[seq.int(cmp_pos +  1, L)]
  } else {
    new('gbFeatureList')
  }
  id1 <- list(ids=index(xHead, FALSE), keys=key(xHead, FALSE))
  nC <- getCompounds(xCmp)
  id2 <- list(ids=rep(index(xCmp, FALSE), nC),
              keys=rep(key(xCmp, FALSE), nC))
  id3 <- list(ids=index(xTail, FALSE), keys=key(xTail, FALSE))
  Map(function(a, b, c) c(a, b, c), a=id1, b=id2, c=Recall(xTail))
}


#' @autoImports
.qualAccess <- function (x, which = "", fixed = FALSE) {
  dbxrefs <- NULL
  dbx <- grepl('db_xref:.+', which)
  if (any(dbx)) {
    dbxrefs <- strsplitN(which[dbx], ':', 2)
    which <- c(which[!dbx], 'db_xref')
  }
  if (fixed) {
    which <- wrap(which, "\\b")
  }
  if (is(x, "gbFeature")) {
    .access(x, which, dbxrefs)
  } else if (is(x, "gbFeatureList")) {
    lapply(x, .access, which, dbxrefs)
  }
}


#' @autoImports
.access <- function (x, which, dbxrefs) {
  q <- x@qualifiers
  n <- length(q)
  els <- c(which[which != 'db_xref' & which != '\\bdb_xref\\b'], dbxrefs)
  
  if (n == 0) {
    return(setNames(rep(NA_character_, length(els)),
                    nm=rmisc::trim(els, "\\\\b")))
  }
  
  if (length(which) == 1) {
    idx <- grepl(which, names(q))
    if (any(idx))
      ans <- q[idx]
    else
      ans <- setNames(rep(NA_character_, length(els)),
                      nm=rmisc::trim(which, "\\\\b"))
  } else {
    idx <- lapply(which, grepl, names(q))
    ans <- lapply(idx, function(i) q[i])
    na <- which(vapply(ans, length, numeric(1)) == 0)
    if (length(na) > 0) {
      for (i in na) {
        ans[[i]] <- setNames(NA_character_, nm=trim(which, "\\\\b")[i])
      }
    }
  }
  ans <- unlist(ans)
  if (length(dbxrefs) > 0) {
    dbx_ <- names(ans) == 'db_xref'
    dbx <- strsplit(ans[dbx_], ':')
    dbx_nm <- vapply(dbx, `[`, 1, FUN.VALUE=character(1))
    dbx_val <- vapply(dbx, `[`, 2, FUN.VALUE=character(1))
    ans <- c(ans[!dbx_], setNames(dbx_val[match(dbxrefs, dbx_nm)], nm=dbxrefs))
  }
  
  return(ans)
}


#' @autoImports
.simplify <- function (x, unlist = TRUE) {
  if (length(len <- unique(unlist(lapply(x, length)))) > 1L) {
    return(x)
  }
  if (len == 1L && unlist) {
    unlist(x, recursive=FALSE)
  } else if (len >= 1L) {
    n <- length(x)
    r <- as.vector(unlist(x, recursive=FALSE))
    if (prod(d <- c(len, n)) == length(r)) {
      return(data.frame(stringsAsFactors=FALSE,
                        matrix(r, nrow=n, byrow=TRUE,
                               dimnames=if (!(is.null(nm <- names(x[[1L]]))))
                                 list(NULL, nm))))
    } else {
      x
    }
  } else {
    x
  }
}


#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings RNAStringSet
#' @importFrom Biostrings BStringSet
#' @importFrom Biostrings xscat
#' @autoImports
.seqAccess <- function (x) {
  
  if (exists("sequence", envir=x@.seqinfo))
    seq <- get("sequence", x@.seqinfo)
  else {
    warning("No sequence associated with this feature", call.=FALSE)
    return(BStringSet())
  }

  if (length(seq) == 0) {
    warning("No sequence associated with this feature", call.=FALSE)
    if (is(seq, "XStringSet"))
      return(seq)
    else
      return(BStringSet())
  }

  SEQF <- match.fun(class(seq))
  if (is(x, "gbFeature")) {
    seq <- merge_seq(seq, x, SEQF)
  } else if (is(x, "gbFeatureList")) {
    seq <- Reduce(append, lapply(x, merge_seq, seq=seq, SEQF=SEQF))
  }
  
  seq@metadata <- list(seqinfo(x))
  seq
}

# merge Sequences
#' @autoImports
merge_seq <- function (seq, x, SEQF) {
  if (length(start(x)) == 1L) {
    outseq <- subseq(x=seq, start=start(x), end=end(x))
  } else {
    outseq <- do.call(xscat, Map(subseq, x=seq, start=start(x), end=end(x)))
  }
  outseq <- SEQF(outseq)
  outseq@ranges@NAMES <- sprintf("%s.%s.%s", accession(x), key(x), index(x))
  outseq
}




