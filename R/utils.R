recycle <- function (val, len) {
  lv <- length(val)
  if (len > lv) {
    val <- c(rep(val, len%/%lv), val[seq_len(len%%lv)])
  }
  val
}


#' @autoImports
merge_lines <- function (lines) {
  if (length(lines) == 1L) {
    trim(lines)
  } else {
    paste0(trim(lines), collapse=" ")
  }
}


#' @autoImports
is_compound <- function (x) {
  if (is(x, "gbFeatureList")) {
    return(vapply(x, function (f) not.na(f@location@compound), logical(1)))
  } else if (is(x, "gbFeature")) {
    return(not.na(x@location@compound))
  } else if (is(x, "gbLocation")) {
    return(not.na(x@compound))
  }
}


#' @autoImports
getCompounds <- function (x) {
  x <- x[which(is_compound(x))]
  if (length(x) == 0) return(NA_real_) 
  cL <- vapply(x, function (f) nrow(f@location@.Data), numeric(1))
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
    .gbFeatureList()
  }
  id1 <- list(ids=index(xHead, FALSE), keys=key(xHead, FALSE))
  nC <- getCompounds(xCmp)
  id2 <- list(ids=rep(index(xCmp, FALSE), nC),
              keys=rep(key(xCmp, FALSE), nC))
  id3 <- list(ids=index(xTail, FALSE), keys=key(xTail, FALSE))
  Map(function(a, b, c) c(a, b, c), a=id1, b=id2, c=Recall(xTail))
}


#' @autoImports
.qualAccess <- function (x, qual = "", fixed = FALSE) {
  
  .access <- function (q) {
    q <- q@qualifiers
    if (fixed) qual <- wrap(qual, "\\b") 
    n <- length(q)
    
    if (n == 0) {
      return(structure(rep(NA_character_, length(qual)),
                       names=trim(qual, "\\\\b")))
    }
    
    idx <- matrix(
      vapply(qual, grepl, names(q),
             USE.NAMES=FALSE, FUN.VALUE=logical(n)),
      nrow=n)
    n_col <- dim(idx)[2]
    if (n_col == 1L) {
      if (any(idx))
        q[idx]
      else
        structure(NA_character_, names=trim(qual, "\\\\b"))
    } else {
      ans <- lapply(seq.int(n_col), function (i) q[ idx[ ,i]])
      if (any(na_idx <- !apply(idx, 2, any))) {
        for (na in which(na_idx)) {
          ans[[na]] <- structure(NA_character_, names=trim(qual, "\\\\b")[na])
        }
      }
      unlist(ans)
    }
  }
  
  if (is(x, "gbFeature")) {
    .access(x)
  } else if (is(x, "gbFeatureList")) {
    lapply(x, .access)
  }
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


#' @autoImports
.seqAccess <- function (s, x, type) {
  
  if (is.null(s))
    stop("No sequence available")
  
  # merge Sequences
  merge_seq <- function (s, x, type) {
    if (length(start(x)) == 1L) {
      seq <- subseq(s, start=start(x), end=end(x))
    } else {
      seq <- do.call(xscat, Map(Biostrings::subseq, s, 
                                start=start(x), end=end(x)))
    }
    seq <- switch(type,
                  DNA=DNAStringSet(seq),
                  AA=AAStringSet(seq),
                  RNA=RNAStringSet(seq))
    seq@ranges@NAMES <- sprintf("%s.%s.%s", x@.ACCN, x@key, x@.ID)
    seq
  }
  
  #' @autoImports
  if (is(x, "gbFeatureList")) {
    ## initiate empty XStringSet
    seq <- switch(type, 
                  DNA=DNAStringSet(),
                  AA=AAStringSet(),
                  RNA=RNAStringSet())
    for (i in seq_along(x)) {
      seq[i] <- merge_seq(s, x[[i]], type)             
    }
  } else if (is(x, "gbFeature")) {
    seq <- merge_seq(s, x, type)
  }
  
  seq@metadata <- list(definition=x@.DEF, database=x@.Dir)
  seq
}

