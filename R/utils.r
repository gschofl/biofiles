#' @importClassesFrom intervals Intervals_full
#' @importClassesFrom intervals Intervals_virtual
#' @importClassesFrom intervals Intervals_virtual_or_numeric
#' @importClassesFrom IRanges IRanges
#' @importClassesFrom IRanges Ranges
#' @importClassesFrom IRanges IntegerList
#' @importClassesFrom IRanges RangesORmissing
#' @importClassesFrom IRanges List
#' @importClassesFrom IRanges AtomicList
#' @importClassesFrom IRanges Vector
#' @importClassesFrom IRanges Annotated
#' @importClassesFrom IRanges DataFrame
#' @importClassesFrom filehash filehashRDS
#' @importClassesFrom filehash filehash
#'  
#' @import rmisc
#' 
#' @importFrom intervals closed
#' 
#' @importFrom stats start
#' @importFrom stats end
#' 
#' @importFrom IRanges DataFrame
#' @importFrom IRanges IRanges
#' @importFrom IRanges precede
#' @importFrom IRanges follow
#' @importFrom IRanges "start<-"
#' @importFrom IRanges "end<-"
#' @importFrom IRanges width
#' @importFrom IRanges shift
#' @importFrom IRanges IntervalTree
#' @importFrom IRanges findOverlaps
#' @importFrom IRanges queryHits
#' @importFrom IRanges subjectHits
#' @importFrom IRanges as.data.frame
#' @importFrom IRanges elementMetadata
#' @importFrom IRanges "elementMetadata<-"
#' @importMethodsFrom IRanges coerce
#' 
#' @importFrom filehash dbCreate
#' @importFrom filehash dbInit
#' @importFrom filehash dbFetch
#' @importFrom filehash dbInsert
#' @importFrom filehash dbDelete
#' 
#' @importFrom Biostrings read.DNAStringSet
#' @importFrom Biostrings read.RNAStringSet
#' @importFrom Biostrings read.AAStringSet
#' @importFrom Biostrings write.XStringSet
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings RNAStringSet
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings xscat
#' @importFrom Biostrings subseq
#' @importFrom Biostrings toString
#' 
#' @importFrom stringr str_extract
#' @importFrom stringr str_detect
#' 
#' @importFrom plyr rbind.fill
#' 
#' @importFrom parallel mcmapply
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
NULL


recycle <- function (x, val) {
  lx <- length(x)
  lv <- length(val)
  if (lx > lv) {
    val <- c(rep(val, lx%/%lv), val[seq_len(lx%%lv)])
  }
  val
}


"%||%" <- function (a, b) {
  if (!a) b else a
}


merge_lines <- function (lines) {
  if (length(lines) == 1L) {
    trim(lines)
  } else {
    paste0(trim(lines), collapse=" ")
  }
}


getCompounds <- function (x) {
  x <- x[which(is.compound(x))]
  if (length(x) == 0) return(NA_real_) 
  cL <- vapply(x, function (f) nrow(f@location@.Data), numeric(1))
  cL
}


expandIds <- function (x) {
  cmp_pos <- Position(is.compound, x)
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


# Access qualifiers from gbFeature or gbFeatureList objects
.qualAccess <- function (x, qual = "", fixed = FALSE) {
  
  .access <- function (q) {
    q <- q@qualifiers
    if (fixed) qual <- wrap(qual, "\\b") 
    n <- length(q)
    
    if (n == 0) {
      return( structure(rep(NA_character_, length(qual)),
                        names=trim(qual, "\\\\b")) )
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


.seqAccess <- function (s, x, type) {
  
  if (is.null(s))
    stop("No sequence available")
  
  # merge Sequences
  merge_seq <- function (s, x, type) {
    if (length(start(x)) == 1L) {
      seq <- subseq(s, start=start(x), end=end(x))
    } else {
      seq <- do.call(xscat, Map(subseq, s, start=start(x), end=end(x)))
    }
    seq <- switch(type,
                  DNA=DNAStringSet(seq),
                  AA=AAStringSet(seq),
                  RNA=RNAStringSet(seq))
    seq@ranges@NAMES <- sprintf("%s.%s.%s", x@.ACCN, x@key, x@.ID)
    seq
  }
  
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

