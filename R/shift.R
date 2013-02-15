.shift_features <- function (x, shift=0L, split=FALSE, order=FALSE,
                             updateDb=FALSE) {
  
  if (is(x, "gbRecord")) {
    len <- x$length
    features <- x$features
  } else if (is(x, "gbFeatureList")) {
    if (!any(hasKey(x, "source"))) {
      stop("No source key in this gbFeatureList")
    }
    len <- end(x["source"])
    features <- x
  }
  
  update_split <- function(x, split_matrix) {
    if (not.na(x@location@compound)) {
      stop("Cannot split a compound location")
    }
    x@location@.Data <- split_matrix
    x@location@compound <- "join"
    x@location@partial <- matrix(c(FALSE, TRUE, TRUE, FALSE), ncol=2)
    x@location@accession <- rep(x@location@accession, 2)
    x@location@remote <- rep(x@location@remote, 2)
    x@location@closed <- matrix(rep(x@location@closed, 2), ncol=2)
    x
  }
  
  src <- features[1]
  f <- features[-1]
  
  start_pos <- start(f)
  end_pos <- end(f)
  
  new_start <- Map("+", start_pos, shift)
  new_end <- Map("+", end_pos, shift)
  
  exceeds_len_start <- which(mapply(function (x) any(x > len), new_start) |
    mapply(function (x) any(x < 0L), new_start)) 
  exceeds_len_end <-  which(mapply(function (x) any(x > len), new_end) |
    mapply(function (x) any(x < 0L), new_end))
  
  if (length(exceeds_len_start) > 0L || length(exceeds_len_end) > 0L) {
    
    start_end <- intersect(exceeds_len_start, exceeds_len_end)
    
    if (length(start_end) > 0L) {
      get_len <- function (x, len) ifelse(x > len, x - len, ifelse(x < 0L, len + x, x))
      new_start[start_end] <- Map(get_len, new_start[start_end], len)
      new_end[start_end] <- Map(get_len, new_end[start_end], len)
    }

    end_only <- setdiff(exceeds_len_end, exceeds_len_start)
    
    if (length(end_only) > 0L) {
      if (split) {
        ss <- mapply("-", new_start[end_only], len)
        se <- mapply("-", new_end[end_only], len)
        ## Split Matrix
        sm <- Map(function(ss, se) matrix(c(len + ss, 1, len, se), ncol=2), ss, se)
        f[end_only] <- Map(update_split, x=f[end_only], split_matrix=sm)
        new_start[end_only] <- Map(function(x) x[,1], sm)
        new_end[end_only] <- Map(function(x) x[,2], sm)
      } else {
        stop("This shiftwidth would split feature(s) ", paste(end_only, collapse=", "))
      }
    }
    
  }
  
  start(f) <- new_start
  end(f) <- new_end
  
  if (order) {
    f <- f[order(mapply("[", new_start, 1))]
  }
  
  f <- new('gbFeatureList', .Data=c(src, f), .Dir=src@.Dir,
           .ACCN=src@.ACCN, .DEF=src@.DEF)
  
  if (updateDb) {
    db <- init_db(src@.Dir, verbose=FALSE)
    dbInsert(db, key="features", value=f)
    
    seq <- dbFetch(db, "sequence")
    
    if (shift > 0) {
      shift_point <- seq@ranges@width - shift + 1L
    } else {
      shift_point <- 0L - shift + 1L
    }
    
    new_seq <- xscat(subseq(seq, start = shift_point), 
                     subseq(seq, start = 1L, end = shift_point - 1L))
    names(new_seq) <- names(seq)
    dbInsert(db, key="sequence", value=new_seq)
  }
  
  return( f )
}

