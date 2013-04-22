.shift_features <- function (x, shift=0L, split=FALSE, order=FALSE) {

  if (is(x, "gbRecord")) {
    len <- seqlengths(x)
    features <- features(x)
  } else if (is(x, "gbFeatureList")) {
    if (all_empty(x["source"])) {
      stop("No source key in this gbFeatureList")
    }
    len <- end(x["source"])
    features <- x
  }
  
  update_split <- function(x, split_matrix) {
    if (not.na(x@location@compound)) {
      stop("Cannot split a compound location")
    }
    x@location@range <- split_matrix
    x@location@compound <- "join"
    x@location@fuzzy <- matrix(c(FALSE, TRUE, TRUE, FALSE), ncol=2)
    x@location@accession <- rep(x@location@accession, 2)
    x@location@remote <- rep(x@location@remote, 2)
    x@location@closed <- matrix(rep(x@location@closed, 2), ncol=2)
    x
  }
  
  merge_split <- function(x) {
    lapply(x, function(x) {
      ld <- location(x)@.Data
      if (ld[-nrow(ld),2] - ld[-1,1] < 1) {
        x@location@.Data <- matrix(range(ld), ncol=2)
        x@location@compound <- NA_character_
        x@location@partial <- matrix(c(FALSE, FALSE), ncol=2)
        x@location@accession <- x@location@accession[1]
        x@location@remote <- x@location@remote[1]
        x@location@closed <- x@location@closed[1,,drop=FALSE]
      }
      x
    })
  }
  
  src <- features["source"]
  if (all_empty(src)) {
    stop("No source key available")
  }
  f <- features[-index(src)]
  
  start_pos <- start(f)
  end_pos <- end(f)
  new_start <- Map("+", start_pos, shift)
  new_end <- Map("+", end_pos, shift)
  
  exceeds <- function(x) any(x > len) | any(x <= 0)
  exceeds_len_start <- which(mapply(exceeds, new_start)) 
  exceeds_len_end <-  which(mapply(exceeds, new_end))
  
  if (length(exceeds_len_start) > 0L || length(exceeds_len_end) > 0L) {
    start_end <- intersect(exceeds_len_start, exceeds_len_end)
    if (length(start_end) > 0L) {
      get_len <- function (x, len) ifelse(x > len, x - len, ifelse(x <= 0L, len + x, x))
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
        stop("This shiftwidth would split features ", paste(end_only, collapse=", "))
      }
    }
  }
  
  start(f, check=FALSE) <- new_start
  end(f, check=FALSE) <- new_end

  cmpnd <- which(is_compound(f))
  if (length(cmpnd) > 0) {
    f[cmpnd] <- merge_split(f[cmpnd])  
  }
  
  if (order) {
    f <- f[order(mapply("[", new_start, 1))]
  }
  
  new('gbFeatureList', .Data=c(src, f), .Info=seqinfo(f))
}

