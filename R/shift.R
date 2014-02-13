#' @importFrom Biostrings reverseComplement xscat subseq
NULL

setAs(from="gbFeatureTable", to="gbRecord", function(from) {
  new_gbRecord(seqinfo = from@.seqinfo, features = from, contig = NULL)
})

update_split <- function(x, split_matrix) {
  if (!is.na(x@location@compound)) {
    stop("Cannot split a compound location")
  }
  x@location@range <- split_matrix
  x@location@fuzzy <- matrix(c(FALSE, TRUE, TRUE, FALSE), ncol=2)
  x@location@compound <- "join"
  x@location@accession <- rep(x@location@accession, 2)
  x@location@remote <- rep(x@location@remote, 2)
  x
}

merge_split <- function(fList, omit_unmergables = FALSE) {
  ## unmergables:
  ## genes that span the start/end of a circular chromosome
  ## genes with trans-splicing structure
  ## e.g.: join(complement(331719..331758),487030..487064)
  assert_that(all(is_compound(fList)))
  lapply(fList, function(f) {
    lRange <- f@location@range 
    if (lRange[-nrow(lRange), 2] - lRange[-1,1] < 1 &&
          length(strand <- unique(f@location@strand)) < 2) {
      f@location@range     <- matrix(range(lRange), ncol=2)
      f@location@fuzzy     <- matrix(c(FALSE, FALSE), ncol=2)
      f@location@strand    <- strand
      f@location@compound  <- NA_character_
      f@location@accession <- f@location@accession[1]
      f@location@remote    <- f@location@remote[1]
      f@location@type      <- f@location@type[1]
      f
    } else {
      if (omit_unmergables)
        NULL
      else
        f
    }
  })
}

.shift <- function(x, shift = 0L, split = FALSE, order = TRUE) {
  assert_that(is(x, "gbRecord") || is(x, "gbFeatureTable"))
  was.gbRecord <- FALSE
  if (is(x, "gbRecord")) {
    x <- .features(x)
    was.gbRecord <- TRUE
  }
  len <- getLength(x)
  if (all_empty(x["source"])) {
    stop("No source key in this gbFeatureTable")
  }  
  src <- x["source"]
  if (all_empty(src)) {
    stop("No source key available")
  }
  f <- x[-index(src)]
  
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
        new_start[end_only] <- Map(function(x) x[, 1], sm)
        new_end[end_only] <- Map(function(x) x[, 2], sm)
      } else {
        stop("This shiftwidth would split features ", paste(end_only, collapse=", "))
      }
    }
  } 
  start(f, check=FALSE) <- new_start
  end(f, check=FALSE) <- new_end
  cmpnd <- which(is_compound(f))
  if (length(cmpnd) > 0) {
    f[cmpnd] <- merge_split(f[cmpnd], omit_unmergables=FALSE)
  }
  f <- if (order) f[order(mapply("[", new_start, 1L))] else f
  ## update sequence
  seqinfo <- f@.seqinfo$clone()
  seq <- .sequence(seqinfo)
  seq_len <- seq@ranges@width
  if (shift > 0) {
    seqinfo$sequence <- xscat(subseq(seq, seq_len - shift + 1),
                              subseq(seq, 1, seq_len - shift))
  } else {
    seqinfo$sequence <- xscat(subseq(seq, abs(shift) + 1),
                              subseq(seq, 1, abs(shift)))
  }
  names(seqinfo$sequence) <- getAccession(f)
  x <- new('gbFeatureTable', .Data = c(src, f), .id = c(src@.id, f@.id), .seqinfo = seqinfo)
  x <- if (was.gbRecord) as(x, "gbRecord") else x
  x
}


.revcomp <- function(x, order = TRUE) {
  assert_that(is(x, "gbRecord") || is(x, "gbFeatureTable"))
  was.gbRecord <- FALSE
  if (is(x, "gbRecord")) {
    x <- .features(x)
    was.gbRecord <- TRUE
  }
  cmpnd <- which(is_compound(x))
  if (length(cmpnd) > 0) {
    x[cmpnd] <- merge_split(fList = x[cmpnd], omit_unmergables=TRUE)
    x <- x[vapply(x, length, 0L) > 0L]
  }
  len <- getLength(x)
  new_end <- Map(function(s) len - s + 1, start(x))
  new_start <- Map(function(e) len - e + 1, end(x))
  new_strand <- Map(`-`, strand(x))
  start(x, check=FALSE) <- new_start
  end(x) <- new_end
  strand(x) <- new_strand
  x <- if (order) x[order(unlist(Map(min, new_start)))] else x
  ## update sequence
  seqinfo <- x@.seqinfo$clone()
  seq <- .sequence(seqinfo)
  seqinfo$sequence <- reverseComplement(seq)
  names(seqinfo$sequence) <- getAccession(x)
  x <- new('gbFeatureTable', .Data = x, .id = x@.id, .seqinfo = seqinfo)
  x <- if (was.gbRecord) as(x, "gbRecord") else x
  x
}

