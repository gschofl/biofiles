#' @keywords internal
#' @importFrom IRanges IRanges ranges IntervalTree findOverlaps subjectHits
.select <- function (x, ..., keys = NULL) {
  args <- parse_args(..., keys = keys)
  if (is.null(args)) {
    return(x)
  }
  # restrict features to the selected indices
  if (!all_empty(args[["idx"]])) {
    x <- x[matchIdx(x, args[["idx"]])]
  }
  # for the selected indices restrict features to
  # the selected range
  if (!all_empty(args[["range"]])) {
    start <- args[["range"]][["start"]]
    end <- args[["range"]][["end"]]
    subject_range <- IntervalTree(ranges(ranges(x, join=TRUE)))
    start[is.na(start)] <- 1
    end[is.na(end)] <- max(end(subject_range))
    query_range <- IntervalTree(IRanges(start, end))
    ovl <- findOverlaps(query_range, subject_range)
    x <- x[subjectHits(ovl)]
  }
  
  # for the selected range restrict features to
  # the selected keys 
  if (!all_empty(args$key)) {
    key_idx <- which(grepl(args[["key"]], key(x)))
    x <- x[key_idx]
  }
  
  # for the selected keys restrict to features with
  # the specified qualifiers
  feature_idx <- is_feature(names(args))
  if (any(feature_idx)) {
    features <- args[feature_idx]
    idx <- Map(function(tag, val) {
      if (isTRUE(val)) {
        ## test for the presence of a qualifier tag among the qualifiers
        ## of each feature
        idx <- vapply(x, function(x_item) {
          any(grepl(tag, names(x_item@qualifiers)))
        }, logical(1L))
      } else {
        ## get indices for matching tags for each feature
        tag_idx_lst <- 
          lapply(x, function(x_item) grep(tag, names(x_item@qualifiers)))
        idx <- mapply(function(x_item, tag_idx) { 
          any(grepl(val[[1]], x_item@qualifiers[tag_idx]))
        }, x, tag_idx_lst, USE.NAMES=FALSE)
      }   
    }, tag=names(features), val=features)    
    x <- x[ base::Reduce("&", idx) ]
  }
  
  # return features
  x  
}

is_range <- function (x) grepl("^range", x)

is_key <- function (x) grepl("^key", x)

is_index <- function (x) grepl("^idx|^index", x)

is_feature <- function (x) !(is_range(x) | is_key(x) | is_index(x))

#' @keywords internal
parse_args <- function (..., keys) {
  args <- list(...)
  args <- merge_dups(c(parse_tags(args), parse_tags(keys)))
  
  if (is.null(args))
    return(args)
  
  arg_names <- names(args)

  # evaluate ranges
  range_idx <- charmatch("range", arg_names)
  if (!is.na(range_idx)) {
    args[range_idx] <- parse_range(args[[range_idx]])
  }
  
  # keys
  key_idx <- charmatch("key", arg_names)
  if (!is.na(key_idx)) {
    args[key_idx] <- paste(args[[key_idx]], collapse="|")
  }

  # features
  feature_idx <- is_feature(arg_names)
  if (sum(feature_idx) > 0) {
    for (idx in which(feature_idx)) {
      if (!is.logical(args[[idx]])) {
        args[idx] <- paste0(args[[idx]], collapse="|")
      } else {
        args[idx] <- unique(args[[idx]])
      }
    }
  }
  
  args
}


#' @keywords internal
parse_tags <- function (keys) {
  if (all_empty(keys)) {
    return(NULL)
  }
  x <- if (is.list(keys)) {
    keys <- compact(keys)
    if (all_empty(names(keys))) {
      stop("Require named arguments")
    }
    l_idx <- vapply(keys, Negate(is.logical), logical(1))
    c(
      unlist(
        Map(paste0, names(keys[l_idx]), "=", keys[l_idx], USE.NAMES = FALSE)
      ),
      unlist(
        Map(paste0, names(keys[!l_idx]), USE.NAMES = FALSE)
      )
    )
  } else {
    usplit(keys, ";", fixed=TRUE)
  }
  
  c(
    parse_index(index=x[is_index(x)]),
    prepare_range(range=x[is_range(x)]),
    prepare_key(key=x[is_key(x)]),
    prepare_features(feature=x[is_feature(x)])
  )
}


#' @keywords internal
parse_index <- function (index) {
  if (all_empty(index)) {
    return(NULL)
  }
  idx <- vapply(strsplit(index, "=", fixed=TRUE), "[", 2, FUN.VALUE="")
  idx <- usplit(idx, ",")
  idx <- unlist(lapply(idx, function(i) eval(parse(text=i))))
  list(idx = unique(idx))
}


#' @keywords internal
prepare_range <- function(range) {
  if (all_empty(range)) {
    return(NULL)
  }
  range <- vapply(strsplit(range, "=", fixed=TRUE), "[", 2, FUN.VALUE="")
  list(range = uusplit(range, ",", fixed=TRUE))
}


#' @keywords internal
prepare_key <- function(key) {
  if (all_empty(key)) {
    return(NULL)
  }
  key <- vapply(strsplit(key, "=", fixed=TRUE), "[", 2, FUN.VALUE="")
  list(key = uusplit(key, ",", fixed=TRUE))
}


#' @keywords internal
prepare_features <- function(feature) {
  if (all_empty(feature)) {
    return(NULL)
  }
  tag_val <- strsplit(feature, "=")
  tags <- vapply(tag_val, "[", 1, FUN.VALUE="")
  vals <- vapply(tag_val, "[", 2, FUN.VALUE="")
  vals <- strsplit(vals, ",", fixed=TRUE)
  vals <- merge_dups(as.list(setNames(vals, nm=tags)))
  lapply(vals, function(v) if (all(is.na(v))) TRUE else v )
}


#' @keywords internal
parse_range <- function (range) {
  if (all_empty(range)) {
    return(NULL)
  }
  start <- end <- numeric(0)
  r <- strsplit(range, "..", fixed=TRUE)
  r <- list(start = as.numeric(vapply(r, "[", 1, FUN.VALUE="")),
            end   = as.numeric(vapply(r, "[", 2, FUN.VALUE="")))
  list(range = r)
}

