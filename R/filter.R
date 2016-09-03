#' Filter gbFeatures from a FeatureTable
#'
#' @param x A gbFeatureTable
#' @param ... key=value pairs interpreted as filters
#' @param .cols gbFeature annotations returned as a data.frame
#' @rdname dot-filter
#' @keywords internal
.filter <- function(x, ..., .cols = NULL) {
  keys <- parse_keys(...)
  if (is.null(keys)) {
    return(x)
  }
  # restrict features to the selected indices
  if (!is.null(keys$idx)) {
    x <- x[ match_idx(x, idx = keys$idx) ]
  }
  # for the selected indices restrict features to
  # the selected range
  if (!is.null(keys$range)) {
    start <- keys$range$start
    end <- keys$range$end
    subject_range <- .IRanges(x)
    start[is.na(start)] <- 1
    end[is.na(end)] <- max(end(subject_range))
    query_range <- IRanges::IRanges(start, end)
    ovl <- IRanges::findOverlaps(query_range, subject_range)
    x <- x[S4Vectors::subjectHits(ovl)]
  }
  # for the selected range restrict features to
  # the selected keys 
  if (!is.null(keys$key)) {
    if (is(keys$key, 'regexp')) {
      x <- x[grep(pattern = keys$key, x = key(x))]
    } else {
      x <- x[keys$key]
    }  
  }
  # for the selected keys restrict to features with
  # the specified qualifiers
  feature_idx <- is_feature(names(keys))
  if (any(feature_idx)) {
    features <- keys[feature_idx]
    idx <- Map(function(tag, val) {
      if (isTRUE(val)) {
        ## test for the presence of a qualifier tag among the qualifiers
        ## of each feature
        idx <- vapply(x, function(x_item) {
          any(grepl(pattern = tag, x = names(x_item@qualifiers)))
        }, FALSE)
      } else {
        ## get indices for matching tags for each feature
        tag_idx_lst <- 
          lapply(x, function(x_item) grep(pattern = tag, x = names(x_item@qualifiers)))
        idx <- mapply(function(x_item, tag_idx) { 
          any(grepl(pattern = val[[1]], x = x_item@qualifiers[tag_idx]))
        }, x, tag_idx_lst, USE.NAMES=FALSE)
      }   
    }, tag=names(features), val=features)    
    x <- x[ Reduce("&", idx) ]
  }
  # return features
  if (!is.null(.cols)) {
    select(x, .cols = .cols)
  } else {
    x
  }
}

is_range <- function(x) grepl("^range", x)

is_key <- function(x) grepl("^key", x)

is_index <- function(x) grepl("^in?de?x", x)

is_feature <- function(x) !(is_range(x) | is_key(x) | is_index(x))

parse_keys <- function(...) {
  keys <- list(...)
  if (all_empty(keys)) {
    return(NULL)
  }
  key_names <- names(keys)
  names(keys)[is_index(key_names)] <- "idx"
  if (!all(nzchar(key_names))) {
    stop("Provide named keys")
  }
  if (!is.na(charmatch("range", key_names))) {
    keys$range <- parse_range(keys$range)
  }
  if (!is.na(charmatch("key", key_names))) {
    keys$key <- if (length(keys$key) > 1) re(collapse(keys$key, "|")) else keys$key
  }
  feature_idx <- is_feature(key_names)
  if (sum(feature_idx) > 0) {
    for (i in which(feature_idx)) {
      if (!is.logical(keys[[i]])) {
        keys[i] <- re(collapse(keys[[i]], "|"))
      } else {
        keys[i] <- unique(keys[[i]])
      }
    }
  }
  keys
}

parse_range <- function(range) {
  if (all_empty(range)) {
    return(NULL)
  }
  range <- uusplit(range, ",", fixed=TRUE)
  start <- end <- numeric(0)
  range <- strsplit(range, "..", fixed=TRUE)
  range <- list(start = as.numeric(vapply(range, `[`, 1L, FUN.VALUE = "")),
                end   = as.numeric(vapply(range, `[`, 2L, FUN.VALUE = "")))
  range
}

match_idx <- function(x, idx) {
  which(index(x) %in% idx)
}


