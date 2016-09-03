#' @keywords internal
find_neighbors <- function(query, subject, n = 5,  direction = 'flanking',
                           include_key = 'any', exclude_key = 'none') {
  
  direction <- match.arg(direction, c("flanking","downstream","upstream"))
  if (!(is(query, "gbFeatureTable") || is(query, "gbFeature"))) {
    stop("class ", sQuote(class(query)), " no supported query")
  }
  if (!(is(subject, "gbRecord") || is(subject, "gbFeatureTable"))) {
    stop("class ", sQuote(class(subject)), " no supported subject")
  }
  if (is(query, "gbFeatureTable") || is(query, "gbFeature")) {
    query <- .IRanges(query)
  }
  if (is(subject, "gbRecord")) {
    subject <- .features(subject)
  }
  if (include_key != "any") {
    subject <- filter(subject, key = include_key)
    if (length(subject) == 0L) {
      message("Include key ", sQuote(paste0(include_key, collapse = ",")), " returns no features")
      return(NULL)
    }
  }
  if ("none" %in% exclude_key) {
    # exclude at least "source" from gbFeatureTable
    exclude_key <- "source"
  } else {
    exclude_key <- c(exclude_key, "source")
  }
  subject <- subject[key(subject) %ni% exclude_key]
  if (length(subject) == 0L) {
    message("Exclude key ", sQuote(paste0(exclude_key, collapse = ",")), " returns no features")
    return(NULL)
  }
  subject_range <- .IRanges(subject)
  subject_idx <- index(subject)
  FUN <- switch(direction,
                upstream   = list(IRanges::follow),
                downstream = list(IRanges::precede),
                flanking   = list(IRanges::follow, IRanges::precede)) 
  # if "flanking" is called we need to reset the query ranges to
  # the original values
  ori_query <- query
  UID <- vector("list", length(query))
  for (f in FUN) {
    for (i in seq_len(n)) {
      ## find the nearest range neighbor following|preceding the query range(s)
      hits <- f(x = query, subject = subject_range, select = "all")
      ## with multiple queries we have to split the hits and loop over them
      ## with a single query this takes no effect anyways.
      split_hits <- S4Vectors::split(hits, S4Vectors::queryHits(hits))
      new_query <- IRanges::IRanges()
      ## loop over potentially multiple queries
      for (j in seq_along(split_hits)) {
        hit_idx <- S4Vectors::subjectHits(split_hits[[j]])
        new_range <- subject_range[hit_idx, ]
        UID[[j]] <- c(UID[[j]], subject_idx[hit_idx])
        new_query <- c(new_query, new_range)
      }
      ## set queries to be the hits from the previous iteration
      query <- new_query
    }
    ## reset query ranges if a second FUN is called
    query <- ori_query 
  }
  ## now we should have a list of integer vectors (feature indices), the lengh
  ## of which depends on the number of query features.
  f_list <- lapply(UID, function(id) {
    filter(subject, idx = unique(sort(id)))
  })
  if (length(f_list) == 1L) {
    f_list[[1L]]
  } else {
    f_list
  }
}

#' Find flanking features.
#' 
#' @param query A \code{\linkS4class{gbFeature}} or \code{\linkS4class{gbFeatureTable}}
#' object.
#' @param subject A \code{\linkS4class{gbRecord}} or \code{\linkS4class{gbFeatureTable}}
#' object within which the n nearest upstream features are found.
#' @param n  The number of upstream features to be returned.
#' @param include_key Which features should be returned. Defaults to
#' \dQuote{all}. 
#' @param exclude_key Which feature(s) should be excluded from the search.
#' Defaults to \dQuote{none}.
#' @return A (list of) \code{\linkS4class{gbFeatureTable}}s.
#' @rdname upstream
#' @export
#' @examples
#' load(system.file("extdata", "S_cerevisiae_mito.rda", package = "biofiles"))
#' cytb <- ft(filter(x, product = "^cytochrome b$"))
#' 
#' ## find the three nearest upstream neighbor CDS to CYTB 
#' upstream(cytb, x["CDS"], n = 3)
upstream <- function(query, subject, n = 5, include_key = "all", exclude_key = "none") {
 find_neighbors(query = query, subject = subject, n = n, include_key = include_key,
                exclude_key = exclude_key, direction = "upstream")
}

#' @rdname upstream
#' @export
downstream <- function(query, subject, n = 5, include_key = "all", exclude_key = "none") {
  find_neighbors(query = query, subject = subject, n = n, include_key = include_key,
                 exclude_key = exclude_key, direction = "downstream")
}

#' @rdname upstream
#' @export
flanking <- function(query, subject, n = 5, include_key = "all", exclude_key = "none") {
  find_neighbors(query = query, subject = subject, n = n, include_key = include_key,
                 exclude_key = exclude_key, direction = "flanking")
}

