#' @importFrom rmisc Curry
NULL

#' @keywords internal
#' @autoImports
find_neighbors <- function (query, subject, n = 5, 
                            direction = 'flanking',
                            include_key = 'any',
                            exclude_key = 'none') {
  
  direction <- match.arg(direction, c("flanking","downstream","upstream"))
  
  if (!(is(query, "gbFeatureList") || is(query, "gbFeature"))) {
    stop("class ", sQuote(class(query)), " no supported query")
  }
  
  if (!(is(subject, "gbRecord") || is(subject, "gbFeatureList"))) {
    stop("class ", sQuote(class(subject)), " no supported subject")
  }
  
  if (is(query, "gbFeatureList") || is(query, "gbFeature")) {
    query <- ranges(ranges(query))
  }
  
  if (is(subject, "gbRecord")) {
    subject <- features(subject)
  }
  
  if (include_key != "any") {
    subject <- select(subject, key=include_key)
    if (length(subject) == 0L) {
      message("Include key ", sQuote(paste(include_key, collapse=",")), " returns no features")
      return(NULL)
    }
  }
  
  if ("none" %in% exclude_key) {
    # exclude at least "source" from gbFeatureList
    exclude_key <- "source"
  } else {
    exclude_key <- c(exclude_key, "source")
  }
  
  subject <- subject[vapply(subject, function (f) {
    f@key %ni% exclude_key 
  }, logical(1))]
  
  if (length(subject) == 0L) {
    message("Exclude key ", sQuote(paste(exclude_key, collapse=",")), " returns no features")
    return(NULL)
  }
  
  subject_range <- ranges(ranges(subject))
  subject_idx <- vapply(subject, function (f) f@.Id, integer(1))
  
  FUN <- switch(direction,
                upstream=list(IRanges::follow),
                downstream=list(IRanges::precede),
                flanking=list(IRanges::follow, IRanges::precede))
  
  # if "flanking" is called we need to reset the query ranges to
  # the original values
  ori_query <- query
  UID <- vector("list", length(query))
  for (i in seq_along(FUN)) {
    for (j in seq_len(n)) {
      hits <- FUN[[i]](query, subject_range, select="all")
      split_hits <- split(hits, queryHits(hits))
      new_query <- IRanges()
      for (k in seq_along(split_hits)) {
        hit_idx <- subjectHits(split_hits[[k]])
        new_range <- subject_range[hit_idx,]
        UID[[K]] <- c(UID[[k]], subject_idx[hit_idx])
        new_query <- c(new_query, new_range)
      }
      query <- new_query
    }
    query <- ori_query # reset query ranges if a second FUN is called
  }
  
  feature_list <- lapply(UID, function (id) {
    select(subject, idx=unique(sort(id)))
  })
  
  if (length(feature_list) == 1L) {
    return(feature_list[[1L]])
  } else {
    return(feature_list)
  }
}


#' Find flanking features
#' 
#' @usage upstream(query, subject, n=5, include_key="all", exclude_key="none")
#' 
#' @param query The query \code{\linkS4class{gbRange}},
#' \code{\linkS4class{gbFeatureList}}, or
#' \code{\linkS4class{gbFeature}}  instance.
#' @param subject The subject \code{\linkS4class{gbRecord}} or
#' \code{\linkS4class{gbFeatureList}} instance within which the
#' n nearest upstream features are searched.
#' @param n  The number of upstream features to be included
#' @param include_key Which features should be returned. Defaults to
#' \dQuote{all}.
#' @param exclude_key Which feature(s) should be excluded from the search.
#' Defaults to \dQuote{none}.
#' 
#' @rdname upstream
#' @export
upstream <- Curry("find_neighbors", direction = "upstream")


#' @usage downstream(query, subject, n=5, include_key="all", exclude_key="none)
#' @export
#' @rdname upstream
downstream <- Curry("find_neighbors", direction = "downstream")


#' @usage flanking(query, subject, n=5, include_key="all", exclude_key="none)
#' @export
#' @rdname upstream
flanking <- Curry("find_neighbors", direction = "flanking")

