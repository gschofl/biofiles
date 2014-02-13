#' Select annotations from a FeatureTable
#'
#' @param x A gbFeatureTable
#' @param ... character keys
#' @param .cols gbFeature annotations returned as a data.frame
#'
#' @keywords internal
#' @importFrom IRanges elementMetadata "elementMetadata<-"
#' @include filter.R
.select <- function (x, ..., .cols = NULL) {
  cols <- unique(c(c(...), .cols))
  if (all_empty(x) || is.empty(cols)) {
    return(x)
  }
  column_names <- character(0)
  i <- k <- l <- q <- NULL
  idx <- is_index(cols)
  if (any(idx)) {
    i <- index(x)
    column_names <- c(column_names, cols[idx])
    cols <- cols[!idx]
    if (all_empty(cols)) {
      return(.return(i, .Names=column_names))
    }
  }
  idx <- is_key(cols)
  if (any(idx)) {
    k <- vapply(x, slot, name = 'key', FUN.VALUE = '')
    column_names <- c(column_names, cols[idx])
    cols <- cols[!idx]
    if (all_empty(cols)) {
      return(.return(i, k, .Names=column_names))
    }
  }
  idx <- grepl("start|end|width|strand", cols)
  if (any(idx)) {
    loc <- cols[idx]
    column_names <- c(column_names, loc)
    l <- lapply(loc, function(loc) {
      switch(loc,
             start  = start(x),
             end    = end(x),
             width  = width(x),
             strand = strand(x))
    })
    cols <- cols[!idx]
    if (all_empty(cols)) {
      return(.return(i, k, l, .Names=column_names))
    }
  }
  q <- .simplify(.qual_access(x, which = cols, fixed = TRUE), unlist = FALSE)
  .return(i, k, l, q, .Names=column_names)
}

.return <- function(..., .Names) {
  cols <- flatten1(compact(list(...)))
  if (is.null(names(cols))) {
    names(cols) <- .Names
  } else {
    names(cols)[names(cols) == ''] <- .Names
  }
  data.frame(stringsAsFactors = FALSE, cols)
}

#' @importFrom plyr rbind.fill
parseDbXref <- function(dbx) {
  n <- nrow(dbx)
  n <- if (is.null(n)) length(dbx) else n
  if (is.atomic(dbx)) {
    structure(list(strsplitN(dbx, ":", 2)),
              names=unique(strsplitN(dbx, ":", 1)))
  } else if (is.data.frame(dbx)) {
    dbs <-`dim<-`(
        vapply(dbx, Compose(unique, strsplitN), ":", 1, FUN.VALUE=character(1)),
        NULL
      )
    structure(lapply(dbx, strsplitN, ":", 2), names = dbs)
  } else if (is.list(dbx)) {
    as.list(
      rbind.fill(
        lapply(dbx, function (x) {
          a <- sapply(x, strsplit, ":")
          a <- setNames(sapply(a, "[", 2), sapply(a, "[", 1))
          data.frame(stringsAsFactors=FALSE, as.list(a))
        })
      )
    )
  }
}

