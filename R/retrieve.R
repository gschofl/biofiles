#' @keywords internal
#' @importFrom IRanges elementMetadata "elementMetadata<-"
.retrieve <- function (x, cols = NULL) {
  
  if (all_empty(x) || is.null(cols))
    return(x)
  
  cols <- gsub("\n|\t", " ", cols)
  cols <- unlist(strsplit(cols, ";"))
  col_names <- character(0)
  i <- k <- l <- q <- NULL
  
  idx <- grepl("idx|index", cols, ignore.case=TRUE)
  if (any(idx)) {
    i <- vapply(x, function (x) x@.id, numeric(1))
    col_names <- c(col_names, cols[idx])
    cols <- cols[!idx]
    if (all_empty(cols))
      return( .return(i, .Names=col_names) )
  }
  
  idx <- grepl("key", cols, ignore.case=TRUE)
  if (any(idx)) {
    k <- vapply(x, function (x) x@key, character(1))
    col_names <- c(col_names, cols[idx])
    cols <- cols[!idx]
    if (all_empty(cols))
      return( .return(i, k, .Names=col_names) )
  }

  idx <- grepl("location|range|start|end|width|strand", cols, ignore.case=TRUE)
  if (any(idx)) {
    if (any(grepl("location|range", cols, ignore.case=TRUE))) {
      l <- ranges(x, join = TRUE)
    } else {
      col_names <- c(col_names, loc <- cols[idx])
      l <- lapply(loc, function (loc) {
        switch(loc,
               start=start(x),
               end=end(x),
               width=width(x),
               strand=strand(x))
      })
    }

    if (all_empty(cols <- cols[!idx]))
      return( .return(i, k, l, .Names=col_names) )
  }
  
  col_names <- c(col_names, cols)
  q <- .simplify( .qual_access(x, cols, TRUE) )
  
  ## special treatment for db_xref
  if ("db_xref" %in% cols) {
    if (length(cols) == 1) {
      dbx_list <- parseDbXref(dbx=q)
      q <- dbx_list
    } else {
      dbx_list <- parseDbXref(q[[which(cols %in% "db_xref")]])
      q <- append(q[-i.dbx], dbx_list, i.dbx - 1)
    }
    dbx_name <- which(col_names %in% "db_xref")
    col_names <- append(col_names[-dbx_name], names(dbx_list), dbx_name - 1)
  }
  
  if (length(q) == 1L) q <- q[[1L]]
  
  .return( i, k, l, q, .Names=col_names )
}

# args = compact(list(i, k, l, q))
#' @keywords internal
.return <- function (..., .Names) {
  args <- compact(list(...))
  if (!any(hasList(args)) && length(args) == 1L) {
    structure(
      args[[1L]],
      names = rep(.Names, length(args[[1L]]))
    )
  } else if (!any(hasListOfLists(args))) {
    args <- rmisc::flatten(args)
    r <- vapply(args, is, "GRanges", FUN.VALUE=logical(1))
    if ( any(r) ) {
      range <- args[r][[1]]
      cols <- DataFrame(args[!r])
      names(cols) <- .Names
      elmd <- IRanges::cbind(elementMetadata(range), cols)
      elementMetadata(range) <- elmd
      range
    } else {
      structure(
        data.frame(stringsAsFactors = FALSE, args),
        names=.Names
      )
    }
  } else {
    args <- rmisc::flatten(args, stop.at = 2)
    r <- vapply(args, is, "GRanges", FUN.VALUE=logical(1))
    if ( any(r) ) {
      .Names <- base::append(.Names, "range", base::which(r) - 1)
    }
    structure(args, names=.Names)
  }
}


hasList <- function (x) vapply(x, typeof, character(1)) == "list"


hasListOfLists <- function (x) unlist(lapply(x[hasList(x)], hasList))


matchIdx <- function (x, idx) {
  x_idx <- vapply(x, function (x) x@.id, numeric(1))
  which(x_idx %in% idx)  
}


#' @importFrom plyr rbind.fill
parseDbXref <- function (dbx) {
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

