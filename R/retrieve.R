#' @keywords internal
#' @autoImports
.retrieve <- function (x, cols = NULL) {
  
  if (is_empty(x) || is.null(cols))
    return(x)
  
  cols <- gsub("\n|\t", " ", cols)
  cols <- unlist(strsplit(cols, ";"))
  col_names <- character(0)
  i <- k <- l <- q <- NULL
  
  if (any(idx <- grepl("idx|index", cols, ignore.case=TRUE))) {
    i <- vapply(x, function (x) x@.ID, numeric(1))
    col_names <- c(col_names, cols[idx])
    if (is_empty(cols <- cols[!idx]))
      return(.return(i, .Names=col_names))
  }
  
  if (any(idx <- grepl("key", cols, ignore.case=TRUE))) {
    k <- vapply(x, function (x) x@key, character(1))
    col_names <- c(col_names, cols[idx])
    if (is_empty(cols <- cols[!idx]))
      return(.return(i, k, .Names=col_names))
  }

  if (any(idx <- grepl("location|range|start|end|width|strand", cols, ignore.case=TRUE))) {
    
    if (any(grepl("location|range", cols, ignore.case=TRUE))) {
      l <- range(x, join = TRUE)
    } else {
      col_names <- c(col_names, loc <- cols[idx])
      l <- lapply(loc, function (loc) {
        switch(loc, start=start(x), end=end(x), width=width(x), strand=strand(x))
      })
    }

    if (is_empty(cols <- cols[!idx]))
      return( .return(i, k, l, .Names=col_names) )
  }

  col_names <- c(col_names, cols)
  q <- lapply(cols, function (qual) {
    lapply(x, function (x) .qualAccess(x, qual, TRUE)) %@% .simplify
  })
  
  ## special treatment for db_xref
  if ("db_xref" %in% cols) {
    i.dbx <- which(cols %in% "db_xref")
    dbx_list <- parseDbXref(q[[i.dbx]])
    q <- append(q[-i.dbx], dbx_list, i.dbx - 1)
    dbx_name <- which(col_names %in% "db_xref")
    col_names <- append(col_names[-dbx_name], names(dbx_list), dbx_name - 1)
  }
  
  if (length(q) == 1L) q <- q[[1L]]
  
  .return( i, k, l, q, .Names=col_names )
}

# args = list(i, k, l, q)
#' @keywords internal
#' @autoImports
.return <- function (..., .Names) {
  args <- list(...)
  zero <- vapply(args, is.null, logical(1L))
  args[zero] <- NULL
  if (!any(hasList(args)) && length(args) == 1L)
    structure(args[[1L]], names = rep(.Names, length(args[[1L]])))
  else if (!any(hasListOfLists(args))) {
    args <- flatten(args)
    if (any(r <- vapply(args, function (x) is(x, "gbRange"), logical(1)))) {
      range <- args[r][[1]]
      df <-  DataFrame(range@elementMetadata, args[!r])
      dimnames(df)[[2]] <- c("strand", .Names)
      range@elementMetadata <- df
      range
    } else {
      structure(
        data.frame(stringsAsFactors = FALSE, args),
        names=.Names
      )
    }
  } else {
    args <- flatten(args, stop.at = 2)
    if (any(r <- vapply(args, function (x) is(x, "gbRange"), logical(1)))) {
      .Names <- append(.Names, "range", which(r) - 1)
    }
    structure(args, names=.Names)
  }
}

hasList <- function (x) vapply(x, typeof, character(1L)) == "list"

hasListOfLists <- function (x) unlist(lapply(x[hasList(x)], hasList))

matchIdx <- function (x, idx) {
  x_idx <- vapply(x, function (f) f@.ID, numeric(1))
  which(x_idx %in% idx)  
}

#' @autoImports
parseDbXref <- function (dbx) {
  n <- if (is.null(n <- nrow(dbx))) length(dbx) else n
  
  if (is.atomic(dbx)) {
    structure(list(strsplitN(dbx, ":", 2)), names=unique(strsplitN(dbx, ":", 1)))
  } else if (is.data.frame(dbx)) {
    dbs <- vapply(dbx, strsplitN, ":", 1, FUN.VALUE=character(1)) %@%
      unique %@% `dim<-`(NULL)
    structure(lapply(dbx, strsplitN, ":", 2), names = dbs)
  } else if (is.list(dbx)) {
    l <- lapply(dbx, function (x) {
      a <- sapply(x, strsplit, ":")
      a <- setNames(sapply(a, "[", 2), sapply(a, "[", 1))
      data.frame(stringsAsFactors=FALSE, as.list(a))
    }) %@% rbind.fill
    as.list(l)
  }
}

