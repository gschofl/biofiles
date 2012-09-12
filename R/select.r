#' @keywords internal
.select <- function (x, ..., keys = "") {
  
  args <- parseArgs(..., keys=keys)
  
  # restrict features to the selected indices
  if (!isEmpty(args$idx)) {
    x <- x[matchIdx(x, args$idx)]
  }
  
  # for the selected indices restrict features to
  # the selected range
  if (!isEmpty(args$range)) {
    start <- args$range$start
    end <- args$range$end
    subject_range <- IntervalTree(range(x, join=TRUE))
    start[is.na(start)] <- 1
    end[is.na(end)] <- max(end(subject_range))
    query_range <- IntervalTree(IRanges(start, end))
    overlaps <- findOverlaps(query_range, subject_range)
    x <- x[overlaps@subjectHits]
  }
  
  # for the selected range restrict features to
  # the selected keys 
  if (!isEmpty(args$key)) {
    key_idx <- which(grepl(args$key, getKey(x, attributes=FALSE)))
    x <- x[key_idx]
  }
  
  # for the selected keys restrict to features with
  # the specified qualifiers 
  if (any(feature_idx <- isFeature(names(args)))) {
    features <- args[feature_idx]
    idx <- Map(function (tag, val) {
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
        idx <- mapply( function(x_item, tag_idx) { 
          any(grepl(val[[1]], x_item@qualifiers[tag_idx]))
        }, x, tag_idx_lst, USE.NAMES=FALSE)
      }   
    }, tag=names(features), val=features)    
    x <- x[ Reduce("&", idx) ]
  }
  
  # return features
  x  
}


#' @keywords internal
parseArgs <- function (..., keys="") {
  
  args <- list(...)
  
  # rename 'index' to 'idx'
  names(args)[charmatch("index", names(args))] <- "idx"
  
  # evaluate ranges
  range_idx <- charmatch("range", names(args))
  if (!is.na(range_idx)) {
    args[range_idx] <- parseRange(paste0("range=", args[[range_idx]]))
  }
  
  # keys
  key_idx <- isKey(names(args))
  if (sum(key_idx) > 0) {
    if (sum(key_idx) > 1) {
      stop("")
    }
    args[key_idx] <- paste(args[[which(key_idx)]], collapse="|")
  }

  # features
  feature_idx <- isFeature(names(args))
  if (sum(feature_idx) > 0) {
    for (idx in which(feature_idx)) {
      args[idx] <- paste(args[[idx]], collapse="|")
    }
  }
  
  merge.list(args, parseKeys(keys))
}

#' @keywords internal
parseKeys <- function (keys) {
  x <- strsplit(keys, ";")[[1]]
  c(parseIndex(x[isIndex(x)]),
    parseRange(x[isRange(x)]),
    parseKey(x[isKey(x)]),
    parseFeatures(feature=x[isFeature(x)]))
}

#' @keywords internal
parseIndex <- function (index) {
  if (length(index) == 0) {
    list(idx = numeric(0))
  } else {
    idx <- vapply(strsplit(index, "="), "[", 2, FUN.VALUE=character(1))
    idx <- unlist(strsplit(idx, ","))
    idx <- unlist(lapply(idx, function (i) {
      eval(parse(text=i)) 
    }))
    list(idx = unique(idx))
  }
}

#' @keywords internal
parseRange <- function (range) {
  if (length(range) == 0) {
    list(range = numeric(0))
  } else {
    start <- end <- numeric(0)
    r <- vapply(strsplit(range, "="), "[", 2, FUN.VALUE=character(1))
    r <- strsplit(unlist(strsplit(r, ",")), ":")
    r <- list(start = as.numeric(vapply(r, "[", 1, FUN.VALUE=character(1))),
              end = as.numeric(vapply(r, "[", 2, FUN.VALUE=character(1))))
    list(range = r)
  }
}

#' @keywords internal
parseKey <- function (key) {
  if (length(key) == 0) {
    list(key = character(0))
  } else {
    key <- vapply(strsplit(key, "="), "[", 2, FUN.VALUE=character(1))
    key <- paste0(gsub(",", "|", key), collapse="|")
    list(key = key)
  }
}

#' @keywords internal
parseFeatures <- function (feature) {
  if (length(feature) == 0) {
    NULL
  } else {
    tag_val <- strsplit(feature, "=")
    tags <- vapply(tag_val, "[", 1, FUN.VALUE=character(1))
    vals <- sub(",", "|", vapply(tag_val, "[", 2, FUN.VALUE=character(1)))
    vals <- lapply(vals, function (v) if (is.na(v)) TRUE else v )
    setNames(vals, nm=tags)
  }
}

#' @keywords internal
merge.list <- function (x, y, ...) {
  if (length(x) == 0) 
    return(y)
  if (length(y) == 0) 
    return(x)
  i = match(names(y), names(x))
  i = is.na(i)
  if (any(i)) 
    x[names(y)[which(i)]] = y[which(i)]
  x
}

isRange <- function (x) grepl("^range", x)

isKey <- function (x) grepl("^key", x)

isIndex <- function (x) grepl("^idx|^index", x)

isFeature <- function (x) !(isRange(x) | isKey(x) | isIndex(x))

isEmpty <- function (x) length(x) == 0L || !nzchar(x)


#' @keywords internal
.retrieve <- function (x, cols = "") {
  
  if (!nzchar(cols) || length(x) == 0)
    return(x)
  
  cols <- gsub("\n|\t", " ", cols)
  cols <- strsplit(cols, ";")[[1]]
  col_names <- character(0)
  i <- k <- l <- q <- NULL
  
  if (any(idx <- grepl("idx|index", cols, ignore.case=TRUE))) {
    i <- vapply(x, function (x) x@.ID, numeric(1))
    col_names <- c(col_names, cols[idx])
    if (isEmpty(cols <- cols[!idx]))
      return(.return(i, .Names=col_names))
  }
  
  if (any(idx <- grepl("key", cols, ignore.case=TRUE))) {
    k <- vapply(x, function (x) x@key, character(1))
    col_names <- c(col_names, cols[idx])
    if (isEmpty(cols <- cols[!idx]))
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

    if (isEmpty(cols <- cols[!idx]))
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
  }
  else {
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

parseDbXref <- function (dbx) {
  n <- if (is.null(n <- nrow(dbx))) length(dbx) else n
  
  if (is.atomic(dbx)) {
    dbs <-  unique(sapply(strsplit(dbx, ":"), "[", 1))
    structure(list(unlist(lapply(strsplit(dbx, ":"), "[", 2), use.names=FALSE)),
              names = dbs)
  } else if (is.data.frame(dbx)) {
    dbs <- sapply(dbx, function (x) sapply(strsplit(x, ":"), "[", 1)) %@%
      unique %@% `dim<-`(NULL)
    structure(lapply(dbx, function (x) sapply(strsplit(x, ":"), "[", 2)),
              names = dbs)
  } else if (is.list(dbx)) {
    l <- lapply(dbx, function (x) {
      a <- sapply(x, strsplit, split=":")
      a <- setNames(sapply(a, "[", 2), sapply(a, "[", 1))
      data.frame(stringsAsFactors=FALSE, as.list(a))
    }) %@% rbind.fill
    as.list(l)
  }
}

# --R-- vim:ft=r:sw=2:sts=2:ts=4:tw=76:
#       vim:fdm=marker:fmr={{{,}}}:fdl=0