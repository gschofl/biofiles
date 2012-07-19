# which="idx=1,2,3,4,8:10,12:14;key=CDS"
# which="idx=1"
# which="idx=1:4,6"
# which="index=2;index=3;idx=40:42"
# which="loc=:20000;key=CDS;product=replication;gene=tnpR;pseudo"
# which="loc=10000:20000,30000:40000,60000:"
# which="location=:10000;key=CDS,gene"
# which="key=CDS,gene"
# which="product"
# which="key=gene;pseudo"
# keys="idx=10:20"
# index=c(1,2,3,8:12)
# subset=571:573

## select features by index, key, location, qualifier values from
## gbData or gbFeatureList objects
.select <- function (x, index = NULL, keys = "") {
  
  indices <- as.numeric(index)
  keys <- gsub("\n|\t", " ", keys)
  keys <- strsplit(keys, ";")[[1]]
  
  #### restrict to selected indices
  if (any(i <- isIndex(keys))) {
    for (term in keys[i]) {
      idx <- strsplit(strsplit(term, "=")[[1]][2], ",")[[1]]
      if (any(r <- grepl(":", idx))) {
        ranges <- strsplit(idx[r], ":")
        low <- sapply(ranges, function(x) as.numeric(x[1]))
        high <- sapply(ranges, function(x) as.numeric(x[2]))
        idx <- c(as.numeric(idx[!r]), unlist(Map(seq, low, high)))
      } else {
        idx <- as.numeric(idx)
      }
      indices <- c(indices, idx)
    }
    indices <- unique(indices)
  }
  if (!isEmpty(indices))
    x <- x[matchIdx(x, indices)]
  
  
  #### restrict to selected location
  if (any(l <- isLocation(keys))) {
    start <- end <- numeric(0)
    for (term in keys[l]) {
      loc <- strsplit(strsplit(strsplit(term, "=")[[1]][2], ",")[[1]], ":")
      start <- c(start, as.numeric(unlist(lapply(loc, "[", 1))))
      end <- c(end, as.numeric(unlist(lapply(loc, "[", 2))))
    }
    
    subject_range <- IntervalTree(range(x, join=TRUE))
    start[is.na(start)] <- min(start(subject_range))
    end[is.na(end)] <- max(end(subject_range))
    query_range <- IntervalTree(IRanges(start, end))
    overlaps <- findOverlaps(query_range, subject_range)
    x <- x[overlaps@subjectHits]
  }
  
  #### within the selected window restrict by key(s) 
  if (any(k <- isKey(keys))) {
    key <- character(0)
    for (term in keys[k]) {
      key <- c(key, gsub(",", "|", strsplit(term, "=")[[1]][2]))
    }
    key <- paste0(key, collapse="|")
    key_idx <- which(grepl(key, getKey(x, attributes=FALSE)))
    x <- x[key_idx]
  }
  
  #### within resticted window and for restricted key,
  #### search pattern for specified qualifiers 
  if (any(q <- !( isIndex(keys) |  isLocation(keys) | isKey(keys) ))) {
    tag_value <- strsplit(keys[q], "=")
    tags <- vapply(tag_value, "[", 1, FUN.VALUE=character(1))
    vals <- vapply(tag_value, "[", 2, FUN.VALUE=character(1))
    vals <- gsub(",", "|", vals)
    idx <- 
      lapply(tags, function (tag, val=vals[charmatch(tag, tags)]) {
        if (is.na(val)) {
          ## test for the presence of a qualifier tag among the qualifiers
          ## of each feature
          idx <- vapply(x, function(x_item) {
            any(grepl(tag, names(x_item@qualifiers)))
          }, logical(1L))
        }
        else {
          ## get indices for matching tags for each feature
          tag_idx_lst <- 
            lapply(x, function(x_item) grep(tag, names(x_item@qualifiers)))
          idx <- mapply( function(x_item, tag_idx) { 
            any(grepl(val, x_item@qualifiers[tag_idx]))
          }, x, tag_idx_lst, USE.NAMES=FALSE)
        }   
      })    
    x <- x[ Reduce("&", idx) ]
  }
  x
}

.retrieve <- function (x, cols = "") {
  
  if (!nzchar(cols))
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

isLocation <- function (x) grepl("^loc=|^location=", x)

isKey <- function (x) grepl("^key=", x)

isIndex <- function (x) grepl("^idx=|^index=", x)

isEmpty <- function (x) length(x) == 0L || !nzchar(x)

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