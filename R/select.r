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


## select features by index, key, location, qualifier values from
## gbData or gbFeatureList objects
.select <- function (x, which="")
{
  which <- gsub("\n|\t", " ", which)
  which_list <- strsplit(which, ";")[[1]]
  
  #### restrict to selected indices
  if (any(i <- isIndex(which_list))) {
    indices <- numeric(0)
    for (term in which_list[i]) {
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
    x <- x[ matchIdx(x, indices) ]
  }
  
  #### restrict to selected location
  if (any(l <- isLocation(which_list))) {
    start <- end <- numeric(0)
    for (term in which_list[l]) {
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
  if (any(k <- isKey(which_list))) {
    key <- character(0)
    for (term in which_list[k]) {
      key <- c(key, gsub(",", "|", strsplit(term, "=")[[1]][2]))
    }
    key <- paste0(key, collapse="|")
    key_idx <- grepl(key, getKey(x, attributes=FALSE))
    x <- x[ key_idx ]
  }
  
  #### within resticted window and for restricted key,
  #### search pattern for specified qualifiers 
  if (any(q <- !( isIndex(which_list) | 
                  isLocation(which_list) | 
                  isKey(which_list)))) {
    
    tag_value <- strsplit(which_list[q], "=")
    tags <- vapply(tag_value, "[", 1, FUN.VALUE=character(1))
    vals <- vapply(tag_value, "[", 2, FUN.VALUE=character(1))
    vals <- gsub(",","|",vals)
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

.retrieve <- function (x, which="") {
  
  if (!nzchar(which)) return( x )
  
  which <- gsub("\n|\t", " ", which)
  which_list <- strsplit(which, ";")[[1]]
  col_names <- character(0)
  
  if (any(idx <- grepl("idx|index", which_list, ignore.case=TRUE))) {
    i <- vapply(x, function (x) x@.ID, FUN.VALUE=numeric(1))
    col_names <- c(col_names, which_list[idx])
    if (isEmpty(which_list <- which_list[!idx]))
      return( .return(i, .Names=col_names) )
  } else {
    i <- NULL
  }
  
  if (any(idx <- grepl("key", which_list, ignore.case=TRUE))) {
    k <- vapply(x, function (x) x@key, FUN.VALUE=character(1))
    col_names <- c(col_names, which_list[idx])
    if (isEmpty(which_list <- which_list[!idx]))
      return( .return(i, k, .Names=col_names) )
  } else {
    k <- NULL
  }

  if (any(idx <- grepl("location|start|end|strand", which_list, ignore.case=TRUE))) {
    
    if (any(grepl("location", which_list, ignore.case=TRUE))) {
      col_names <- c(col_names, c("start", "end", "width", "strand"))
      l <- range(x)
    } else {
      col_names <- c(col_names, loc <- which_list[idx])
      l <- lapply(loc, function (loc) {
        switch(loc, start=start(x), end=end(x), strand=strand(x))
      })
    }

    if (isEmpty(which_list <- which_list[!idx]))
      return( .return(i, k, l, .Names=col_names) )
  } else {
    l <- NULL
  }

  col_names <- c(col_names, which_list)
  q <- lapply(which_list, function (which_qual) {
    .simplify(lapply(x, function (x) 
      .access(x, which_qual, fixed=TRUE)))
  })
  if (length(q) == 1L) q <- q[[1L]]
  
  .return( i, k, l, q, .Names=col_names )
}


.return <- function (..., .Names) {
  args <- list(...)
  zero <- vapply(args, is.null, logical(1L))
  args[zero] <- NULL
  if (!any(hasList(args)) && length(args) == 1L)
    return(args[[1L]])
  else if (!any(hasListOfLists(args))) {
    args <- flatten(args)
    structure(
      data.frame(stringsAsFactors=FALSE, args),
      names=.Names)
  }
  else {
    args <- flatten(args, stop.at=2)
    structure(args, names=.Names)
  }
}

isLocation <- function (x) grepl("^loc=|^location=", x)

isKey <- function (x) grepl("^key=", x)

isIndex <- function (x) grepl("^idx=|^index=", x)

isEmpty <- function (x) length(x) == 0L || !nzchar(x)

hasList <- function (x) vapply(x, typeof, character(1L)) == "list"

hasListOfLists <- function (x) unlist(lapply(x[hasList(x)], hasList))

matchIdx <- function(x, idx) {
  x_idx <- vapply(x, function(f) f@.ID, FUN.VALUE=numeric(1))
  which(x_idx %in% idx)  
}


# --R-- vim:ft=r:sw=2:sts=2:ts=4:tw=76:
#       vim:fdm=marker:fmr={{{,}}}:fdl=0