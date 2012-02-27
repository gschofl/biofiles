## select features by key, location, qualifier values from
## gbData or gbFeatureList objects
.select <- function (x, which=c(""))
{
  #### restict location to the selected window
  if (any(isLocation(which))) {
    start <- as.numeric(which[["start"]])
    end <- as.numeric(which[["end"]])
    strand <- as.numeric(which[["strand"]])
    loc <- lapply(x, getLocation, check=FALSE)
    feature_idx <-
      (vapply(loc, min.start, numeric(1L)) >= if (isEmpty(start)) 1L else start) &
      (vapply(loc, max.end, numeric(1L)) <= if (isEmpty(end)) max(loc[,"end"]) else end) &
      if (isEmpty(strand)) TRUE else vapply(loc, grep.strand, numeric(1L)) == strand
    x <- x[feature_idx] 
  } 
  
  #### within the selected window restrict by key(s) 
  if (any(isKey(which))) {
    key <- paste(which[["key"]], collapse="|")
    key_idx <- grepl(key, vapply(x, function(x_item) x_item@key, character(1L)),
                     ignore.case=TRUE)
    x <- x[key_idx]
  }
  
  #### within resticted window and for restricted key,
  #### search pattern for specified qualifiers 
  if (any(q_idx <- !(isLocation(which) | isKey(which)))) {
    
    idx <- 
      lapply(names(which[q_idx]), function(tag, value=which[q_idx][[tag]]) {
        if (is.logical(value)) {
          ## test for the presence of a qualifier tag among the qualifiers
          ## of each feature
          idx <- vapply(x, function(x_item) {
            any(grepl(tag, names(x_item@qualifiers)))
          }, logical(1L))
          if (isTRUE(value)) idx else !idx
        }
        else {
          ## get indeces for matching tags for each feature
          tag_idx_lst <- 
            lapply(x, function(x_item) grep(tag, names(x_item@qualifiers)))
          idx <- mapply( function(x_item, tag_idx) { 
            any(grepl(value, x_item@qualifiers[tag_idx]))
          }, x, tag_idx_lst, USE.NAMES=FALSE)
        }   
      })
    
    x <- x[Reduce("&", idx)]
  }
  x
}

.retrieve <- function (x, what) {
  
  if (missing(what))
    return(x)
  
  col_names <- character(0)
  
  if (any(idx <- grepl("index", what, ignore.case=TRUE))) {
    i <- .simplify(lapply(x, function (x) x@.ID))
    col_names <- c(col_names, what[idx])
    if (isEmpty(what <- what[!idx]))
      return(.return(i, .Names=col_names))
  } else {
    i <- NULL
  }
  
  if (any(idx <- grepl("key", what, ignore.case=TRUE))) {
    k <- .simplify(lapply(x, function (x) .access(x, where="key")))
    col_names <- c(col_names, what[idx])
    if (isEmpty(what <- what[!idx]))
      return(.return(i, k, .Names=col_names))
  } else {
    k <- NULL
  }

  if (any(idx <- grepl("location|start|end|strand", what, ignore.case=TRUE))) {
    
    if (any(grepl("location", what, ignore.case=TRUE)))
      col_names <- c(col_names, loc <- c("start", "end", "strand"))
    else
      col_names <- c(col_names, loc <- what[idx])
    
    l <- lapply(loc, function (which_loc) {
      .simplify(lapply(x, function (x) .access(x, "location", which_loc)))
    })
    if (length(l) == 1L) l <- l[[1L]]
    
    if (isEmpty(what <- what[!idx]))
      return(.return(i, k, l, .Names=col_names))
  } else {
    l <- NULL
  }

  col_names <- c(col_names, what)
  q <- lapply(what, function (which_qual) {
    .simplify(lapply(x, function (x) 
      .access(x, "qualifiers", which_qual, fixed=TRUE)))
  })
  if (length(q) == 1L) q <- q[[1L]]
  
  .return(i, k, l, q, .Names=col_names)
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


isLocation <- function (x) !is.na(charmatch(names(x), c("start", "end", "strand")))

isKey <- function (x) !is.na(charmatch(names(x), "key"))

isEmpty <- function (x) length(x) == 0L

min.start <- function (x) min(x[grep("start", names(x))])

max.end <- function (x) max(x[grep("end", names(x))])

grep.strand <- function (x) x[grep("strand", names(x))]

hasList <- function (x) vapply(x, typeof, character(1L)) == "list"

hasListOfLists <- function (x) unlist(lapply(x[hasList(x)], hasList))


# --R-- vim:ft=r:sw=2:sts=2:ts=4:tw=76:
#       vim:fdm=marker:fmr={{{,}}}:fdl=0