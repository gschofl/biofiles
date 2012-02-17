## select features by key, location, qualifier values from
## gbData or gbFeatureList objects
.select <- function (x, key, qualifier, location)
{
  
  if (!is.character(key))
    key <- as.character(key)
  
  if (!is.character(qualifier))
    key <- as.character(qualifier)
  
  if (!is.character(location))
    key <- as.character(location)
  
  is.any <- function(x) identical(x, "ANY")
  
  #### first => restict location to the selected window(s)
  if (!grepl("\\d+", location)) location <- "ANY"
  
  if (!is.any(location)) {
    ## parse location
    m <- gregexpr("(?:\\[(?<start>\\d*),(?<end>\\d*)\\])+", location, perl=TRUE)
    start <- as.numeric(mapply( function (str, start, len) {
      substring(str, start, start + len - 1)
    }, location, lapply(m, attr, "capture.start")[[1]][,"start"], 
                                lapply(m, attr, "capture.length")[[1]][,"start"], USE.NAMES=FALSE))
    end <- as.numeric(mapply( function (str, start, len) {
      substring(str, start, start + len - 1)
    }, location, lapply(m, attr, "capture.start")[[1]][,"end"], 
                              lapply(m, attr, "capture.length")[[1]][,"end"], USE.NAMES=FALSE))
    
    if (length(start) != length(end))
      stop("Unequal locations")
    
    ## restrict location
    feature_start_pos <- lapply(x@.Data, .location, "start")
    feature_end_pos <- lapply(x@.Data, .location, "end")
    loc_idx <- rep(seq(length(feature_start_pos)), 
                   vapply(feature_start_pos, length, FUN.VALUE=0))
    loc <- cbind(start=unlist(feature_start_pos), 
                 end=unlist(feature_end_pos), loc_idx)
    feature_idx <- unique(unlist(
      Map(function (start, end) {
        loc[(loc[,"start"] >= if (is.na(start)) 1 else start) & 
          (loc[,"end"] < if (is.na(end)) max(loc[,"end"]) else end),
            "loc_idx"]
      }, start, end))
                          )
    
    x_features <- x[feature_idx] 
  } else { 
    x_features <- x
  }
  
  #### next => if not ANY restrict to key(s) within selected location
  if (!is.any(key)) {
    ## parse key
    key_term <- gsub(",", "|", key)
    
    ## restrict
    key_idx <- grepl(key_term, vapply(x_features, function(x) x@key, FUN.VALUE=""))
    x_features <- x_features[key_idx]
  }
  
  #### next => search pattern for specified qualifiers within resticted
  #### location and for restricted key
  if (!is.any(qualifier)) {
    ## get the qualifier search pattern by matching anything between
    ## square brackets, remove whitespace join comma-separated items or
    ## separate items by |
    m <- strmatch("\\[([^][]*)\\]", qualifier, perl=TRUE)$capture[[1]]
    qual_term <- paste(gsub(",", "|", gsub("[[:space:]]+", "", m)), collapse="|")
    qual_term <- if (nchar(qual_term) == 0) "ANY" else qual_term
    ## get the text value search pattern => anything not surounded by square
    ## brackets
    m <- strmatch("(?<!\\[)(\\b\\w[^][]+\\b)(?!\\])", qualifier, perl=TRUE)$capture[[1]]
    qual_pat <- paste(gsub(",", "|", gsub("[[:space:]]+", " ", m)), collapse="|")
    qual_pat <- if (nchar(qual_pat) == 0) "ANY" else qual_pat
  } else {
    qual_term <- "ANY"
    qual_pat <- "ANY"
  }
  
  ## show which search conditions we are useing
  message(gettextf("Selecting pattern '%s' at qualifier(s) '%s' for key '%s' in location '%s'\n",
                   qual_pat, qual_term, key, location))
  
  ## restrict the set of searched qualifiers
  if (is.any(qual_term) && !is.any(qual_pat)) {
    idx <- vapply(x_features, function(x) any(grepl(qual_pat, x@qualifiers)),
                  FUN.VALUE=TRUE)
    return(x_features[idx])
  } else if (!is.any(qual_term) && !is.any(qual_pat)) {
    qual_idx_lst <- 
      lapply(x_features, function(x) grep(qual_term, names(x@qualifiers)))
    idx <- mapply( function(x, qual_idx) { 
      any(grepl(qual_pat, x@qualifiers[qual_idx])) }, x_features, qual_idx_lst,
                   USE.NAMES=FALSE)
    return(x_features[idx])
  } else {
    return(x_features)
  }
}

# --R-- vim:ft=r:sw=2:sts=2:ts=4:tw=76:
#       vim:fdm=marker:fmr={{{,}}}:fdl=0