#### gbFeature-class ####

setClassUnion("charOrNull", c("character", "NULL"))

##' gbFeature class
##' 
##' \sQuote{gbFeature} is an S4 class that extends the
##' \code{\link{gbFeatureList-class}}. This class provides a container
##' for feature data retrived from GenBank flat files.
##' 
##' \describe{
##'    \item{.Dir}{The path to the database file containing the GenBank
##'    record the feature is part of.}
##'    \item{.ACCN}{Accession number of the GenBank record that the
##'    feature is part of.}
##'    \item{.DEF}{The definition line (brief description of the sequence)
##'    of the GenBank record the feature is part of.}
##'    \item{.ID}{Identifier (sequential index) of the feature in the
##'    GenBank record the feature is part of.}
##'    \item{key}{Feature key (e.g. Source, CDS, gene, etc.)}
##'    \item{location}{Named numeric vector. Name attributes: 'start<N>',
##'    'end<N>', 'length<N>', 'strand', 'join', 'order'.}
##'    \item{qualifiers}{Named character vector. Name attributes
##'    correspond to GenBank qualifier tags.}       
##' }
##' 
##' @name gbFeature-class
##' @rdname gbFeature-class
##' @exportClass gbFeature
##' @aliases show,gbFeature-method
##' @aliases getIndex,gbFeature-method
##' @aliases getKey,gbFeature-method
##' @aliases getLocation,gbFeature-method
##' @aliases getQualifier,gbFeature-method
##' @aliases getSequence,gbFeature-method
##' @aliases hasKey,gbFeature-method
##' @aliases hasQualifier,gbFeature-method
##' @aliases [[,gbFeature-method
##' @aliases $,gbFeature-method
.gbFeature <- 
  setClass("gbFeature",
           representation(.Dir="character",
                          .ACCN="character",
                          .DEF="character",
                          .ID="integer",
                          key="character",
                          location="gbLocation",
                          qualifiers="charOrNull"))

##' @export
setMethod("show", 
          #### show-method ####
          signature("gbFeature"),
          function (object) {
            op <- options("useFancyQuotes")
            options(useFancyQuotes=FALSE)
            indent <- 16
            pad <- blanks(indent)
            len_feat <- nchar(object@key)
            pad_feat <- blanks(indent - len_feat)
            
            # if necessary wrap the lines
            qua <- names(object@qualifiers)
            loc <- linebreak(as(object@location, "character"),
                             offset=indent+1, indent=0, split=",", FORCE=TRUE)
            val <- linebreak(dQuote(object@qualifiers), offset=indent+1, 
                             indent=-(nchar(qua) + 2), FORCE=TRUE)
            
            cat("Feature:         Location/Qualifiers:\n",
                sprintf("%s%s%s\n", object@key, pad_feat, loc),
                sprintf("%s/%s=%s\n", pad, qua, val))
            
            options(op)
            invisible(object)
          })


# Constructor ---------------------------------------------------------

gbFeature <- function (db_dir, accession, definition, id, key, location, qualifiers) 
{
  .gbFeature(.Dir=as.character(db_dir),
             .ACCN=as.character(accession),
             .DEF=as.character(definition),
             .ID=as.integer(id),
             key=as.character(key),
             location=.getLocationS4(location),
             qualifiers=qualifiers)
}

## Extract location information from a genBank base span line
.getLocation <- function(gb_base_span)
{
  # transforms location information in the GenBank format (e.g. 1..23,
  # or complement(join(345..543,567..567)) into a named vector
  # start(1, start2), end(1, end2), length(1, length2), strand (1,-1), length
  strand <- ifelse(grepl("complement", gb_base_span), -1, 1)
  join <- ifelse(grepl("join", gb_base_span), 1, 0)
  order <- ifelse(grepl("order", gb_base_span), 1, 0)
  split_loc <- strsplit(unlist(strsplit(gsub("[^0-9\\.,]+", "", gb_base_span), ",")), "\\.\\.")
  start <- as.numeric(lapply(split_loc, "[", 1))
  end <- as.numeric(lapply(split_loc, "[", 2))
  end[is.na(end)] <- start[is.na(end)]
  length <- end - start + 1
  loc <- c(start=start, end=end, length=length, strand=strand, join=join, order=order)
  loc
}

# test cases
# gb_base_span <- "340..565"
# gb_base_span <- "<340..565"
# gb_base_span <- "566..>567"
# gb_base_span <- "102.110"
# gb_base_span <- "123^124"
# gb_base_span <- "J00194.1:100..202"
# gb_base_span <- "complement(565..>567)"
# gb_base_span <- "join(345..543,567..>590)"
# gb_base_span <- "order(<345..543,<567..>569,666..7000)"
# gb_base_span <- "join(complement(4918..5163),complement(2691..4571),7665..7899)"
# gb_base_span <- "complement(join(345..543,AL121804.2:567..>569,AL121804.2:<600..603))"

.getLocationS4 <- function(gb_base_span)
{                       
  # single location
  sil <- "^\\d+$"
  
  # between location
  bl <- "\\d+\\^\\d+"
  # within location
  wl <- "[<]?\\d+\\.[>]?\\d+"
  # paired location
  pl <- "[<]?\\d+\\.\\.[>]?\\d+"
  
  # simple location
  sl <- sprintf("([a-zA-z][a-zA-Z0-9]*(\\.[a-zA-Z0-9]+)?\\:)?(%s|%s|%s)", bl, wl, pl)
  # complemented simple location
  csl <- sprintf("complement\\(%s\\)", sl)
  # possibly complemented simplex location
  pcsl <- sprintf("(%s|%s)", sl, csl)
  
  
  # remote accession
  ra <- "([a-zA-Z][a-zA-Z0-9]*(\\.[a-zA-Z0-9]+)?)" 
  # compound location
  cl <- sprintf("(join|order)\\(%s(,%s)*\\)", pcsl, pcsl)
  # complemented compound location
  ccl <- sprintf("complement\\(%s\\)", cl)
  
  .parseSimpleSpan <- function (base_span) {  ## test for strand
    strand <- ifelse(grepl(csl, base_span), -1L, 1L)
    ## get span string
    span_str <- str_extract(base_span,  sl)
    ## get remote accession number
    accn <- str_extract(span_str, ra)
    remote <- ifelse(!is.na(accn), TRUE, FALSE)
    ## get closed and span
    span <- gsub(paste0(ra, "\\:"), "", span_str)
    closed <- ifelse(grepl(wl, span), FALSE, TRUE)
    span <- do.call(rbind, strsplit(span, "\\.\\.|\\.|\\^"))
    ## get partial
    partial <- matrix(grepl("^(<|>)", span), ncol=2)
    span <- matrix(as.integer(gsub("^(<|>)", "", span)), ncol=2)
    return(list(span=span, strand=strand, partial=partial, accn=accn,
                remote=remote, closed=closed))
  }
  
  # test for possibly complemented simple location first
  if (str_detect(gb_base_span, sprintf("^%s$", pcsl))) {
    l <- .parseSimpleSpan(gb_base_span)
    return(.gbLocation(.Data=l$span, strand=l$strand,
                       compound=NA_character_, partial=l$partial,
                       accession=l$accn, remote=l$remote,
                       closed=l$closed))
  }
  # test for possibly complemented compound
  else if (str_detect(gb_base_span, cl)) {
    ## test for complementary strand
    strand <- ifelse(grepl(ccl, gb_base_span), -1L, 1L)
    ## get compound
    cmpnd_str <- str_extract(gb_base_span, cl)
    compound <- str_extract(cmpnd_str, "(join|order)")
    ## get span strings
    span_str <- strsplit(str_extract(cmpnd_str, sprintf("%s(,%s)*", pcsl, pcsl)), ",")[[1L]]
    l <- lapply(span_str, .parseSimpleSpan)
    
    if (any(vapply(l, "[[", "strand", FUN.VALUE=integer(1)) < 1L)) {
      strand <- vapply(l, "[[", "strand", FUN.VALUE=integer(1))
    }
    
    return(.gbLocation(.Data=do.call(rbind, lapply(l, "[[", "span")), 
                       strand=strand, compound=compound,
                       partial=do.call(rbind, lapply(l, "[[", "partial")),
                       accession=vapply(l, "[[", "accn", FUN.VALUE=character(1)),
                       remote=vapply(l, "[[", "remote", FUN.VALUE=logical(1))))
  }
  else {
    return(.gbLocation())
  }
}


# Accessors -----------------------------------------------------------

.access <- function (data, where, which, fixed=FALSE)
{
  if (identical(where, "key"))
    return(structure(data@key, names=NULL))
  
  if (identical(where, "qualifiers")) {
    .data <- data@qualifiers
    NA_ <- NA_character_
  }
  
  if (fixed) which <- paste0("\\<", which, "\\>")
  
  n_row <- length(.data)
  idx <- matrix(
    vapply(which, grepl, x=names(.data),
           USE.NAMES=FALSE, FUN.VALUE=logical(n_row)),
    nrow=n_row)

  if (ncol(idx) == 1L) {
    ans <- structure(.data[idx], names=NULL)
    
    if (length(ans) > 0L) 
      return(ans)
    else
      return(structure(NA_, names=NULL))
  } 
  else if (ncol(idx) > 1L) {
    ans <- lapply(seq.int(ncol(idx)), function (i) .data[idx[,i]])
    
    if (any(na_idx <- vapply(ans, length, FUN.VALUE=integer(1L)) == 0L))
      for (na in which(na_idx))
        ans[[na]] <- structure(NA_, names=which[na])
    
    return(unlist(ans))
  }
}

## some internal shortcuts
.start <- function (x, drop=FALSE) {
  x@location@.Data[, 1, drop=drop]
}

.end <- function (x, drop=FALSE) {
  x@location@.Data[, 2, drop=drop]
}

.strand <- function (x) {
  x@location@strand
}

## access sequence
.seqAccess <- function(s, x, type)
{
  ## initiate empty XStringSet
  template <- switch(type, DNA=DNAStringSet(), AA=AAStringSet(), RNA=RNAStringSet())
  
  if (is(x, "gbFeatureList")) {
    i <- 1
    for (item in x) {
      template <- append(template, subseq(s,
                                          .start(item, drop=TRUE),
                                          .end(item, drop=TRUE)))
      template[i]@ranges@NAMES <- paste0(item@.ACCN, ".", item@.ID)
      i <- i + 1
    }
  }
  else if (is(x, "gbFeature")) {
    template <- append(template, subseq(s,
                                        .start(item, drop=TRUE),
                                        .end(item, drop=TRUE)))
    template@ranges@NAMES <- paste0(x@.ACCN, ".", x@.ID)
  }
  
  template@metadata <- list(definition=x@.DEF, database=x@.Dir)
  template
}

##' Retrieve index of a GenBank feature.
##'
##' @usage index(data)
##'
##' @param data An instance of \code{\link{gbFeature-class}} or
##' \code{\link{gbFeatureList-class}}
##'
##'@family accessors
##' @return A numeric vector of feature indeces
##'
##' @docType methods
##' @rdname gbFeature-class
##' @export
setGeneric("getIndex",
           #### getIndex-generic ####
           function(x, ...) {
             standardGeneric("getIndex")
           })


##' @export
setMethod("getIndex",
          #### getIndex-method ####
          signature(x="gbFeature"),
          function (x, attributes=FALSE) {
            ans <- x@.ID
            if (attributes)
                structure(ans,
                          accession=x@.ACCN,
                          definition=x@.DEF,
                          database=x@.Dir)
            else ans
          })

##' @docType methods
##' @rdname gbFeature-class
##' @export
setGeneric("getKey", 
           #### getKey-generic ####
           function(x, ...) {
             standardGeneric("getKey")
             })


##' @export
setMethod("getKey",
          #### getKey-method ####
          signature(x="gbFeature"), 
          function (x, attributes=FALSE) {
            ans <- .access(x, where="key")
            if (attributes)
              structure(ans, id=x@.ID,
                        accession=x@.ACCN,
                        definition=x@.DEF,
                        database=x@.Dir)
            else ans
          })

##' @docType methods
##' @rdname gbFeature-class
##' @export
setGeneric("getLocation",
           #### getLocation-generic ####
           function(x, which=c("start", "end", "strand"), ...) {
             standardGeneric("getLocation")
           })

##' @export
setMethod("getLocation",
          #### getLocation-method ####
          signature(x="gbFeature"),
          function (x, attributes=FALSE) {  
            ans <- cbind(x@location@.Data, x@location@strand)
            colnames(ans) <- c("start", "end", "strand")
            if (attributes) {
              structure(ans, key=x@key, id=x@.ID,
                        accession=x@.ACCN,
                        definition=x@.DEF,
                        database=x@.Dir)
            } else {
              ans
            }
          })

##' @docType methods
##' @rdname gbFeature-class
##' @export
setGeneric("getQualifier",
           #### getQualifier-generic ####
           function(x, which=c(""), ...) {
             standardGeneric("getQualifier")
           })

##' @export
setMethod("getQualifier",
          #### getQualifier-method ####
          signature(x="gbFeature"), 
          function (x, which=c(""), attributes=FALSE, fixed=FALSE) {
            ans <- .access(x, "qualifiers", which, fixed)
            if (attributes)
              structure(ans, key=x@key, id=x@.ID,
                accession=x@.ACCN,
                definition=x@.DEF,
                database=x@.Dir)
            else ans
          })

##' @docType methods
##' @rdname gbFeature-class
##' @export
setGeneric("getSequence",
           #### getSequence-generic ####
           function(x, ...) {
             standardGeneric("getSequence")
             })

##' @export
setMethod("getSequence",
          #### getSequence-method ####
          signature(x="gbFeature"),
          function (x) {
            stopifnot(hasValidDb(x))
            db <- initGenBank(x@.Dir, verbose=FALSE)
            ans <- .seqAccess(dbFetch(db, "sequence"), x, dbFetch(db, "type"))
            ans
          })

##' @docType methods
##' @rdname gbFeature-class
##' @export
setGeneric("hasKey",
           #### hasKey-generic ####
           function(x, key, ...) {
             standardGeneric("hasKey")
           })

##' @export
setMethod("hasKey", 
          #### hasKey-method ####
          signature(x="gbFeature"), 
          function (x, key) {
            !is.na(charmatch(key, x@key))
          })

##' @docType methods
##' @rdname gbFeature-class
##' @export
setGeneric("hasQualifier", 
           #### hasQualifier-generic ####
           function(x, qualifier, ...) {
             standardGeneric("hasQualifier")
           })

##' @export
setMethod("hasQualifier",
          #### hasQualifier-method ####
          signature("gbFeature"),
          function (x, qualifier) {
            !is.na(charmatch(qualifier,
                             names(.access(x, "qualifiers", what=""))))
          })


# Subsetting ----------------------------------------------------------


##' @export
setMethod("[[",
          #### [[-method ####
          signature(x = "gbFeature", i = "character", j = "missing"),
          function(x, i, j) {
            return(slot(object, i))
          })

##' @export
setMethod("$",
          #### $-method ####
          signature(x = "gbFeature"),
          function(x, name) {
            return(slot(x, name))
          })
