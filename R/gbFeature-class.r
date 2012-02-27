#### gbFeature objects
setClassUnion("charOrNull", c("character", "NULL"))

##' GenBank feature objects
##'
##' A gbFeature object represents a features parsed from a GenBank flat file
##' record.
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
##' 
##' @examples
##' getSlots("gbFeature")
setClass("gbFeature",
         representation(.Dir="character",
                        .ACCN="character",
                        .DEF="character",
                        .ID="integer",
                        key="character",
                        location="numeric",
                        qualifiers="charOrNull"))

#### Constructor ###########################################################
gbFeature <- function (db_dir, accession, definition, id, key, location, qualifiers)
{
  f <- new("gbFeature", .Dir=as.character(db_dir), .ACCN=as.character(accession),
           .DEF=as.character(definition), .ID=as.integer(id), key=as.character(key),
           location=.getLocation(location), qualifiers=qualifiers)
  f
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

### Accessor methods #######################################################
.access <- function (data, where, which, fixed=FALSE)
{
  if (identical(where, "key"))
    return(structure(data@key, names=NULL))
  
  if (identical(where, "location")) {
    .data <- data@location
    NA_ <- NA_integer_
  } 
  else if (identical(where, "qualifiers")) {
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
.start <- function (data, where="location", what="start") {
  .access(data, where, what)
}

.end <- function (data, where="location", what="end") {
  .access(data, where, what)
}

.strand <- function (data, where="location", what="strand") {
  .access(data, where, what)
}


## get sequence
.seqAccess <- function(s, x, type)
{
  ## initiate empty XStringSet
  template <- switch(type, DNA=DNAStringSet(), AA=AAStringSet(), RNA=RNAStringSet())
  
  if (is(x, "gbFeatureList")) {
    i <- 1
    for (item in x) {
      template <- append(template, subseq(s,
        .access(item, "location", "start"),
        .access(item, "location", "end")))
      template[i]@ranges@NAMES <- paste0(item@.ACCN, ".", item@.ID)
      i <- i + 1
    }
  }
  else if (is(x, "gbFeature")) {
    template <- append(template, subseq(s,
      .access(x, "location", "start"),
      .access(x, "location", "end")))
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
##' @return A numeric vector of feature indeces
##'
##' @export
##' @docType methods
##' @rdname accessor-methods
##'
##' @examples
##' ##
setGeneric("getIndex", function(x, ...) standardGeneric("getIndex"))

##' @aliases index,gbFeature,gbFeature-method
##' @rdname accessor-methods
setMethod("getIndex",
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

### getKey #################################################################
##' @export
##' @docType methods
##' @rdname accessor-methods
setGeneric("getKey", function(x, ...) standardGeneric("getKey"))

##' @aliases getKey,gbFeature,gbFeature-method
##' @rdname accessor-methods
setMethod("getKey",
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

### getLocation ##################################################################
##' @export
##' @docType methods
##' @rdname accessor-methods
setGeneric("getLocation",
           function(x, which=c("start", "end", "strand"), ...)
             standardGeneric("getLocation")
           )

##' @aliases getLocation,gbFeature,gbFeature-method
##' @rdname accessor-methods
setMethod("getLocation",
          signature(x="gbFeature"), 
          function (x, which=c("start", "end", "strand"),
                    attributes=FALSE, check=TRUE)
          {
            if (check && !all(grepl("start|end|strand", which, ignore.case=TRUE)))
              stop("Invalid location identifier. Use 'start', 'end', or 'strand'")

            ans <- .access(x, "location", which)
            if (attributes)
              structure(ans, key=x@key, id=x@.ID,
                        accession=x@.ACCN,
                        definition=x@.DEF,
                        database=x@.Dir)
            else ans
          })

### getQualifier #########################################################
##' @export
##' @docType methods
##' @rdname accessor-methods
setGeneric("getQualifier",
           function(x, which=c(""), ...)
             standardGeneric("getQualifier")
           )

##' @aliases getQualifier,gbFeature,gbFeature-method
##' @rdname accessor-methods
setMethod("getQualifier",
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

### getSequence ############################################################
##' @export
##' @docType methods
##' @rdname accessor-methods
setGeneric("getSequence", function(x, ...) standardGeneric("getSequence"))

##' @aliases getSequence,gbFeature,gbFeature-method
##' @rdname accessor-methods
setMethod("getSequence",
          signature(x="gbFeature"),
          function (x) {
            stopifnot(hasValidDb(x))
            db <- initGenBank(x@.Dir, verbose=FALSE)
            ans <- .seqAccess(dbFetch(db, "sequence"), x, dbFetch(db, "type"))
            ans
          })


### hasKey #################################################################
##' @export
##' @docType methods
##' @rdname accessor-methods
setGeneric("hasKey", function(x, key, ...) standardGeneric("hasKey"))

##' @aliases getKey,gbFeature,gbFeature-method
##' @rdname accessor-methods
setMethod("hasKey", "gbFeature", 
            function (x, key) {
              !is.na(charmatch(key, x@key))
            })

### hasQualifier ############################################################
##' @export
##' @docType methods
##' @rdname accessor-methods
setGeneric("hasQualifier", function(x, qualifier, ...) standardGeneric("hasQualifier"))

##' @aliases hasQualifier,gbFeature,gbFeature-method
##' @rdname accessor-methods
setMethod("hasQualifier", "gbFeature",
          function (x, qualifier) {
              !is.na(charmatch(qualifier,
                               names(.access(x, "qualifiers", what=""))))
          })

################################################################

### show method
##' @export
##' @aliases show,gbFeature,gbFeature-method
##' @rdname gbFeature-class
setMethod("show", signature="gbFeature",
          function (object) {
            op <- options("useFancyQuotes")
            options(useFancyQuotes=FALSE)
            indent <- 16
            pad <- blanks(indent)
            len_feat <- nchar(object@key)
            pad_feat <- blanks(indent - len_feat)
            
            loc <- paste(ifelse(.start(object) == .end(object), 
                                sprintf("%i", .start(object)),
                                sprintf("%i..%i", .start(object), .end(object))),
                         collapse=",")
            
            # if there are more than one start postions use the join(i..i)
            # or order(i..i) syntax
            has_loc_op <- FALSE
            if (object@location["order"] == 1) {
              loc_op <- "order"
              has_loc_op <- TRUE
            } else if (object@location["join"] == 1) {
              loc_op <- "join"
              has_loc_op <- TRUE
            }
            
            if (has_loc_op)
              loc <- sprintf("%s(%s)", loc_op, loc)
            
            # if on the minus strand use the complement(i..i) syntax
            if (object@location["strand"] == -1)
              loc <- sprintf("complement(%s)", loc)
            
            # if necessary wrap the lines
            qua <- names(object@qualifiers)
            loc <- linebreak(loc, offset=indent+1, indent=0,
                             split=",", FORCE=TRUE)
            val <- linebreak(dQuote(object@qualifiers), offset=indent+1, 
                             indent=-(nchar(qua) + 2), FORCE=TRUE)
            
            cat("Feature:         Location/Qualifiers:\n",
                sprintf("%s%s%s\n", object@key, pad_feat, loc),
                sprintf("%s/%s=%s\n", pad, qua, val))
            
            options(op)
            invisible(object)
          })

### Subsetting
##' @export
##' @aliases [[,gbFeature,gbFeature-method
##' @rdname extract-methods
setMethod("[[", signature(x = "gbFeature", i = "character", j = "missing"),
          function(x, i, j) {
            return(slot(object, i))
          })

##' @export
##' @aliases $,gbFeature,gbFeature-method
##' @rdname extract-methods
setMethod("$", signature(x = "gbFeature"),
          function(x, name) {
            return(slot(x, name))
          })
