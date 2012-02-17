#### gbFeature objects
setClassUnion("charOrNull", c("character", "NULL"))

setClass("gbFeature",
         representation(.Dir="character",
                        .ACCN="character",
                        .DEF="character",
                        .ID="integer",
                        key="character",
                        location="numeric",
                        qualifiers="charOrNull"))

#### Constructor
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

#### Accesssor methods
## index
setGeneric("index", function(object, ...) standardGeneric("index"))

setMethod("index", "gbFeature",
          function (object) {
            idx <- object@.ID
            attr(idx, "accession") <- object@.ACCN
            attr(idx, "definition") <- object@.DEF
            attr(idx, "db_dir") <- object@.Dir
            idx
          })

## key
setGeneric("key", function(object, ...) standardGeneric("key"))

setMethod("key", "gbFeature",
          function (object) {
            object@key
          })

## location
setGeneric("start", function(object, ...) standardGeneric("start"))
setGeneric("end", function(object, ...) standardGeneric("end"))
setGeneric("strand", function(object, ...) standardGeneric("strand"))

setMethod("start", "gbFeature", 
          function (object) 
            .location(object, "start")
          )

setMethod("end", "gbFeature", 
          function (object) 
            .location(object, "end")
          )

setMethod("strand", "gbFeature", 
          function (object) 
            .location(object, "strand")
          )

setMethod("length", "gbFeature", 
          function (x) 
            .location(x, "length")
          )

.location <- function (object, what="start")
{
  loc <- object@location
  loc <- loc[grep(what, attr(loc, "names"))]
  names(loc) <- NULL
  loc
}

## qualifiers
setGeneric("qualifiers", function(object, ...) standardGeneric("qualifiers"))
setGeneric("locusTag", function(object, ...) standardGeneric("locusTag"))
setGeneric("gene", function(object, ...) standardGeneric("gene"))
setGeneric("product", function(object, ...) standardGeneric("product"))
setGeneric("note", function(object, ...) standardGeneric("note"))
setGeneric("proteinId", function(object, ...) standardGeneric("proteinId"))
setGeneric("dbxref", function(object, ...) standardGeneric("dbxref"))
setGeneric("translation", function(object, ...) standardGeneric("translation"))
setGeneric("is.pseudo", function(object, ...) standardGeneric("is.pseudo"))

setMethod("qualifiers", "gbFeature", 
          function (object) {
            names(object@qualifiers)
          })

setMethod("locusTag", "gbFeature", 
          function (object) {
            idx <- charmatch("locus_tag", names(object@qualifiers))
            if (is.na(idx))
              return(NA_character_)
            else
              return(object@qualifiers[[idx]])
          })

setMethod("gene", "gbFeature", 
          function (object) {
            idx <- charmatch("gene", names(object@qualifiers))
            if (is.na(idx))
              return(NA_character_)
            else
              return(object@qualifiers[[idx]])
          })

setMethod("product", "gbFeature", 
          function (object) {
            idx <- charmatch("product", names(object@qualifiers))
            if (is.na(idx))
              return(NA_character_)
            else
              return(object@qualifiers[[idx]])
          })

setMethod("note", "gbFeature", 
          function (object) {
            idx <- charmatch("note", names(object@qualifiers))
            if (is.na(idx))
              return(NA_character_)
            else
              return(object@qualifiers[[idx]])
          })

setMethod("proteinId", "gbFeature", 
          function (object) {
            idx <- charmatch("protein_id", names(object@qualifiers))
            if (is.na(idx))
              return(NA_character_)
            else
              return(object@qualifiers[[idx]])
          })

setMethod("dbxref", "gbFeature", 
          function (object) {
            idx <- grep("db_xref", names(object@qualifiers))
            if (any(is.na(idx)))
              return(NA_character_)
            else {
              ans <- object@qualifiers[idx]
              names(ans) <- NULL
              ans
            }
          })

setMethod("translation", "gbFeature", 
          function (object) {
            idx <- charmatch("translation", names(object@qualifiers))
            if (is.na(idx))
              return(NA_character_)
            else
              return(AAString(object@qualifiers[[idx]]))
          })

setMethod("is.pseudo", "gbFeature", 
          function (object) {
            !is.na(charmatch("pseudo", names(object@qualifiers)))
          })

#### show method
setMethod("show", signature="gbFeature",
          function (object) {
            op <- options("useFancyQuotes")
            options(useFancyQuotes=FALSE)
            indent <- 16
            pad <- blanks(indent)
            len_feat <- nchar(object@key)
            pad_feat <- blanks(indent - len_feat)
            
            loc <- paste(ifelse(start(object) == end(object), 
                                sprintf("%i", start(object)),
                                sprintf("%i..%i", start(object), end(object))),
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

## Subsetting
setMethod("[[", signature(x = "gbFeature", i = "character", j = "missing"),
          function(x, i, j) {
            return(slot(object, i))
          })

setMethod("$", signature(x = "gbFeature"),
          function(x, name) {
            return(slot(x, name))
          })
