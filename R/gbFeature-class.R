#' @include gbLocation-class.R
NULL

#' Class \code{"gbFeature"}
#' 
#' \dQuote{gbFeature} is an S4 class that provides a container
#' for GenBank feature tables.
#' 
#' @slot .seqinfo An \code{\linkS4class{seqinfo}} object containing the
#' full-lenght sequence of the GenBank record that the feature is part
#' of as an \code{\linkS4class{XStringSet}} object, and sequence metadata
#' as a \code{\linkS4class{gbHeader}} object.
#' @slot .id Identifier (index) of the feature in the GenBank record
#' the feature is part of.
#' @slot key The feature key.
#' @slot location A \code{\linkS4class{gbLocation}} object.
#' @slot qualifiers A named character vector. Name attributes
#' correspond to GenBank qualifier tags.
#' 
#' @section Accessor functions:
#' \code{\link{getHeader}}, \code{\link{getSequence}},
#' \code{\link{ranges}}, \code{\link{key}}, \code{\link{index}},
#' \code{\link{qualif}}
#' 
#' @seealso
#'    \code{\linkS4class{gbFeatureTable}}, \code{\linkS4class{gbRecord}}
#'    
#' @export
setClass(
  "gbFeature",
  slots = list(
    .seqinfo  = "seqinfo",
    .id       = "integer",
    key       = "character",
    location  = "gbLocation",
    qualifiers = "character"
  )
)


setValidity2("gbFeature", function(object) {
  TRUE
})


# show -------------------------------------------------------------------


show_gbFeature <- function(object, showInfo = TRUE, write_to_file = FALSE) {
  op <- options("useFancyQuotes")
  options(useFancyQuotes = FALSE)
  on.exit(options(op))
  
  if (write_to_file) {
    ws <- 5       ## added whitespace if we write to file
    width <- 80
  } else {
    ws <- 0
    width <- getOption("width") - 4
    cat("Feature:        Location/Qualifiers:\n")
  }
  
  loc_fmt <- paste0("%s%-16s%s")
  qua_fmt <- paste0("%-16s%s%s = %s")
  loc <- linebreak(as(location(object), "character"), width = width,
                   offset = 17 + ws, indent = 0, split = ",", FORCE = FALSE)
  loc_line <- sprintf(loc_fmt, dup(' ', ws), key(object), loc)
  if (all_empty(object@qualifiers)) {
    qua_line <- ""
  } else {
    qua <- names(object@qualifiers)
    indent <- -(nchar(qua) + 17 + ws + 2)
    val <- unlist(.mapply(linebreak,
                          list(s = dQuote(object@qualifiers), indent = indent),
                          list(width = width, offset = 16 + ws, FORCE = TRUE)))
    qua_line <- sprintf(qua_fmt, "", paste0(dup(' ', ws), "/"), qua, val)
  }
  ft <- paste0(loc_line, "\n", paste0(qua_line, collapse = "\n"))
  
  if (!write_to_file) {
    cat(ft, sep = "\n")
    if (showInfo) {
      show(.seqinfo(object))
    }
  }
  invisible(ft)
}


#' @export
setMethod("show", "gbFeature", function(object) {
  show_gbFeature(object, showInfo = TRUE, write_to_file = FALSE)
})


# summary ----------------------------------------------------------------


#' @export
#' @rdname summary-methods
setMethod("summary", "gbFeature", function(object, ...) {
  idx  <- c("Id", index(object))
  key  <- c("Feature", key(object))
  loc  <- c("Location", as(location(object), "character"))
  gene <- c("GeneId", geneID(object))
  prod <- c("Product", product(object))
  note <- c("Note", collapse(as.list(note(object)), '; '))
  max_idx_len    <- max(nchar(idx))
  max_key_len    <- max(nchar(key))
  max_loc_len    <- max(nchar(loc))
  max_geneid_len <- max(nchar(gene))
  max_prod_len   <- max(nchar(prod))
  fmt <- paste0('%+', max_idx_len + 1, 's %-', max_key_len + 1, 's%-',
                max_loc_len + 1, 's%-', max_geneid_len + 1, 's%-',
                max_prod_len + 1, 's%s')
  showme <- ellipsize(sprintf(fmt, idx, key, loc, gene, prod, note),
                      width = getOption("width") - 3)
  cat(showme, sep = "\n")
  return(invisible(NULL))
})


# Internal getters ----------------------------------------------------------


setMethod('.seqinfo', 'gbFeature', function(x) {
  x@.seqinfo
})

setMethod('.locus', 'gbFeature', function(x) {
  .locus(.seqinfo(x))
})

setMethod('.header', 'gbFeature', function(x) {
  .header(.seqinfo(x))
})

setMethod('.sequence', 'gbFeature', function(x) {
  .sequence(.seqinfo(x))
})

setMethod('.dbSource', 'gbFeature', function(x) {
  parse_dbsource(getDBSource(x))
})

setMethod(".defline", "gbFeature", function(x) {
  paste0("lcl|", key(x), '.', index(x), .dbSource(x), getAccession(x), ' ',
         getDefinition(x))
})


# getters ----------------------------------------------------------------

#' @export
setMethod("getLocus", "gbFeature", function(x) getLocus(.seqinfo(x)) )
#' @export
setMethod("getLength", "gbFeature", function(x) getLength(.seqinfo(x)) )
#' @export
setMethod("getMoltype", "gbFeature", function(x) getMoltype(.seqinfo(x)) )
#' @export
setMethod("getTopology", "gbFeature", function(x) getTopology(.seqinfo(x)) )
#' @export
setMethod("getDivision", "gbFeature", function(x) getDivision(.seqinfo(x)) )
#' @export
setMethod("getDate", "gbFeature", function(x) getDate(.seqinfo(x)) )
#' @export
setMethod("getDefinition", "gbFeature", function(x) getDefinition(.seqinfo(x)) )
#' @export
setMethod("getAccession", "gbFeature", function(x) getAccession(.seqinfo(x)) )
#' @export
setMethod("getVersion", "gbFeature", function(x) getVersion(.seqinfo(x)) )
#' @export
setMethod("getGeneID", "gbFeature", function(x, db = 'gi') getGeneID(.seqinfo(x), db = db) )
#' @export
setMethod("getDBLink", "gbFeature", function(x) getDBLink(.seqinfo(x)) )
#' @export
setMethod("getDBSource", "gbFeature", function(x) getDBSource(.seqinfo(x)) )
#' @export
setMethod("getSource", "gbFeature", function(x) getSource(.seqinfo(x)) )
#' @export
setMethod("getOrganism", "gbFeature", function(x) getOrganism(.seqinfo(x)) )
#' @export
setMethod("getTaxonomy", "gbFeature", function(x) getTaxonomy(.seqinfo(x)) )
#' @export
setMethod("getReference", "gbFeature", function(x) getReference(.seqinfo(x)) )
#' @export
setMethod("getKeywords", "gbFeature", function(x) getKeywords(.seqinfo(x)) )
#' @export
setMethod("getComment", "gbFeature", function(x) getComment(.seqinfo(x)) )

#' @export
#' @rdname getHeader-methods
setMethod("header", "gbFeature", function(x) .header(.seqinfo(x)))

#' @export
#' @rdname getHeader-methods
setMethod("getHeader", "gbFeature", function(x) .header(.seqinfo(x)))

#' @export
#' @rdname getSequence-methods
setMethod("getSequence", "gbFeature", function(x) .seq_access(x))


#' @export
setMethod("ranges", "gbFeature", function(x, include = "none", exclude = "", join = FALSE) {
  .GRanges(x, include = include, exclude = exclude, join = join)
})


#' @export
#' @rdname start-methods
setMethod("start", "gbFeature", function(x, join = FALSE, drop = TRUE) {
  start(x@location, join = join, drop = drop)
})

#' @name start<-
#' @export
#' @rdname start-methods
setReplaceMethod("start", "gbFeature", function(x, check = TRUE, value) {
  start(x@location, check = check) <- value
  if (check) {
    validObject(x)
  }
  x
})

#' @export
#' @rdname end-methods
setMethod("end", "gbFeature", function(x, join = FALSE, drop = TRUE) { 
  end(x@location, join = join, drop = drop)
})

#' @name end<-
#' @export
#' @rdname end-methods
setReplaceMethod("end", "gbFeature", function(x, check = TRUE, value) {
  end(x@location, check = check) <- value
  if (check)
    validObject(x)
  x
})

#' @export
#' @rdname strand-methods
setMethod("strand", "gbFeature", function(x, join = FALSE) {
  strand(x@location, join = join)
})

#' @name strand<-
#' @export
#' @rdname strand-methods
setReplaceMethod("strand", "gbFeature", function(x, value) { 
  strand(x@location) <- value
  x
})

#' @export
#' @rdname width-methods
setMethod("width", "gbFeature", function(x) {
  width(x@location)
})

#' @export
#' @rdname width-methods
setMethod("joint_width", "gbFeature", function(x) {
  joint_width(x@location)
})

#' @export
#' @rdname width-methods
setMethod("joint_range", "gbFeature", function(x) {
  joint_range(x@location)
})

#' @export
#' @rdname dbxref-methods
setMethod("dbxref", "gbFeature", function(x, db = NULL, ...) {
  dbx <- "db_xref"
  if (!is.null(db)) {
    dbx <- paste0(dbx, ".", db)
  }
  .qual_access(x, which = dbx, ...)
})

#' @export
#' @rdname location-methods
setMethod("location", "gbFeature", function(x) x@location)

#' @export
#' @rdname fuzzy-methods
setMethod("fuzzy", "gbFeature", function(x) fuzzy(x@location))

#' @export
#' @rdname index-methods
setMethod("index", "gbFeature", function(x) x@.id)

#' @export
#' @rdname key-methods
setMethod("key", "gbFeature", function(x) structure(x@key, names = NULL) )

#' @name key<-
#' @export
#' @rdname qualif-methods
setReplaceMethod("key", "gbFeature", function(x, check = TRUE, value) {
  x <- initialize(x, key = value)
  if (check)
    validObject(x)
  x
})

#' @export
#' @rdname qualif-methods
setMethod("qualif", "gbFeature", function(x, which, fixed = FALSE, use.names = TRUE) {
  if (missing(which)) {
    x@qualifiers
  } else {
    .qual_access(x, which, fixed, use.names)
  }
})

#' @name qualif<-
#' @export
#' @rdname qualif-methods
setReplaceMethod("qualif", "gbFeature", function(x, which, check = TRUE, value) {
  x@qualifiers[which] <- value
  if (check)
    validObject(x)
  x
})


# listers ----------------------------------------------------------------


#' @export
#' @rdname qualifList-methods
setMethod("qualifList", "gbFeature", function(x) {
  names(x@qualifiers)
})


# testers ----------------------------------------------------------------


#' @export
#' @rdname hasKey-methods
setMethod("hasKey", "gbFeature", function(x, key) {
  !is.na(charmatch(key, x@key))
})


#' @export
#' @rdname hasQualif-methods
setMethod("hasQualif", "gbFeature", function(x, qualifier) {
  !is.na(charmatch(qualifier, names(x@qualifiers)))
})


# shift ---------------------------------------------------------------


#' @export
#' @rdname shift-methods
setMethod("shift", "gbFeature", function(x, shift = 0L, ...) {
  x@location <- shift(x@location, shift)
  x
})


# subsetting ----------------------------------------------------------


#' @export
setMethod("[[", c("gbFeature", "character", "missing"), function(x, i, j) {
  if (i %in% c("key","location",".id")) {
    slot(x, i)
  } else {
    x@qualifiers[i]
  }
})


#' @export
setMethod("$", "gbFeature", function(x, name) {
  if (name %in% c("key","location",".id")) {
    slot(x, name)
  } else {
    x@qualifiers[name]
  }
})

