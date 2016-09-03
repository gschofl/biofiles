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
    .seqinfo   = "seqinfo",
    .id        = "integer",
    key        = "character",
    location   = "gbLocation",
    qualifiers = "character"
  )
)

S4Vectors::setValidity2("gbFeature", function(object) {
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


setMethod("show", "gbFeature", function(object) {
  show_gbFeature(object, showInfo = TRUE, write_to_file = FALSE)
})


# summary ----------------------------------------------------------------


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


#' @rdname accessors
setMethod("getLocus", "gbFeature", function(x) getLocus(.seqinfo(x)) )

#' @rdname accessors
setMethod("getLength", "gbFeature", function(x) getLength(.seqinfo(x)) )

#' @rdname accessors
setMethod("getMoltype", "gbFeature", function(x) getMoltype(.seqinfo(x)) )

#' @rdname accessors
setMethod("getTopology", "gbFeature", function(x) getTopology(.seqinfo(x)) )

#' @rdname accessors
setMethod("getDivision", "gbFeature", function(x) getDivision(.seqinfo(x)) )

#' @rdname accessors
setMethod("getDate", "gbFeature", function(x) getDate(.seqinfo(x)) )

#' @rdname accessors
setMethod("getDefinition", "gbFeature", function(x) getDefinition(.seqinfo(x)) )

#' @rdname accessors
setMethod("getAccession", "gbFeature", function(x) getAccession(.seqinfo(x)) )

#' @rdname accessors
setMethod("getVersion", "gbFeature", function(x) getVersion(.seqinfo(x)) )

#' @param db Which database identifier (default: 'gi')
#' @rdname accessors
setMethod("getGeneID", "gbFeature", function(x, db = 'gi') getGeneID(.seqinfo(x), db = db) )

#' @rdname accessors
setMethod("getDBLink", "gbFeature", function(x) getDBLink(.seqinfo(x)) )

#' @rdname accessors
setMethod("getDBSource", "gbFeature", function(x) getDBSource(.seqinfo(x)) )

#' @rdname accessors
setMethod("getSource", "gbFeature", function(x) getSource(.seqinfo(x)) )

#' @rdname accessors
setMethod("getOrganism", "gbFeature", function(x) getOrganism(.seqinfo(x)) )

#' @rdname accessors
setMethod("getTaxonomy", "gbFeature", function(x) getTaxonomy(.seqinfo(x)) )

#' @rdname accessors
setMethod("getReference", "gbFeature", function(x) getReference(.seqinfo(x)) )

#' @rdname accessors
setMethod("getKeywords", "gbFeature", function(x) getKeywords(.seqinfo(x)) )

#' @rdname accessors
setMethod("getComment", "gbFeature", function(x) getComment(.seqinfo(x)) )

#' @rdname getHeader-methods
setMethod("header", "gbFeature", function(x) .header(.seqinfo(x)))

#' @rdname getHeader-methods
setMethod("getHeader", "gbFeature", function(x) .header(.seqinfo(x)))

#' @rdname getSequence-methods
setMethod("getSequence", "gbFeature", function(x) .seq_access(x))

#' @rdname ranges
setMethod("ranges", "gbFeature", function(x, include = "none", exclude = "", join = FALSE) {
  .GRanges(x, include = include, exclude = exclude, join = join)
})

#' @rdname start
setMethod("start", "gbFeature", function(x, join = FALSE) {
  start(x@location, join = join)
})


.gbFeature_replace_start <- function(x, check = TRUE, value) {
  start(x@location, check = check) <- value
  if (check) {
    validObject(x)
  }
  x
}

#' @rdname start
setReplaceMethod("start", "gbFeature", function(x, ..., value)
  .gbFeature_replace_start(x, ..., value = value)
)

#' @rdname end
setMethod("end", "gbFeature", function(x, join = FALSE) { 
  end(x@location, join = join)
})

.gbFeature_replace_end <- function(x, check = TRUE, value) {
  end(x@location, check = check) <- value
  if (check) {
    validObject(x)
  }
  x
}

#' @rdname end
setReplaceMethod("end", "gbFeature", function(x, ..., value)
  .gbFeature_replace_end(x, ..., value = value)
)

#' @rdname strand
setMethod("strand", "gbFeature", function(x, join = FALSE) {
  strand(x@location, join = join)
})

#' @rdname strand
setReplaceMethod("strand", "gbFeature", function(x, ..., value) { 
  strand(x@location, ...) <- value
  x
})

#' @rdname span
setMethod("span", "gbFeature", function(x, join = FALSE) {
  span(x@location, join = join)
})

#' @rdname span
setMethod("joint_range", "gbFeature", function(x) {
  joint_range(x@location)
})

#' @rdname dbxref-methods
setMethod("dbxref", "gbFeature", function(x, db = NULL, ...) {
  dbx <- "db_xref"
  if (!is.null(db)) {
    dbx <- paste0(dbx, ".", db)
  }
  .qual_access(x, which = dbx, ...)
})

#' @rdname location-methods
setMethod("location", "gbFeature", function(x) x@location)

#' @rdname fuzzy
setMethod("fuzzy", "gbFeature", function(x) fuzzy(x@location))

#' @rdname index-methods
setMethod("index", "gbFeature", function(x) x@.id)

#' @rdname key-methods
setMethod("key", "gbFeature", function(x) structure(x@key, names = NULL) )

#' @rdname key-methods
setReplaceMethod("key", "gbFeature", function(x, check = TRUE, value) {
  x <- initialize(x, key = value)
  if (check)
    validObject(x)
  x
})

#' @rdname qualif-methods
setMethod("qualif", "gbFeature", function(x, which, fixed = FALSE, use.names = TRUE) {
  if (missing(which)) {
    x@qualifiers
  } else {
    .qual_access(x, which, fixed, use.names)
  }
})

#' @rdname qualif-methods
setReplaceMethod("qualif", "gbFeature", function(x, which, check = TRUE, value) {
  assertthat::assert_that(!missing(which))
  x@qualifiers[which] <- value
  if (check)
    validObject(x)
  x
})


# listers ----------------------------------------------------------------


#' @rdname qualifList-methods
setMethod("qualifList", "gbFeature", function(x) {
  names(x@qualifiers)
})


# testers ----------------------------------------------------------------


#' @rdname hasKey-methods
setMethod("hasKey", "gbFeature", function(x, key) {
  !is.na(charmatch(key, x@key))
})


#' @rdname hasQualif-methods
setMethod("hasQualif", "gbFeature", function(x, qualifier) {
  !is.na(charmatch(qualifier, names(x@qualifiers)))
})


# shift ---------------------------------------------------------------


#' @rdname shift
setMethod("shift", "gbFeature", function(x, shift = 0L, ...) {
  x@location <- shift(x@location, shift)
  x
})


# subsetting ----------------------------------------------------------


#' @rdname extract-methods
setMethod("[[", c("gbFeature", "character", "missing"), function(x, i, j) {
  if (i %in% c("key", "location", ".id")) {
    slot(x, i)
  } else {
    x@qualifiers[i]
  }
})

#' @param name The name of the element to extract.
#' @rdname extract-methods
setMethod("$", "gbFeature", function(x, name) {
  if (name %in% c("key", "location", ".id")) {
    slot(x, name)
  } else {
    x@qualifiers[name]
  }
})

