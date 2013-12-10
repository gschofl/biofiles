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
#'    \code{\linkS4class{gbFeatureList}}, \code{\linkS4class{gbRecord}}   
#' 
#' @name gbFeature-class
#' @rdname gbFeature-class
#' @exportClass gbFeature
setClass(
  "gbFeature",
  slots = list(
    .seqinfo="seqinfo",
    .id="integer",
    key="character",
    location="gbLocation",
    qualifiers="character"
  )
)


setValidity2("gbFeature", function(object) {
  TRUE
})

# show -------------------------------------------------------------------


show_gbFeature <- function(object, showInfo=TRUE, write_to_file=FALSE) {
  op <- options("useFancyQuotes")
  options(useFancyQuotes=FALSE)
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
  qua_fmt <- paste0("%s%+17s%s=%s")
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
    qua_line <- sprintf(qua_fmt, dup(' ', ws), "/", qua, val)
  }
  ft <- paste0(loc_line, "\n", collapse(qua_line, "\n"))
  
  if (!write_to_file) {
    cat(ft, sep="\n")
    if (showInfo) {
      show(.seqinfo(object))
    }
  }
  invisible(ft)
}


setMethod("show", "gbFeature", function(object) {
  show_gbFeature(object, showInfo=TRUE, write_to_file=FALSE)
})


# summary ----------------------------------------------------------------


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
                      width=getOption("width") - 3)
  cat(showme, sep="\n")
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


setMethod("getLocus", "gbFeature", function(x) getLocus(.seqinfo(x)) )

setMethod("getLength", "gbFeature", function(x) getLength(.seqinfo(x)) )

setMethod("getMoltype", "gbFeature", function(x) getMoltype(.seqinfo(x)) )

setMethod("getTopology", "gbFeature", function(x) getTopology(.seqinfo(x)) )

setMethod("getDivision", "gbFeature", function(x) getDivision(.seqinfo(x)) )

setMethod("getDate", "gbFeature", function(x) getDate(.seqinfo(x)) )

setMethod("getDefinition", "gbFeature", function(x) getDefinition(.seqinfo(x)) )

setMethod("getAccession", "gbFeature", function(x) getAccession(.seqinfo(x)) )

setMethod("getVersion", "gbFeature", function(x) getVersion(.seqinfo(x)) )

setMethod("getGeneID", "gbFeature", function(x, db='gi') getGeneID(.seqinfo(x), db=db) )

setMethod("getDBLink", "gbFeature", function(x) getDBLink(.seqinfo(x)) )

setMethod("getDBSource", "gbFeature", function(x) getDBSource(.seqinfo(x)) )

setMethod("getSource", "gbFeature", function(x) getSource(.seqinfo(x)) )

setMethod("getOrganism", "gbFeature", function(x) getOrganism(.seqinfo(x)) )

setMethod("getTaxonomy", "gbFeature", function(x) getTaxonomy(.seqinfo(x)) )

setMethod("getReference", "gbFeature", function(x) getReference(.seqinfo(x)) )

setMethod("getKeywords", "gbFeature", function(x) getKeywords(.seqinfo(x)) )

setMethod("getComment", "gbFeature", function(x) getComment(.seqinfo(x)) )

#' @export
#' @aliases getHeader,gbFeature-method
#' @rdname getHeader-methods
setMethod("header", "gbFeature", function(x) .header(.seqinfo(x)))

#' @export
#' @aliases getHeader,gbFeature-method
#' @rdname getHeader-methods
setMethod("getHeader", "gbFeature", function(x) .header(.seqinfo(x)))

#' @export
#' @aliases getSequence,gbFeature-method
#' @rdname getSequence-methods
setMethod("getSequence", "gbFeature", function(x) .seq_access(x))


# ' @export
# ' @aliases ranges,gbFeature-method
# ' @rdname ranges-methods
setMethod("ranges", "gbFeature", function(x, include = "none", exclude = "", join = FALSE) {
  .make_GRanges(x, include = include, exclude = exclude, join = join)
})


#' @export
#' @aliases start,gbFeature-method
#' @rdname start-methods
setMethod("start", "gbFeature", function(x, join = FALSE, drop = TRUE) {
  start(x@location, join = join, drop = drop)
})

# ' @export
# ' @aliases start<-,gbFeature-method
# ' @rdname start-methods
setReplaceMethod("start", "gbFeature", function(x, check=TRUE, value) {
  start(x@location, check=check) <- value
  if (check) {
    validObject(x)
  }
  x
})

#' @export
#' @aliases end,gbFeature-method
#' @rdname end-methods
setMethod("end", "gbFeature", function(x, join = FALSE, drop = TRUE) { 
  end(x@location, join = join, drop = drop)
})

# ' @export
# ' @aliases end<-,gbFeature-method
# ' @rdname end-methods
setReplaceMethod("end", "gbFeature", function(x, check=TRUE, value) {
  end(x@location, check=check) <- value
  if (check)
    validObject(x)
  x
})

#' @export
#' @aliases strand,gbFeature-method
#' @rdname strand-methods
setMethod("strand", "gbFeature", function(x, join = FALSE) {
  strand(x@location, join = join)
})

# ' @export
# ' @aliases strand<-,gbFeature-method
# ' @rdname strand-methods
setReplaceMethod("strand", "gbFeature", function(x, value) { 
  strand(x@location) <- value
  x
})

#' @export
#' @aliases width,gbFeature-method
#' @rdname width-methods
setMethod("width", "gbFeature", function(x) {
  width(x@location)
})

#' @export
#' @aliases width,gbFeature-method
#' @rdname width-methods
setMethod("joint_width", "gbFeature", function(x) {
  joint_width(x@location)
})

#' @export
#' @aliases dbxref,gbFeature-method
#' @rdname dbxref-methods
setMethod("dbxref", "gbFeature", function(x, db = NULL, ...) {     
  ans <- .qual_access(x, "db_xref")
  if (all(is.na(ans))) {
    return( NA_character_ )
  } else {
    dbs <- strsplitN(unname(ans), ":", 1, fixed=TRUE)
    ids <- strsplitN(unname(ans), ":", 2, fixed=TRUE)
    if (is.null(db)) {
      setNames(ids, dbs)
    } else {
      db_pattern <- paste0(wrap(db, "\\b"), collapse="|")
      db_pos <- grep(db_pattern, dbs, ignore.case=TRUE)
      if (all_empty(db_pos)) {
        return( NA_character_ )
      } else {
        setNames(ids[db_pos], dbs[db_pos])
      }
    }
  }
})

#' @export
#' @aliases location,gbFeature-method
#' @rdname location-methods
setMethod("location", "gbFeature", function(x) x@location)

#' @export
#' @aliases fuzzy,gbFeature-method
#' @rdname fuzzy-methods
setMethod("fuzzy", "gbFeature", function(x) fuzzy(x@location))

#' @export
#' @aliases index,gbFeature-method
#' @rdname index-methods
setMethod("index", "gbFeature", function(x) x@.id)

#' @export
#' @aliases key,gbFeature-method
#' @rdname key-methods
setMethod("key", "gbFeature", function(x) structure(x@key, names=NULL) )

setReplaceMethod("key", "gbFeature", function(x, check=TRUE, value) {
  x <- initialize(x, key=value)
  if (check)
    validObject(x)
  x
})

#' @export
#' @aliases qualif,gbFeature-method
#' @rdname qualif-methods
setMethod("qualif", "gbFeature", function(x, which, fixed=FALSE, use.names=TRUE) {
  if (missing(which)) {
    x@qualifiers
  } else {
    .qual_access(x, which, fixed, use.names)
  }
})

setReplaceMethod("qualif", "gbFeature", function(x, which, check=TRUE, value) {
  x@qualifiers[which] <- value
  if (check)
    validObject(x)
  x
})


# listers ----------------------------------------------------------------


#' @export
#' @aliases listQualif,gbFeature-method
#' @rdname listQualif-methods
setMethod("listQualif", "gbFeature", function(x) {
  names(x@qualifiers)
})

#' @export
#' @aliases tableQualif,gbFeature-method
#' @rdname tableQualif-methods
setMethod("tableQualif", "gbFeature", function(x) {
  table(names(x@qualifiers))
})

# testers ----------------------------------------------------------------


#' @export
#' @aliases hasKey,gbFeature-method
#' @rdname hasKey-methods
setMethod("hasKey", "gbFeature", function(x, key) {
  !is.na(charmatch(key, x@key))
})


#' @export
#' @aliases hasQualif,gbFeature-method
#' @rdname hasQualif-methods
setMethod("hasQualif", "gbFeature", function(x, qualifier) {
  !is.na(charmatch(qualifier, names(x@qualifiers)))
})


# shift ---------------------------------------------------------------


setMethod("shift", "gbFeature",
          function(x, shift=0L, ...) {
            x@location <- shift(x@location, shift)
            x
          })


# subsetting ----------------------------------------------------------


#' @export
setMethod("[[", c("gbFeature", "character", "missing"),
          function(x, i, j) {
            if (i %in% c("key","location",".id")) {
              slot(x, i)
            } else {
              x@qualifiers[i]
            }
          })


#' @export
setMethod("$", "gbFeature",
          function(x, name) {
            if (name %in% c("key","location",".id")) {
              slot(x, name)
            } else {
              x@qualifiers[name]
            }
          })

