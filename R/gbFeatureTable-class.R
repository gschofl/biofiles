

# gbFeatureTable-class -------------------------------------------------


#' @include gbFeature-class.R
NULL

setOldClass("list")

#' Class \code{"gbFeatureTable"}
#'
#' \dQuote{gbFeatureTable} is an S4 class that provides a container for
#' \dQuote{\linkS4class{gbFeature}}s retrived from GenBank flat files.
#'
#' @slot .seqinfo A \code{\linkS4class{seqinfo}} object containing the
#' genome sequence as an \code{\linkS4class{XStringSet}} object and
#' sequence metadata as a \code{\linkS4class{gbHeader}} object.
#' @slot .id An integer vector the indices of the \code{\linkS4class{gbFeature}}s
#' contained within a \code{gbFeatureTable} object.
#' @slot .Data A list of \code{\linkS4class{gbFeature}} objects.
#' 
#' @export
setClass("gbFeatureTable",
         slots = list(
           .seqinfo = "seqinfo",
           .id = 'integer'
         ),
         contains = "list")


S4Vectors::setValidity2("gbFeatureTable", function(object) {
  if (!all(vapply(object@.Data, is, 'gbFeature', FUN.VALUE = logical(1))))
    return("All elements in a 'gbFeatureTable' must be 'gbFeature' instances")

  TRUE
})


# show -------------------------------------------------------------------


setMethod("show", "gbFeatureTable", function(object) {
  lo <- length(object)
  cat(sprintf("%s with %i features:\n",  sQuote(class(object)), lo))
  if (lo > 0L) {
    show_gbFeature(object[[1L]], showInfo = FALSE, write_to_file = FALSE)
    if (lo > 1L) {
      cat("...\n")
      show_gbFeature(object[[lo]], showInfo = FALSE, write_to_file = FALSE)
    }
  }
  show(.seqinfo(object))
  invisible()
})


# summary ----------------------------------------------------------------


#' @rdname summary-methods
setMethod("summary", "gbFeatureTable", function(object, n = 8, ...) {
  olen <- length(object)
  if (olen > 2*n) {
    hd  <- object[seq_len(n), check = FALSE]
    tl  <- object[seq.int(to = olen, length.out = min(n, olen)), check = FALSE]
    idx  <- c("Id", index(hd), "...", index(tl))
    key  <- c("Feature", key(hd), "...", key(tl))
    loc  <- c(location(hd), "...", location(tl))
    loc  <- c("Location", vapply(loc, as, "character", FUN.VALUE = ""))
    gene <- c("GeneId", geneID(hd), "...", geneID(tl)) %|NA|% "NA"
    prod <- c("Product", product(hd), "...", product(tl)) %|NA|% "NA"
    note <- c("Note", collapse(as.list(note(hd)), '; '), "...", collapse(as.list(note(tl)), '; ')) %|NA|% "NA"
  } else {
    idx <- c("Id", index(object))
    key <- c("Feature", key(object))
    loc <- c("Location", vapply(location(object), as, "character", FUN.VALUE = ""))
    gene <- c("GeneId", geneID(object)) %|NA|% "NA"
    prod <- c("Product", product(object)) %|NA|% "NA"
    note <- c("Note", collapse(as.list(note(object)), '; ')) %|NA|% "NA"
  }
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


setMethod('.seqinfo', 'gbFeatureTable', function(x) {
  x@.seqinfo
})

setMethod('.locus', 'gbFeatureTable', function(x) {
  .locus(.seqinfo(x))
})

setMethod('.header', 'gbFeatureTable', function(x) {
  .header(.seqinfo(x))
})

setMethod('.sequence', 'gbFeatureTable', function(x) {
  .sequence(.seqinfo(x))
})

setMethod('.dbSource', 'gbFeatureTable', function(x) {
  parse_dbsource(getDBSource(x))
})

setMethod(".defline", "gbFeatureTable", function(x) {
  vapply(x, .defline, "", USE.NAMES = FALSE)
})

# getters ----------------------------------------------------------------

#' @rdname accessors
setMethod("getLocus", "gbFeatureTable", function(x) getLocus(.seqinfo(x)) )

#' @rdname accessors
setMethod("getLength", "gbFeatureTable", function(x) getLength(.seqinfo(x)) )

#' @rdname accessors
setMethod("getMoltype", "gbFeatureTable", function(x) getMoltype(.seqinfo(x)) )

#' @rdname accessors
setMethod("getTopology", "gbFeatureTable", function(x) getTopology(.seqinfo(x)) )

#' @rdname accessors
setMethod("getDivision", "gbFeatureTable", function(x) getDivision(.seqinfo(x)) )

#' @rdname accessors
setMethod("getDate", "gbFeatureTable", function(x) getDate(.seqinfo(x)) )

#' @rdname accessors
setMethod("getDefinition", "gbFeatureTable", function(x) getDefinition(.seqinfo(x)) )

#' @rdname accessors
setMethod("getAccession", "gbFeatureTable", function(x) getAccession(.seqinfo(x)) )

#' @rdname accessors
setMethod("getVersion", "gbFeatureTable", function(x) getVersion(.seqinfo(x)) )

#' @rdname accessors
setMethod("getGeneID", "gbFeatureTable", function(x, db = 'gi') getGeneID(.seqinfo(x), db = db) )

#' @rdname accessors
setMethod("getDBLink", "gbFeatureTable", function(x) getDBLink(.seqinfo(x)) )

#' @rdname accessors
setMethod("getDBSource", "gbFeatureTable", function(x) getDBSource(.seqinfo(x)) )

#' @rdname accessors
setMethod("getSource", "gbFeatureTable", function(x) getSource(.seqinfo(x)) )

#' @rdname accessors
setMethod("getOrganism", "gbFeatureTable", function(x) getOrganism(.seqinfo(x)))

#' @rdname accessors
setMethod("getTaxonomy", "gbFeatureTable", function(x) getTaxonomy(.seqinfo(x)))

#' @rdname accessors
setMethod("getReference", "gbFeatureTable", function(x) getReference(.seqinfo(x)))

#' @rdname accessors
setMethod("getKeywords", "gbFeatureTable", function(x) getKeywords(.seqinfo(x)))

#' @rdname accessors
setMethod("getComment", "gbFeatureTable", function(x) getComment(.seqinfo(x)))

#' @rdname getHeader-methods
setMethod("getHeader", "gbFeatureTable", function(x) .header(.seqinfo(x)))

#' @rdname getHeader-methods
setMethod("header", "gbFeatureTable", function(x) .header(.seqinfo(x)))

#' @rdname getSequence-methods
setMethod("getSequence", "gbFeatureTable", function(x) .seq_access(x))

#' @rdname ranges
setMethod("ranges", "gbFeatureTable", function(x, join = FALSE, key = TRUE,
                                              include = "none", exclude = "", ...) {
  .GRanges(x, join = join, include = include, exclude = exclude, key = key)
})

#' @rdname start
setMethod("start", "gbFeatureTable", function(x, join = FALSE) {
  ans <- lapply(x, start, join = join)
  if (join || all(vapply(ans, length, 0) == 1L)) {
    unlist(ans)
  } else {
    ans
  }
})

.gbFeatureTable_replace_start <- function(x, check = TRUE, value) {
  value <- recycle(val = value, len = length(x))
  new_x <- Map(function(Feature, check, val) { 
    start(Feature, check = check) <- val
    Feature
  }, Feature = x, check = list(check), val = value)
  new('gbFeatureTable', .Data = new_x, .id = x@.id, .seqinfo = x@.seqinfo)
}

#' @rdname start
setReplaceMethod("start", "gbFeatureTable", function(x, ..., value)
  .gbFeatureTable_replace_start(x, ..., value = value)
)

#' @rdname end
setMethod("end", "gbFeatureTable", function(x, join = FALSE) {
  ans <- lapply(x, end, join = join)
  if (join || all(vapply(ans, length, numeric(1)) == 1L)) {
    unlist(ans)
  } else {
    ans
  }
})

.gbFeatureTable_replace_end <- function(x, check = TRUE, value) {
  value <- recycle(val = value, len = length(x))
  new_x <- Map(function(Feature, check, val) { 
    end(Feature, check = check) <- val
    Feature
  }, Feature = x, check = list(check), val = value)
  new('gbFeatureTable', .Data = new_x, .id = x@.id, .seqinfo = x@.seqinfo)
}

#' @rdname end
setReplaceMethod("end", "gbFeatureTable", function(x, ..., value)
  .gbFeatureTable_replace_end(x, ..., value = value)
)

#' @rdname strand
setMethod("strand", "gbFeatureTable", function(x, join = FALSE) {
  ans <- lapply(x, strand, join = join)        
  if (join || all(vapply(ans, length, numeric(1)) == 1L)) {
    unlist(ans)
  } else {
    ans
  }
})

#' @rdname strand
setReplaceMethod("strand", "gbFeatureTable", function(x, ..., value) {
  value <- recycle(val = value, len = length(x))
  new_x <- Map(function(Feature, val) { 
    strand(Feature) <- val
    Feature
  }, Feature = x, val = value)
  
  new('gbFeatureTable', .Data = new_x, .id = x@.id, .seqinfo = x@.seqinfo)
})

#' @rdname span
setMethod("span", "gbFeatureTable", function(x, join = FALSE) {
  ans <- lapply(x, span, join = join)
  if (all(vapply(ans, length, 0) == 1L)) {
    unlist(ans)
  } else {
    ans
  }
})

#' @rdname span
setMethod("joint_range", "gbFeatureTable", function(x) {
  do.call("rbind", lapply(x, joint_range))
})

#' @rdname dbxref-methods
setMethod("dbxref", "gbFeatureTable", function(x, db = NULL, ...) {
  dbx <- "db_xref"
  if (!is.null(db)) {
    dbx <- paste0(dbx, ".", db)
  }
  .simplify(.qual_access(x = x, which = dbx, ...), unlist = FALSE)
})

#' @rdname location-methods
setMethod("location", "gbFeatureTable", function(x, join = FALSE) {
  lapply(x, location)
})

#' @rdname fuzzy
setMethod("fuzzy", "gbFeatureTable", function(x) {
  do.call(rbind, lapply(x, fuzzy))
})

#' @rdname index-methods
setMethod("index", "gbFeatureTable", function(x) x@.id)

#' @rdname key-methods
setMethod("key", "gbFeatureTable", function(x) {
  vapply(x, slot, name = 'key', FUN.VALUE = '')
})

#' @rdname qualif-methods
setMethod("qualif", "gbFeatureTable", function(x, which = "", fixed = FALSE, use.names = TRUE) {
  ans <- .qual_access(x, which, fixed, use.names)
  if (use.names) {
    .simplify(ans, unlist = FALSE)
  } else {
    .simplify(ans, unlist = TRUE)
  }
})


# listers ----------------------------------------------------------------


#' @rdname qualifList-methods
setMethod("qualifList", "gbFeatureTable", function(x) {
  lapply(x, qualifList)
})

tbl_qual <- function(x) {
  tblqnm <- Compose(Partial(table, deparse.level = 0), "names", 
                    Partial(slot, name = "qualifiers"))
  lapply(x, tblqnm)
}

tbl_merge <- function(a, b) {
  n <- intersect(names(a), names(b)) 
  res <- c(a[!(names(a) %in% n)], b[!(names(b) %in% n)], a[n] + b[n])
  res[order(names(res))]
}

#' @rdname qualifTable-methods
setMethod("qualifTable", "gbFeatureTable", function(x) {
  tbls <- tbl_qual(x)
  Reduce(tbl_merge, tbls)
})

#' @rdname featureTable-methods
setMethod("featureTable", "gbFeatureTable", function(x) {
  table(key(x))
})


# testers ----------------------------------------------------------------


#' @rdname hasKey-methods
setMethod("hasKey", "gbFeatureTable", function(x, key) {
  vapply(x, hasKey, key, FUN.VALUE = FALSE)
})

#' @rdname hasQualif-methods
setMethod("hasQualif", "gbFeatureTable", function(x, qualifier) {
  vapply(x, hasQualif, qualifier, FUN.VALUE = FALSE)
})


# subsetting ----------------------------------------------------------


#' @rdname extract-methods
setMethod("[", c("gbFeatureTable", "character", "missing", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            check <- list(...)$check %||% TRUE
            i <- which(vapply(x@.Data, slot, name = 'key', FUN.VALUE = '') == i)
            S4Vectors::new2('gbFeatureTable', .Data = x@.Data[i], .id = x@.id[i], 
                            .seqinfo = x@.seqinfo, check = check)
          })

#' @rdname extract-methods
setMethod("[", c("gbFeatureTable", "numeric", "missing", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            check <- list(...)$check %||% TRUE
            S4Vectors::new2('gbFeatureTable', .Data = x@.Data[i], .id = x@.id[i],
                            .seqinfo = x@.seqinfo, check = check)
          })

#' @rdname extract-methods
setMethod("[", c("gbFeatureTable", "logical", "missing", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            check <- list(...)$check %||% TRUE
            S4Vectors::new2('gbFeatureTable', .Data = x@.Data[i], .id = x@.id[i],
                            .seqinfo = x@.seqinfo, check = check)
          })

#' @rdname extract-methods
setMethod("[", c("gbFeatureTable", "missing", "missing", "ANY"),
          function(x, i, j, ..., drop = TRUE) x
          )

#' @rdname extract-methods
setMethod("[[", "gbFeatureTable",
          function(x, i, j, ...) {
            if (missing(i)) {
              stop("subscript is missing")
            }
            if (!is.numeric(i)) {
              stop("invalid subscript type '", class(i), "'")
            }
            if (length(i) > 1) {
              stop("attempt to extract more than one element")
            } 
            res <- callNextMethod()
            ## inject seqinfo
            res@.seqinfo$header <- .header(x)
            res@.seqinfo$sequence <- .sequence(x)
            validObject(res)
            res
          })


# view -------------------------------------------------------------------


#' @rdname view-methods
setMethod("view", "gbFeatureTable", function(x, n)  {
  for (i in x[seq(if (missing(n)) length(x) else n)]) {
    show(i)
    cat("\n")
  }
})


# filter, select, shift, revcomp ----------------------------------------------


#' @rdname manip-methods
setMethod("filter", "gbFeatureTable", function(x, ..., .cols = NULL) {
  .filter(x, ..., .cols = .cols)
})


#' @rdname manip-methods
setMethod("select", "gbFeatureTable", function(x, ..., .cols = NULL) {
  .select(x, ..., .cols = .cols)
})


#' @rdname shift
setMethod("shift", "gbFeatureTable", function(x, shift = 0L, split = FALSE, order = FALSE) {
  .shift(x = x, shift = shift, split = split, order = order)
})


#' @rdname revcomp-methods
setMethod("revcomp", "gbFeatureTable", function(x, order = TRUE) {
  .revcomp(x = x, order = order)
})

