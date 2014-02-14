
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


setValidity2("gbFeatureTable", function(object) {
  if (!all(vapply(object@.Data, is, 'gbFeature', FUN.VALUE = logical(1))))
    return("All elements in a 'gbFeatureTable' must be 'gbFeature' instances")

  TRUE
})


# show -------------------------------------------------------------------


#' @export
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


#' @export
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
    gene <- c("GeneId", geneID(hd), "...", geneID(tl))
    prod <- c("Product", product(hd), "...", product(tl))
    note <- c("Note", collapse(as.list(note(hd)), '; '), "...", collapse(as.list(note(tl)), '; '))
  } else {
    idx <- c("Id", index(object))
    key <- c("Feature", key(object))
    loc <- c("Location", vapply(location(object), as, "character", FUN.VALUE = ""))
    gene <- c("GeneId", geneID(object))
    prod <- c("Product", product(object))
    note <- c("Note", collapse(as.list(note(object)), '; '))
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

#' @export
setMethod("getLocus", "gbFeatureTable", function(x) getLocus(.seqinfo(x)) )
#' @export
setMethod("getLength", "gbFeatureTable", function(x) getLength(.seqinfo(x)) )
#' @export
setMethod("getMoltype", "gbFeatureTable", function(x) getMoltype(.seqinfo(x)) )
#' @export
setMethod("getTopology", "gbFeatureTable", function(x) getTopology(.seqinfo(x)) )
#' @export
setMethod("getDivision", "gbFeatureTable", function(x) getDivision(.seqinfo(x)) )
#' @export
setMethod("getDate", "gbFeatureTable", function(x) getDate(.seqinfo(x)) )
#' @export
setMethod("getDefinition", "gbFeatureTable", function(x) getDefinition(.seqinfo(x)) )
#' @export
setMethod("getAccession", "gbFeatureTable", function(x) getAccession(.seqinfo(x)) )
#' @export
setMethod("getVersion", "gbFeatureTable", function(x) getVersion(.seqinfo(x)) )
#' @export
setMethod("getGeneID", "gbFeatureTable", function(x, db = 'gi') getGeneID(.seqinfo(x), db = db) )
#' @export
setMethod("getDBLink", "gbFeatureTable", function(x) getDBLink(.seqinfo(x)) )
#' @export
setMethod("getDBSource", "gbFeatureTable", function(x) getDBSource(.seqinfo(x)) )
#' @export
setMethod("getSource", "gbFeatureTable", function(x) getSource(.seqinfo(x)) )
#' @export
setMethod("getOrganism", "gbFeatureTable", function(x) getOrganism(.seqinfo(x)))
#' @export
setMethod("getTaxonomy", "gbFeatureTable", function(x) getTaxonomy(.seqinfo(x)))
#' @export
setMethod("getReference", "gbFeatureTable", function(x) getReference(.seqinfo(x)))
#' @export
setMethod("getKeywords", "gbFeatureTable", function(x) getKeywords(.seqinfo(x)))
#' @export
setMethod("getComment", "gbFeatureTable", function(x) getComment(.seqinfo(x)))

#' @export
#' @rdname getHeader-methods
setMethod("getHeader", "gbFeatureTable", function(x) .header(.seqinfo(x)))

#' @export
#' @rdname getHeader-methods
setMethod("header", "gbFeatureTable", function(x) .header(.seqinfo(x)))


#' @export
#' @rdname getSequence-methods
setMethod("getSequence", "gbFeatureTable", function(x) .seq_access(x))


#' @export
setMethod("ranges", "gbFeatureTable", function(x, join = FALSE, key = TRUE,
                                              include = "none", exclude = "") {
  .GRanges(x, join = join, include = include, exclude = exclude, key = key)
})

#' @export
#' @rdname start-methods
setMethod("start", "gbFeatureTable", function(x, join = FALSE, drop = TRUE) {
  ans <- lapply(x, start, join = join, drop = drop)
  if (drop) {
    if (join || all(vapply(ans, length, 0) == 1L)) {
      unlist(ans)
    } else {
      ans
    }
  } else {
    setNames(data.frame(do.call(rbind, ans)), "start")
  }
})

#' @name start<-
#' @export
#' @rdname start-methods
setReplaceMethod("start", "gbFeatureTable", function(x, check = TRUE, value) {
  value <- recycle(value, length(x))
  new_x <- Map(function(Feature, check, val) { 
    start(Feature, check = check) <- val
    Feature
  }, Feature = x, check = list(check), val = value)
  
  new('gbFeatureTable', .Data = new_x, .id = x@.id, .seqinfo = x@.seqinfo)
})

#' @export
#' @rdname end-methods
setMethod("end", "gbFeatureTable", function(x, join = FALSE, drop = TRUE) {
  ans <- lapply(x, end, join = join, drop = drop)
  if (drop) {
    if (join || all(vapply(ans, length, numeric(1)) == 1L)) {
      unlist(ans)
    } else {
      ans
    }
  } else {
    setNames(data.frame(do.call(rbind, ans)), "end")
  }
})

#' @name end<-
#' @export
#' @rdname end-methods
setReplaceMethod("end", "gbFeatureTable", function(x, check = TRUE, value) {
  value <- recycle(value, length(x))
  new_x <- Map(function(Feature, check, val) { 
    end(Feature, check = check) <- val
    Feature
  }, Feature = x, check = list(check), val = value)
  
  new('gbFeatureTable', .Data = new_x, .id = x@.id, .seqinfo = x@.seqinfo)
})

#' @export
#' @rdname strand-methods
setMethod("strand", "gbFeatureTable", function(x, join = FALSE) {
  ans <- lapply(x, strand, join = join)        
  if (join || all(vapply(ans, length, numeric(1)) == 1L)) {
    unlist(ans)
  } else {
    ans
  }
})

#' @name strand<-
#' @export
#' @rdname strand-methods
setReplaceMethod("strand", "gbFeatureTable", function(x, value) {
  value <- recycle(value, length(x))
  new_x <- Map(function(Feature, val) { 
    strand(Feature) <- val
    Feature
  }, Feature = x, val = value)
  
  new('gbFeatureTable', .Data = new_x, .id = x@.id, .seqinfo = x@.seqinfo)
})

#' @export
#' @rdname width-methods
setMethod("width", "gbFeatureTable", function(x) {
  ans <- lapply(x, width)
  if (all(vapply(ans, length, 0) == 1L)) {
    unlist(ans)
  } else {
    ans
  }
})

#' @export
#' @rdname width-methods
setMethod("joint_width", "gbFeatureTable", function(x) {
  unlist(lapply(x, joint_width))
})

#' @export
#' @rdname width-methods
setMethod("joint_range", "gbFeatureTable", function(x) {
  do.call("rbind", lapply(x, joint_range))
})

#' @export
#' @rdname dbxref-methods
setMethod("dbxref", "gbFeatureTable", function(x, db = NULL, ...) {
  dbx <- "db_xref"
  if (!is.null(db)) {
    dbx <- paste0(dbx, ".", db)
  }
  .simplify(.qual_access(x = x, which = dbx, ...), unlist = FALSE)
})

#' @export
#' @rdname location-methods
setMethod("location", "gbFeatureTable", function(x, join = FALSE) {
  lapply(x, location)
})

#' @export
#' @rdname fuzzy-methods
setMethod("fuzzy", "gbFeatureTable", function(x) {
  do.call(rbind, lapply(x, fuzzy))
})

#' @export
#' @rdname index-methods
setMethod("index", "gbFeatureTable", function(x) x@.id)

#' @export
#' @rdname key-methods
setMethod("key", "gbFeatureTable", function(x) {
  vapply(x, slot, name = 'key', FUN.VALUE = '')
})

#' @export
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


#' @export
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

#' @export
#' @rdname qualifTable-methods
setMethod("qualifTable", "gbFeatureTable", function(x) {
  tbls <- tbl_qual(x)
  Reduce(tbl_merge, tbls)
})

#' @export
#' @rdname featureTable-methods
setMethod("featureTable", "gbFeatureTable", function(x) {
  table(key(x))
})


# testers ----------------------------------------------------------------


#' @export
#' @rdname hasKey-methods
setMethod("hasKey", "gbFeatureTable", function(x, key) {
  vapply(x, hasKey, key, FUN.VALUE = FALSE)
})

#' @export
#' @rdname hasQualif-methods
setMethod("hasQualif", "gbFeatureTable", function(x, qualifier) {
  vapply(x, hasQualif, qualifier, FUN.VALUE = FALSE)
})


# subsetting ----------------------------------------------------------


#' @export
#' @rdname extract-methods
setMethod("[", c("gbFeatureTable", "character", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            check <- list(...)$check %||% TRUE
            i <- which(vapply(x@.Data, slot, name = 'key', FUN.VALUE = '') == i)
            IRanges::new2('gbFeatureTable', .Data = x@.Data[i], .id = x@.id[i], 
                          .seqinfo = x@.seqinfo, check = check)
          })

#' @export
#' @rdname extract-methods
setMethod("[", c("gbFeatureTable", "numeric", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            check <- list(...)$check %||% TRUE
            IRanges::new2('gbFeatureTable', .Data = x@.Data[i], .id = x@.id[i],
                          .seqinfo = x@.seqinfo, check = check)
          })

#' @export
#' @rdname extract-methods
setMethod("[", c("gbFeatureTable", "logical", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            check <- list(...)$check %||% TRUE
            IRanges::new2('gbFeatureTable', .Data = x@.Data[i], .id = x@.id[i],
                          .seqinfo = x@.seqinfo, check = check)
          })

#' @export
#' @rdname extract-methods
setMethod("[", c("gbFeatureTable", "missing", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) x
          )


#' @export
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


setMethod("view", "gbFeatureTable", function(x, n)  {
  for (i in x[seq(if (missing(n)) length(x) else n)]) {
    show(i)
    cat("\n")
  }
})


# filter, select, shift, revcomp ----------------------------------------------


#' @export
#' @rdname manip-methods
setMethod("filter", "gbFeatureTable", function(x, ..., .cols = NULL) {
  .filter(x, ..., .cols = .cols)
})


#' @export
#' @rdname manip-methods
setMethod("select", "gbFeatureTable", function(x, ..., .cols = NULL) {
  .select(x, ..., .cols = .cols)
})


#' @export
#' @rdname shift-methods
setMethod("shift", "gbFeatureTable", function(x, shift = 0L, split = FALSE, order = FALSE) {
  .shift(x = x, shift = shift, split = split, order = order)
})


#' @export
#' @rdname revcomp-methods
setMethod("revcomp", "gbFeatureTable", function(x, order = TRUE) {
  .revcomp(x = x, order = order)
})

