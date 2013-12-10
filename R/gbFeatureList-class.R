
# gbFeatureList-class -------------------------------------------------

#' @include gbFeature-class.R
NULL

setOldClass("list")

#' Class \code{"gbFeatureList"}
#'
#' \dQuote{gbFeatureList} is an S4 class that provides a container for
#' \dQuote{\linkS4class{gbFeature}}s retrived from GenBank flat files.
#'
#' @slot .seqinfo A \code{\linkS4class{seqinfo}} object containing the
#' genome sequence as an \code{\linkS4class{XStringSet}} object and
#' sequence metadata as a \code{\linkS4class{gbHeader}} object.
#' @slot .Data A list of \code{\linkS4class{gbFeature}} objects.
#' 
#' @export
setClass("gbFeatureList",
         slots=list(.seqinfo="seqinfo"),
         contains="list")


setValidity2("gbFeatureList", function(object) {
  if (!all(vapply(object@.Data, is, 'gbFeature', FUN.VALUE=logical(1))))
    return("All elements in a 'gbFeatureList' must be 'gbFeature' instances")

  TRUE
})


# show -------------------------------------------------------------------


#' @export
setMethod("show", "gbFeatureList", function(object) {
  lo <- length(object)
  cat(sprintf("%s with %i features:\n",  sQuote(class(object)), lo))
  if (lo > 0L) {
    show_gbFeature(object[[1L]], showInfo=FALSE, write_to_file=FALSE)
    if (lo > 1L) {
      cat("...\n")
      show_gbFeature(object[[lo]], showInfo=FALSE, write_to_file=FALSE)
    }
  }
  show(.seqinfo(object))
  invisible()
})


# summary ----------------------------------------------------------------


#' @export
#' @rdname summary-methods
setMethod("summary", "gbFeatureList", function(object, n=8, ...) {
  olen <- length(object)
  if (olen > 2*n) {
    hd  <- object[seq_len(n), check=FALSE]
    tl  <- object[seq.int(to = olen, length.out = min(n, olen)), check=FALSE]
    idx  <- c("Id", index(hd), "...", index(tl))
    key  <- c("Feature", key(hd), "...", key(tl))
    loc  <- c(location(hd), "...", location(tl))
    loc  <- c("Location", vapply(loc, as, "character", FUN.VALUE=""))
    gene <- c("GeneId", geneID(hd), "...", geneID(tl))
    prod <- c("Product", product(hd), "...", product(tl))
    note <- c("Note", collapse(as.list(note(hd)), '; '), "...", collapse(as.list(note(tl)), '; '))
  } else {
    idx <- c("Id", index(object))
    key <- c("Feature", key(object))
    loc <- c("Location", vapply(location(object), as, "character", FUN.VALUE=""))
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
                      width=getOption("width") - 3)
  cat(showme, sep="\n")
  return(invisible(NULL))
})


# Internal getters ----------------------------------------------------------


setMethod('.seqinfo', 'gbFeatureList', function(x) {
  x@.seqinfo
})

setMethod('.locus', 'gbFeatureList', function(x) {
  .locus(.seqinfo(x))
})

setMethod('.header', 'gbFeatureList', function(x) {
  .header(.seqinfo(x))
})

setMethod('.sequence', 'gbFeatureList', function(x) {
  .sequence(.seqinfo(x))
})

setMethod('.dbSource', 'gbFeatureList', function(x) {
  parse_dbsource(getDBSource(x))
})

setMethod(".defline", "gbFeatureList", function(x) {
  vapply(x, .defline, "", USE.NAMES=FALSE)
})


# getters ----------------------------------------------------------------

#' @export
setMethod("getLocus", "gbFeatureList", function(x) getLocus(.seqinfo(x)) )
#' @export
setMethod("getLength", "gbFeatureList", function(x) getLength(.seqinfo(x)) )
#' @export
setMethod("getMoltype", "gbFeatureList", function(x) getMoltype(.seqinfo(x)) )
#' @export
setMethod("getTopology", "gbFeatureList", function(x) getTopology(.seqinfo(x)) )
#' @export
setMethod("getDivision", "gbFeatureList", function(x) getDivision(.seqinfo(x)) )
#' @export
setMethod("getDate", "gbFeatureList", function(x) getDate(.seqinfo(x)) )
#' @export
setMethod("getDefinition", "gbFeatureList", function(x) getDefinition(.seqinfo(x)) )
#' @export
setMethod("getAccession", "gbFeatureList", function(x) getAccession(.seqinfo(x)) )
#' @export
setMethod("getVersion", "gbFeatureList", function(x) getVersion(.seqinfo(x)) )
#' @export
setMethod("getGeneID", "gbFeatureList", function(x, db='gi') getGeneID(.seqinfo(x), db=db) )
#' @export
setMethod("getDBLink", "gbFeatureList", function(x) getDBLink(.seqinfo(x)) )
#' @export
setMethod("getDBSource", "gbFeatureList", function(x) getDBSource(.seqinfo(x)) )
#' @export
setMethod("getSource", "gbFeatureList", function(x) getSource(.seqinfo(x)) )
#' @export
setMethod("getOrganism", "gbFeatureList", function(x) getOrganism(.seqinfo(x)))
#' @export
setMethod("getTaxonomy", "gbFeatureList", function(x) getTaxonomy(.seqinfo(x)))
#' @export
setMethod("getReference", "gbFeatureList", function(x) getReference(.seqinfo(x)))
#' @export
setMethod("getKeywords", "gbFeatureList", function(x) getKeywords(.seqinfo(x)))
#' @export
setMethod("getComment", "gbFeatureList", function(x) getComment(.seqinfo(x)))

#' @export
#' @rdname getHeader-methods
setMethod("getHeader", "gbFeatureList", function(x) .header(.seqinfo(x)))

#' @export
#' @rdname getHeader-methods
setMethod("header", "gbFeatureList", function(x) .header(.seqinfo(x)))


#' @export
#' @rdname getSequence-methods
setMethod("getSequence", "gbFeatureList", function(x) .seq_access(x))


#' @export
setMethod("ranges", "gbFeatureList", function(x, join = FALSE, key = TRUE,
                                              include = "none", exclude = "") {
  .make_GRanges(x, join = join, include = include, exclude = exclude, key = key)
})

#' @export
#' @rdname start-methods
setMethod("start", "gbFeatureList", function(x, join = FALSE, drop = TRUE) {
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
setReplaceMethod("start", "gbFeatureList", function(x, check=TRUE, value) {
  value <- recycle(value, length(x))
  new_x <- Map(function(Feature, check, val) { 
    start(Feature, check=check) <- val
    Feature
  }, Feature=x, check=list(check), val=value)
  
  new('gbFeatureList', .Data=new_x, .seqinfo=x@.seqinfo)
})

#' @export
#' @rdname end-methods
setMethod("end", "gbFeatureList", function(x, join = FALSE, drop = TRUE) {
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
setReplaceMethod("end", "gbFeatureList", function(x, check=TRUE, value) {
  value <- recycle(value, length(x))
  new_x <- Map(function(Feature, check, val) { 
    end(Feature, check=check) <- val
    Feature
  }, Feature=x, check=list(check), val=value)
  
  new('gbFeatureList', .Data=new_x, .seqinfo=x@.seqinfo)
})

#' @export
#' @rdname strand-methods
setMethod("strand", "gbFeatureList", function(x, join = FALSE) {
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
setReplaceMethod("strand", "gbFeatureList", function(x, value) {
  value <- recycle(value, length(x))
  new_x <- Map(function(Feature, val) { 
    strand(Feature) <- val
    Feature
  }, Feature=x, val=value)
  
  new('gbFeatureList', .Data=new_x, .seqinfo=x@.seqinfo)
})

#' @export
#' @rdname width-methods
setMethod("width", "gbFeatureList", function(x) {
  ans <- lapply(x, width)
  if (all(vapply(ans, length, 0) == 1L)) {
    unlist(ans)
  } else {
    ans
  }
})

#' @export
#' @rdname width-methods
setMethod("joint_width", "gbFeatureList", function(x) {
  unlist(lapply(x, joint_width))
})

#' @export
#' @rdname dbxref-methods
setMethod("dbxref", "gbFeatureList", function(x, db = NULL, na.rm = TRUE, ...) {     
  ans <- lapply(x, dbxref, db=db)
  names(ans) <- lapply(x, function(f) sprintf("%s.%s", f@key, f@.id))
  if (na.rm)
    ans <- compactNA(ans)
  if (all_empty(ans)) {
    return(NA_character_)
  }
  ans
})

#' @export
#' @rdname location-methods
setMethod("location", "gbFeatureList", function(x, join = FALSE) {
  lapply(x, location)
})

#' @export
#' @rdname fuzzy-methods
setMethod("fuzzy", "gbFeatureList", function(x) {
  do.call(rbind, lapply(x, fuzzy))
})

#' @export
#' @rdname index-methods
setMethod("index", "gbFeatureList", function(x) {
  vapply(x, function(x) x@.id, numeric(1))
})

#' @export
#' @rdname key-methods
setMethod("key", "gbFeatureList", function(x) {
  vapply(x, function(f) f@key, character(1))
})

#' @export
#' @rdname qualif-methods
setMethod("qualif", "gbFeatureList", function(x, which = "", fixed = FALSE, use.names = TRUE) {
  ans <- .qual_access(x, which, fixed, use.names)
  if (use.names) {
    .simplify(ans, unlist=FALSE)
  } else {
    .simplify(ans, unlist=TRUE)
  }
})


# listers ----------------------------------------------------------------


#' @export
#' @rdname listQualif-methods
setMethod("listQualif", "gbFeatureList", function(x) {
  lapply(x, listQualif)
})

tbl_qual <- function(x) {
  tblqnm <- Compose(Partial(table, deparse.level=0), "names", 
                    Partial(slot, name="qualifiers"))
  lapply(x, tblqnm)
}

tbl_merge <- function(a, b) {
  n <- intersect(names(a), names(b)) 
  res <- c(a[!(names(a) %in% n)], b[!(names(b) %in% n)], a[n] + b[n])
  res[order(names(res))]
}

#' @export
#' @rdname tableQualif-methods
setMethod("tableQualif", "gbFeatureList", function(x) {
  tbls <- tbl_qual(x)
  Reduce(tbl_merge, tbls)
})


# testers ----------------------------------------------------------------


#' @export
#' @rdname hasKey-methods
setMethod("hasKey", "gbFeatureList", function(x, key) {
  vapply(x, hasKey, key, FUN.VALUE=FALSE)
})

#' @export
#' @rdname hasQualif-methods
setMethod("hasQualif", "gbFeatureList", function(x, qualifier) {
  vapply(x, hasQualif, qualifier, FUN.VALUE=FALSE)
})


# subsetting ----------------------------------------------------------


#' @export
#' @rdname extract-methods
setMethod("[", c("gbFeatureList", "character", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            check <- list(...)$check %||% TRUE
            idx <- which(vapply(x@.Data, function(f) f@key, "") == i)
            IRanges::new2('gbFeatureList', .Data=x@.Data[idx], .seqinfo=x@.seqinfo,
                          check=check)
          })

#' @export
#' @rdname extract-methods
setMethod("[", c("gbFeatureList", "numeric", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            check <- list(...)$check %||% TRUE
            IRanges::new2('gbFeatureList', .Data=x@.Data[i], .seqinfo=x@.seqinfo,
                          check=check)
          })

#' @export
#' @rdname extract-methods
setMethod("[", c("gbFeatureList", "logical", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
            check <- list(...)$check %||% TRUE
            IRanges::new2('gbFeatureList', .Data=x@.Data[i], .seqinfo=x@.seqinfo, 
                          check=check)
          })

#' @export
#' @rdname extract-methods
setMethod("[", c("gbFeatureList", "missing", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) x
          )


#' @export
#' @rdname extract-methods
setMethod("[[", "gbFeatureList",
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


setMethod("view", "gbFeatureList", function(x, n)  {
  for (i in x[seq(if (missing(n)) length(x) else n)]) {
    show(i)
    cat("\n")
  }
})


# select, shift, revcomp ----------------------------------------------------


#' @export
#' @rdname select-methods
setMethod("select", "gbFeatureList", function(x, ..., keys = NULL, cols = NULL) {
  .retrieve(.select(x, ..., keys = keys), cols = cols)
})


#' @export
#' @rdname shift-methods
setMethod("shift", "gbFeatureList", function(x, shift=0L, split=FALSE, order=FALSE) {
  .shift(x=x, shift=shift, split=split, order=order)
})


#' @export
#' @rdname revcomp-methods
setMethod("revcomp", "gbFeatureList", function(x, order=TRUE) {
  .revcomp(x=x, order=order)
})

