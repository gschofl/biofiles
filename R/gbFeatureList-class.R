
# gbFeatureList-class -------------------------------------------------

#' @include gbFeature-class.R
NULL

setOldClass("list")

#' gbFeatureList
#' 
#' \dQuote{gbFeatureList} is an S4 class that provides a container for
#' \dQuote{\linkS4class{gbFeature}}s retrived from GenBank flat files.
#'
#' @slot .seqinfo An \code{environment} containing the genome sequence as
#' an \code{\linkS4class{XStringSet}} object and sequence metadata
#' as a \code{\linkS4class{Seqinfo}} object.
#' @slot .Data A list of \code{\linkS4class{gbFeature}} objects.
#' 
#' @rdname gbFeatureList
#' @export
#' @classHierarchy
#' @classMethods
setClass("gbFeatureList", 
         representation(.seqinfo="environment"),
         prototype(.seqinfo=new.env(parent=emptyenv())),
         contains="list")


setValidity2("gbFeatureList", function (object) {
  if (!all(vapply(object@.Data, is, 'gbFeature', FUN.VALUE=logical(1))))
    return("All elements in a 'gbFeatureList' must be 'gbFeature' instances")
  seq <- get("sequence", object@.seqinfo)
  if (names(seq) != getAccession(object))
    return("Names of 'seqinfo' and 'sequence' do not match")
  if (seq@ranges@width != unname(getLength(object)))
    return("Length of 'seqinfo' and 'sequence' do not match")
                 
  TRUE
})


# show -------------------------------------------------------------------


setMethod("show", "gbFeatureList", 
          function (object) {
            lo <- length(object)
            cat(sprintf("%s with %i features:\n",  sQuote(class(object)), lo))
            if (lo > 0L) {
              .showGbFeature(object[[1L]], showInfo=FALSE)
              if (lo > 1L) {
                cat("...\n")
                .showGbFeature(object[[lo]], showInfo=FALSE)
              }
            }
            cat("Seqinfo:\n")
            showInfo(seqinfo(object))
            return(invisible(object))
          })


# summary ----------------------------------------------------------------


setMethod("summary", "gbFeatureList",
          function (object, n=8, ...) {
            olen <- length(object)
            if (olen > 2*n) {
              hd <- object[seq_len(n), check=FALSE]
              tl <- object[seq.int(to = olen, length.out = min(n, olen)), check=FALSE]
              idx <- c("N", index(hd), "...", index(tl))
              key <- c("Key", key(hd), "...", key(tl))
              loc <- c(location(hd), "...", location(tl))
              loc <- c("Location", unlist(lapply(loc, as, "character")))
              prod <- c("Product", unlist(product(hd), use.names=FALSE),
                        "...", unlist(product(tl), use.names=FALSE))
            } else {
              idx <- c("N", index(object))
              key <- c("Key", key(object))
              loc <- c("Location", unlist(lapply(location(object), as, "character")))
              prod <- c("Product", unlist(product(object), use.names=FALSE))
            }
            setoff <- rep(dup(" ", list(...)$setoff %||% 0), length(idx))
            max_idx_len <- max(nchar(idx))
            max_key_len <- max(nchar(key))
            max_loc_len <- max(nchar(loc))
            idx <- pad(idx, max_idx_len + 1, "right")
            key <- pad(key, max_key_len + 1, "right")
            loc <- pad(loc, max_loc_len + 1, "right")
            showme <- ellipsize(sprintf("%s%s%s%s%s", setoff, idx, key, loc, prod),
                                width=getOption("width") - 1)
            cat(showme, sep="\n")
            return(invisible(NULL))
          })


# getters ----------------------------------------------------------------


setMethod("start", "gbFeatureList",
          function (x, join = FALSE, drop = TRUE) {
            ans <- lapply(x, start, join = join, drop = drop)
            if (drop) {
              if (join || all(vapply(ans, length, numeric(1)) == 1L)) {
                unlist(ans)
              } else {
                ans
              }
            } else {
              setNames(data.frame(do.call(rbind, ans)), "start")
            }
          })


setMethod("end", "gbFeatureList",
          function (x, join = FALSE, drop = TRUE) {
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


setMethod("strand", "gbFeatureList",
          function (x, join = FALSE) {
            ans <- lapply(x, strand, join = join)        
            if (join || all(vapply(ans, length, numeric(1)) == 1L)) {
              unlist(ans)
            } else {
              ans
            }
          })


setMethod("width", "gbFeatureList",
          function (x, join = FALSE) {
            ans <- lapply(x, width, join = join)
            if (join || all(vapply(ans, length, numeric(1)) == 1L)) {
              unlist(ans)
            } else {
              ans
            }
          })


setMethod("seqinfo", "gbFeatureList",
          function (x) tryCatch(get("seqinfo", x@.seqinfo),
                                error = function (e) Seqinfo() ))


setMethod("getLength", "gbFeatureList",
          function (x) unname(seqlengths(seqinfo(x))))


setMethod("getAccession", "gbFeatureList",
          function (x) seqnames(seqinfo(x)))


setMethod("getDefinition", "gbFeatureList",
          function (x) unname(genome(seqinfo(x))))


setMethod("ranges", "gbFeatureList",
          function (x, join = FALSE, key = TRUE, include = "none", exclude = "") {
            .make_GRanges(x, join = join, include = include, exclude = exclude, key = key)
          })


setMethod("location", "gbFeatureList",
          function (x, join = FALSE) {
            lapply(x, location)
          })


setMethod("index", "gbFeatureList",
          function (x) {
            vapply(x, function(x) x@.id, numeric(1))
          })


setMethod("key", "gbFeatureList",
          function (x) {
            vapply(x, function(f) f@key, character(1))
          })


setMethod("qualif", "gbFeatureList",
          function (x, which = "", fixed = FALSE) {
            .simplify(.qualAccess(x, which, fixed), unlist=FALSE)
          })


setMethod("dbxref", "gbFeatureList",
          function (x, db = NULL, na.rm = TRUE, ...) {     
            ans <- lapply(x, dbxref, db=db)
            names(ans) <- lapply(x, function(f) sprintf("%s.%s", f@key, f@.id))
            if (na.rm)
              ans <- compact(ans, filter=function(x) all(is.na(x)))
            if (all_empty(ans)) {
              return(NA_character_)
            }
            
            ans
          })


setMethod("getSequence", "gbFeatureList", function (x) .seqAccess(x))


setMethod('.dbSource', 'gbFeatureList', function (x) {
  parse_dbsource(get("dbsource", x@.seqinfo))
})


setMethod(".defline", "gbFeatureList", function (x) {
  vapply(x, .defline, character(1), USE.NAMES=FALSE)
})


# setters ----------------------------------------------------------------


setReplaceMethod("start", "gbFeatureList",
                 function (x, check=TRUE, value) {
                   value <- recycle(value, length(x))
                   new_x <- Map(function(Feature, check, val) { 
                     start(Feature, check=check) <- val
                     Feature
                   }, Feature=x, check=list(check), val=value)
                   
                   new('gbFeatureList', .Data=new_x, .seqinfo=x@.seqinfo)
                 })


setReplaceMethod("end", "gbFeatureList",
                 function(x, check=TRUE, value) {
                   value <- recycle(value, length(x))
                   new_x <- Map(function(Feature, check, val) { 
                     end(Feature, check=check) <- val
                     Feature
                   }, Feature=x, check=list(check), val=value)
                   
                   new('gbFeatureList', .Data=new_x, .seqinfo=x@.seqinfo)
                 })


setReplaceMethod("strand", "gbFeatureList",
                 function(x, check=TRUE, value) {
                   value <- recycle(value, length(x))
                   new_x <- Map(function(Feature, check, val) { 
                     strand(Feature, check=check) <- val
                     Feature
                   }, Feature=x, check=list(check), val=value)
                   
                   new('gbFeatureList', .Data=new_x, .seqinfo=x@.seqinfo)
                 })


# listers ----------------------------------------------------------------


setMethod("listQualif", "gbFeatureList", 
          function (x) {
            lapply(x, listQualif)
          })


# testers ----------------------------------------------------------------


setMethod("hasKey", "gbFeatureList", 
          function (x, key) {
            vapply(x, hasKey, key, FUN.VALUE=logical(1))
          })


setMethod("hasQualif", "gbFeatureList", 
          function (x, qualifier) {
            vapply(x, hasQualif, qualifier, FUN.VALUE=logical(1))
          })


# subsetting ----------------------------------------------------------


#' @export
setMethod("[", c("gbFeatureList", "character", "missing", "ANY"),
          function (x, i, j, ..., drop = TRUE) {
            check <- list(...)$check %||% TRUE
            idx <- which(vapply(x@.Data, function(f) f@key, character(1L)) == i)
            IRanges::new2('gbFeatureList', .Data=x@.Data[idx], .seqinfo=x@.seqinfo,
                          check=check)
          })

#' @export
setMethod("[", c("gbFeatureList", "numeric", "missing", "ANY"),
          function (x, i, j, ..., drop = TRUE) {
            check <- list(...)$check %||% TRUE
            IRanges::new2('gbFeatureList', .Data=x@.Data[i], .seqinfo=x@.seqinfo,
                          check=check)
          })


#' @export
setMethod("[", c("gbFeatureList", "logical", "missing", "ANY"),
          function (x, i, j, ..., drop = TRUE) {
            check <- list(...)$check %||% TRUE
            IRanges::new2('gbFeatureList', .Data=x@.Data[i], .seqinfo=x@.seqinfo, 
                          check=check)
          })


#' @export
setMethod("[", c("gbFeatureList", "missing", "missing", "ANY"),
          function (x, i, j, ..., drop = TRUE) x
          )


# select -----------------------------------------------------------------


setMethod("select", "gbFeatureList", 
          function (x, ..., keys = NULL, cols = NULL) {
            ans <- .select(x, ..., keys = keys)
            ans <- .retrieve(ans, cols = cols)
            ans
          })

    
# view -------------------------------------------------------------------


setMethod("view", "gbFeatureList", 
          function (x, n)  {
            for (i in x[seq(if (missing(n)) length(x) else n)]){
              show(i)
              cat("\n")
            }
          })


# shift ------------------------------------------------------------------


setMethod("shift", "gbFeatureList",
          function(x, shift=0L, split=FALSE, order=FALSE) {
            .shift_features(x=x, shift=shift, split=split, order=order)
          })

