
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
  if (names(seq) != accession(object))
    return("Names of 'seqinfo' and 'sequence' do not match")
  if (seq@ranges@width != unname(seqlengths(x)))
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
            if (length(object) > 2*n) {
              head <- head(object, n=n)
              tail <- tail(object, n=n)
              x <- lapply(head, summary)
              cat("...\n")
              x <- lapply(tail, summary)  
            } else  {
              x <- lapply(object, summary)
            }

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


setMethod("seqlengths", "gbFeatureList",
          function (x) seqlengths(seqinfo(x)))


setMethod("accession", "gbFeatureList",
          function (x) seqnames(seqinfo(x)))


setMethod("definition", "gbFeatureList",
          function (x) genome(seqinfo(x)))


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
            vapply(x, function(f) f@.id, numeric(1))
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
            names(ans) <- lapply(x, function(f) sprintf("%s.%s", f@key, f@.Id))
            if (na.rm)
              ans <- compact(ans, filter=function(x) all(is.na(x)))
            if (all_empty(ans)) {
              return(NA_character_)
            }
            
            ans
          })


setMethod("sequence", "gbFeatureList",
          function (x) .seqAccess(x))


# setters ----------------------------------------------------------------


setReplaceMethod("start", "gbFeatureList",
                 function (x, check=TRUE, value) {
                   value <- recycle(value, length(x))
                   new_x <- Map(function(Feature, check, val) { 
                     start(Feature, check=check) <- val
                     Feature
                   }, Feature=x, check=list(check), val=value)
                   
                   new('gbFeatureList', .Data=new_x, .Info=seqinfo(x))
                 })


setReplaceMethod("end", "gbFeatureList",
                 function(x, check=TRUE, value) {
                   value <- recycle(value, length(x))
                   new_x <- Map(function(Feature, check, val) { 
                     end(Feature, check=check) <- val
                     Feature
                   }, Feature=x, check=list(check), val=value)
                   
                   new('gbFeatureList', .Data=new_x, .Info=seqinfo(x))
                 })


setReplaceMethod("strand", "gbFeatureList",
                 function(x, check=TRUE, value) {
                   value <- recycle(value, length(x))
                   new_x <- Map(function(Feature, check, val) { 
                     strand(Feature, check=check) <- val
                     Feature
                   }, Feature=x, check=list(check), val=value)
                   
                   new('gbFeatureList', .Data=new_x, .Info=seqinfo(x))
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
            idx <- which(vapply(x@.Data, function(f) f@key, character(1L)) == i)
            new('gbFeatureList', .Data=x@.Data[idx], .Info=seqinfo(x))
          })

#' @export
setMethod("[", c("gbFeatureList", "numeric", "missing", "ANY"),
          function (x, i, j, ..., drop = TRUE) {
            new('gbFeatureList', .Data=x@.Data[i], .Info=seqinfo(x))
          })


#' @export
setMethod("[", c("gbFeatureList", "logical", "missing", "ANY"),
          function (x, i, j, ..., drop = TRUE) {
            new('gbFeatureList', .Data=x@.Data[i], .Info=seqinfo(x))
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

