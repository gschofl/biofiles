
# gbFeatureList-class -------------------------------------------------

#' @include gbFeature-class.R
NULL

setOldClass("list")

#' gbFeatureList
#' 
#' \dQuote{gbFeatureList} is an S4 class that provides a container for
#' \dQuote{\linkS4class{gbFeature}}s retrived from GenBank flat files.
#'
#' @slot .Info A \code{\linkS4class{gbInfo}} instance.
#' @slot .Data A list of \code{\linkS4class{gbFeature}} objects.
#' 
#' @rdname gbFeatureList
#' @export
#' @classHierarchy
#' @classMethods
setClass("gbFeatureList", 
         representation(.Info="gbInfo"),
         contains="list")


setValidity2("gbFeatureList", function (object) {
  if (any(vapply(object@.Data, class, character(1)) != "gbFeature")) {
    return("All elements in a 'gbFeatureList' must be 'gbFeature' instances")
  }
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
          function (object, ...) {
            x <- lapply(object, summary)
            return(invisible(x))
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
          function (x) x@.Info)


setMethod("seqlengths", "gbFeatureList",
          function (x) seqlengths(seqinfo(x)))


setMethod("accession", "gbFeatureList",
          function (x) seqnames(seqinfo(x)))


setMethod("definition", "gbFeatureList",
          function (x) genome(seqinfo(x)))


setMethod("ranges", "gbFeatureList",
          function (x, join = TRUE, with_qual = "none", without_qual = "") {
            .make_GRanges(x, join = join, with_qual = with_qual,
                          without_qual = without_qual)
          })



setMethod("location", "gbFeatureList",
          function (x, seqinfo = FALSE, join = FALSE) {
            ans <- lapply(x, location)
            if (seqinfo) {
              structure(ans,
                        accession=accession(x),
                        definition=unname(definition(x)),
                        dir=seqinfo(x)@db@dir)
            }
            else {
              ans
            }
          })


setMethod("index", "gbFeatureList",
          function (x, seqinfo = FALSE) { 
            ans <- vapply(x, function(f) f@.Id, numeric(1))
            if (seqinfo) {
              structure(ans,
                        accession=accession(x),
                        definition=unname(definition(x)),
                        dir=seqinfo(x)@db@dir)
            } else {
              ans
            }
          })


setMethod("key", "gbFeatureList",
          function (x, seqinfo = FALSE) {
            ans <- vapply(x, function(f) f@key, character(1))
            if (seqinfo) {
              structure(ans,
                        accession=accession(x),
                        definition=unname(definition(x)),
                        dir=seqinfo(x)@db@dir)
            } else {
              ans  
            }
          })


setMethod("qualif", "gbFeatureList",
          function (x, which = "", seqinfo = FALSE, fixed = FALSE) {
            ans <- .qualAccess(x, which, fixed) %@% .simplify(unlist=FALSE)
            if (seqinfo) {
              ans <- structure(ans, 
                               id=vapply(x, function(f) f@.Id, numeric(1)),
                               accession=accession(x),
                               definition=unname(definition(x)),
                               dir=seqinfo(x)@db@dir)
            } else {
              ans            
            }
          })


setMethod("dbxref", "gbFeatureList",
          function (x, db = NULL, na.rm = TRUE, ...) {     
            ans <- lapply(x, dbxref, db=db)
            names(ans) <- lapply(x, function(f) sprintf("%s.%s", f@key, f@.Id))
            if (na.rm)
              ans <- compactNA(ans)
            if (all_empty(ans)) {
              return(NA_character_)
            }
            
            ans
          })


setMethod("sequence", "gbFeatureList",
          function (x, db = NULL) {
            stopifnot(hasValidDb(x))
            db <- slot(seqinfo(x), "db")
            .seqAccess(dbFetch(db, "sequence"), x, dbFetch(db, "type"))
          })


# setters ----------------------------------------------------------------


setReplaceMethod("start", "gbFeatureList",
                 function (x, value) {
                   value <- recycle(value, length(x))
                   new_x <- Map(function(Feature, val) { 
                     start(Feature) <- val
                     Feature
                   }, Feature=x, val=value)
                   
                   new('gbFeatureList', .Data=new_x, .Info=seqinfo(x))
                 })


setReplaceMethod("end", "gbFeatureList",
                 function(x, value) {
                   value <- recycle(value, length(x))
                   new_x <- Map(function(Feature, val) { 
                     end(Feature) <- val
                     Feature
                   }, Feature=x, val=value)
                   
                   new('gbFeatureList', .Data=new_x, .Info=seqinfo(x))
                 })


setReplaceMethod("strand", "gbFeatureList",
                 function(x, value) {
                   value <- recycle(value, length(x))
                   new_x <- Map(function(Feature, val) { 
                     strand(Feature) <- val
                     Feature
                   }, Feature=x, val=value)
                   
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
          function(x, shift=0L, split=FALSE, order=FALSE, updateDb=FALSE) {
            .shift_features(x=x, shift=shift, split=split,
                            order=order, updateDb=updateDb)
          })


# revcomp ----------------------------------------------------------------


setMethod("revcomp", "gbFeatureList",
          function(x, order=FALSE, updateDb=FALSE) {
            .revcomp_features(x=x, order=order, updateDb=updateDb)
          })


