
# gbFeatureList-class -------------------------------------------------

#' @include gbFeature-class.r
NULL

setOldClass("list")

#' gbFeatureList class
#' 
#' \dQuote{\code{gbFeatureList}} is an S4 class that extends the
#' \code{\linkS4class{gbRecord}} class. This class provides a container for
#' feature data retrived from GenBank flat files.
#'
#' \dQuote{\code{gbFeatureList}} instances have four slots in addition to the7
#' slots provided by \code{\linkS4class{gbRecord}} objects:
#' 
#' @slot .Dir The path to the database file containing the GenBank
#' record the feature list is part of.
#' @slot .ACCN Accession number of the GenBank record that the
#' feature list is part of.
#' @slot .DEF The definition line (brief description of the sequence)
#' of the GenBank record the feature list is part of.
#' @slot .Data A list of \dQuote{\code{gbFeature}} instances.
#' 
#' @rdname gbFeatureList
#' @exportClass gbFeatureList
#' @classHierarchy
#' @classMethods
.gbFeatureList <- setClass("gbFeatureList", 
                           representation(.Dir="character",
                                          .ACCN="character",
                                          .DEF="character"),
                           contains="list")


setValidity("gbFeatureList", function (object) {
  if (any(vapply(object@.Data, class, character(1)) != "gbFeature")) {
    return("All elements in a 'gbFeatureList' must be 'gbFeature' instances")
  }
  if (any(vapply(object@.Data, function(x) x@.Dir, character(1)) != object@.Dir)) {
    return("All elements in a 'gbFeatureList' must stem from a valid 'gbRecord' instance")
  }
  if (any(vapply(object@.Data, function(x) x@.ACCN, character(1)) != object@.ACCN)) {
    return("All elements in a 'gbFeatureList' must stem from the same 'gbRecord' instance")
  }
  
  TRUE
})


# show -------------------------------------------------------------------


setMethod("show", "gbFeatureList", 
          function (object) {
            n_f <- length(object)
            cat(sprintf("'%s' with %i features:\n\n", 
                        class(object), n_f))
            if (n_f > 0L) {
              show(object[[1L]])
              if (n_f > 1L) {
                cat("\n...\n")
                show(object[[n_f]])
              }
            }
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


setMethod("accession", "gbFeatureList",
          function (x) x@.ACCN)


setMethod("definition", "gbFeatureList",
          function (x) x@.DEF)


setMethod("range", "gbFeatureList",
          function (x, join = FALSE) {
            start <- as.integer(unlist(start(x, join = join)))
            width <- as.integer(unlist(end(x, join = join))) - start + 1L
            strand <- unlist(strand(x, join = join))
            .gbRange(start, width, strand)
          })


setMethod("location", "gbFeatureList",
          function (x, attributes = FALSE, join = FALSE) {
            ans <- range(x, join = join)
            keys <- key(x, attributes=FALSE)
            ids <- index(x, attributes=FALSE)
            if (join || length(ans) == length(keys)) {
              ans@elementMetadata$feature <- keys
              ans@elementMetadata$id <- ids
            } else {
              exp <- expandIds(x)
              ans@elementMetadata$feature <- exp$keys
              ans@elementMetadata$id <- exp$ids
            }
            if (attributes) {
              structure(ans,
                        accession=x@.ACCN,
                        definition=x@.DEF, 
                        database=x@.Dir)
            }
            else {
              ans
            }
          })


setMethod("index", "gbFeatureList",
          function (x, attributes = FALSE) { 
            ans <- vapply(x, function(f) f@.ID, numeric(1))
            if (attributes) {
              ans <- structure(ans, accession=x@.ACCN,
                               definition=x@.DEF,
                               database=x@.Dir)
            }
            ans
          })


setMethod("key", "gbFeatureList",
          function (x, attributes = FALSE) {
            ans <- vapply(x, function(f) f@key, character(1))
            if (attributes) {
              ans <- structure(ans,
                               id=vapply(x, function(f) f@.ID, numeric(1)),
                               accession=x@.ACCN,
                               definition=x@.DEF,
                               database=x@.Dir)
            }
            ans
          })


setMethod("qualif", "gbFeatureList",
          function (x, which, attributes = FALSE, fixed = FALSE) {
            if (missing(which))
              which <- ""
            ans <- .qualAccess(x, which, fixed) %@% .simplify(unlist=FALSE)
            if (attributes) {
              ans <- structure(ans, 
                               id=vapply(x, function(f) f@.ID, numeric(1)),
                               accession=x@.ACCN,
                               definition=x@.DEF,
                               database=x@.Dir)
            }
            
            ans
          })


setMethod("dbxref", "gbFeatureList",
          function (x, db = NULL, na.rm = TRUE, ...) {     
            ans <- lapply(x, dbxref, db=db)
            names(ans) <- lapply(x, function(f) sprintf("%s.%s", f@key, f@.ID))
            if (na.rm)
              ans <- ans[!is.na(ans)]
            if (is_empty(ans)) {
              return(NA_character_)
            }
            
            ans
          })


setMethod("sequence", "gbFeatureList",
          function (x, db = NULL) {
            stopifnot(hasValidDb(x))
            db <- initGB(x@.Dir, verbose=FALSE)
            .seqAccess(dbFetch(db, "sequence"), x, dbFetch(db, "type"))
          })


# setters ----------------------------------------------------------------


setReplaceMethod("start", "gbFeatureList",
                 function (x, value) {
                   value <- recycle(x, value)
                   new_x <- Map(function(Feature, val) { 
                     start(Feature) <- val
                     Feature
                   }, Feature=x, val=value)
                   
                   .gbFeatureList(.Data=new_x, .Dir=x@.Dir,
                                  .ACCN=x@.ACCN, .DEF=x@.DEF)
                 })


setReplaceMethod("end", "gbFeatureList",
                 function(x, value) {
                   value <- recycle(x, value)
                   new_x <- Map(function(Feature, val) { 
                     end(Feature) <- val
                     Feature
                   }, Feature=x, val=value)
                   
                   .gbFeatureList(.Data=new_x, .Dir=x@.Dir,
                                  .ACCN=x@.ACCN, .DEF=x@.DEF)
                 })


setReplaceMethod("strand", "gbFeatureList",
                 function(x, value) {
                   value <- recycle(x, value)
                   new_x <- Map(function(Feature, val) { 
                     strand(Feature) <- val
                     Feature
                   }, Feature=x, val=value)
                   
                   .gbFeatureList(.Data=new_x, .Dir=x@.Dir,
                                  .ACCN=x@.ACCN, .DEF=x@.DEF)
                 })


# testers ----------------------------------------------------------------


setMethod("hasKey", "gbFeatureList", 
          function (x, key)
            vapply(x, hasKey, key, FUN.VALUE=logical(1)))


setMethod("hasQualif", "gbFeatureList", 
          function (x, qualifier)
            vapply(x, hasQualif, qualifier, FUN.VALUE=logical(1)))


# subsetting ----------------------------------------------------------


#' @export
setMethod("[", c("gbFeatureList", "character", "missing", "ANY"),
          function (x, i, j, ..., drop = TRUE) {
            idx <- vapply(x@.Data, function(f) f@key, character(1L)) == i       
            .gbFeatureList(.Data=x@.Data[idx], .Dir=x@.Dir, .ACCN=x@.ACCN,
                           .DEF=x@.DEF)
          })


#' @export
setMethod("[", c("gbFeatureList", "numeric", "missing", "ANY"),
          function (x, i, j, ..., drop = TRUE) {
            .gbFeatureList(.Data=x@.Data[i], .Dir=x@.Dir, .ACCN=x@.ACCN,
                           .DEF=x@.DEF)
          })


#' @export
setMethod("[", c("gbFeatureList", "logical", "missing", "ANY"),
          function (x, i, j, ..., drop = TRUE) {
            .gbFeatureList(.Data=x@.Data[i], .Dir=x@.Dir, .ACCN=x@.ACCN,
                           .DEF=x@.DEF)
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


