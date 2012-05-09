
# gbFeatureList-class -------------------------------------------------

##' @include gbFeature-class.r
NULL

setOldClass("list")

##' gbFeatureList class
##' 
##' gbFeatureList is an S4 class that extends the
##' \code{\link{gbRecord-class}}. This class provides a container for
##' feature data retrived from GenBank flat files.
##'
##' gbFeatureList objects have four slots in addition to the slots provided
##' by \code{\link{gbRecord-class}} objects:
##' \describe{
##'    \item{.Dir}{The path to the database file containing the GenBank
##'    record the feature list is part of.}
##'    \item{.ACCN}{Accession number of the GenBank record that the
##'    feature list is part of.}
##'    \item{.DEF}{The definition line (brief description of the sequence)
##'    of the GenBank record the feature list is part of.}
##'    \item{.Data}{A list of gbFeature objects}
##' }
##' 
##' @importFrom Biostrings xscat
##' @importFrom Biostrings reverseComplement
##' 
##' @name gbFeatureList-class
##' @rdname gbFeatureList-class
##' @exportClass gbFeatureList
##' @aliases show,gbFeatureList-method
##' @aliases start,gbFeatureList-method
##' @aliases end,gbFeatureList-method
##' @aliases strand,gbFeatureList-method
##' @aliases start<-,gbFeatureList-method
##' @aliases end<-,gbFeatureList-method
##' @aliases strand<-,gbFeatureList-method
##' @aliases width,gbFeatureList-method
##' @aliases range,gbFeatureList-method
##' @aliases getIndex,gbFeatureList-method
##' @aliases getKey,gbFeatureList-method
##' @aliases getLocation,gbFeatureList-method
##' @aliases getQualifier,gbFeatureList-method
##' @aliases dbXref,gbFeatureList-method
##' @aliases getSequence,gbFeatureList-method
##' @aliases hasKey,gbFeatureList-method
##' @aliases hasQualifier,gbFeatureList-method
##' @aliases [,gbFeatureList-method
##' @aliases select,select-method,gbFeatureList-method
##' @aliases view,view-method,gbFeatureList-method
.gbFeatureList <- 
  #### gbFeatureList ####
  setClass("gbFeatureList", 
           representation(.Dir="character",
                          .ACCN="character",
                          .DEF="character"),
           contains="list")


# show-method ---------------------------------------------------------


##' @export
setMethod("show",
          signature="gbFeatureList", 
          function (object) {
            n_f <- length(object)
            cat(sprintf("'%s' with %i features:\n\n", 
                        class(object), n_f))
            if (n_f > 0L) {
              biofiles::show(object[[1L]])
              if (n_f > 1L) {
                cat("\n...\n")
                biofiles::show(object[[n_f]])
              }
            }
            return(invisible(object))
          })


# Constructor ---------------------------------------------------------


#' @keywords internal
gbFeatureList <- function(db_dir, accession, definition, features)
{
  if (!is.list(features))
    stop("'features' must be a list")
  if (!all(vapply(features, is, "gbFeature", FUN.VALUE=TRUE)))
    stop("all elements in 'features' must be gbFeature objects")
  if (!all(vapply(features, function(f) f@.Dir, FUN.VALUE="") == db_dir))
    stop("all elements in 'features' must be from a valid gbData object")
  if (!all(vapply(features, function(f) f@.ACCN, FUN.VALUE="") == accession))
    stop("all elements in 'features' must be from the same gbData object")
  
  .gbFeatureList(.Data=features, .Dir=as.character(db_dir),
                 .ACCN=as.character(accession), .DEF=as.character(definition)) 
}


# Accessor-generics ---------------------------------------------------


##' Select elements from a GenBank Record
##' 
##' @usage select(x, subset = c(""), select =c(""))
##' 
##' @param x A \sQuote{\code{gbRecord}} or \sQuote{\code{gbFeatureList}}
##' object
##' @param subset Which elements to select from indicated by index, key,
##' location, or qualifier value.
##' @param select Which information to be retrieved from the subsetted
##' elements.
##' 
##' @export
##' @docType methods
setGeneric( "select", function(x, subset = "", select = "") standardGeneric("select") )


# Accessor-methods ----------------------------------------------------


#' @keywords internal
.simplify <- function (x, unlist=TRUE)
{
  if (length(len <- unique(unlist(lapply(x, length)))) > 1L) {
    return(x)
  }
  if (len == 1L && unlist) {
    unlist(x, recursive=FALSE)
  } else if (len >= 1L) {
    n <- length(x)
    r <- as.vector(unlist(x, recursive=FALSE))
    if (prod(d <- c(len, n)) == length(r)) {
      return( data.frame(stringsAsFactors=FALSE,
                         matrix(r, nrow=n, byrow=TRUE,
                                dimnames=if (!(is.null(nm <- names(x[[1L]]))))
                                  list(NULL, nm))) )
    } else {
      x
    }
  }
  else {
    x
  }
}

##' @export
setMethod("start", "gbFeatureList",
          function(x, drop=TRUE) {
            if (drop) {
              if (max(unlist(lapply(ans <- lapply(x, start), length))) > 1L) {
                ans
              } else {
                unlist(ans)
              }
            } else {
              ans <- data.frame(do.call(rbind, lapply(x, start, drop=FALSE)))
              names(ans) <- "start"
              ans   
            }
          })


##' @export
setMethod("start<-", "gbFeatureList",
          function(x, value){
            if (length(value) < length(x)) {
              value <- c(rep(value, length(x)%/%length(value)),
                         value[seq_len(length(x)%%length(value))])
            }
            new_x <- Map(function(Feature, val) { 
              start(Feature) <- val
              Feature }, Feature=x, val=value)
            
            .gbFeatureList(.Data=new_x, .Dir=x@.Dir,
                           .ACCN=x@.ACCN, .DEF=x@.DEF)
          })


##' @export
setMethod("end", "gbFeatureList",
          function(x, drop=TRUE)  {
            if (drop) {
              if (max(unlist(lapply(ans <- lapply(x, end), length))) > 1L) {
                ans
              } else {
                unlist(ans)
              }
            } else {
              ans <- data.frame(do.call(rbind, lapply(x, end, drop=FALSE)))
              names(ans) <- "end"
              ans   
            }
          })

##' @export
setMethod("end<-", "gbFeatureList",
          function(x, value){
            if (length(value) != length(x)) {
              value <- c(rep(value, length(x)%/%length(value)),
                         value[seq_len(length(x)%%length(value))])
            }
            new_x <- Map(function(Feature, val) { 
              end(Feature) <- val
              Feature }, Feature=x, val=value)
            
            .gbFeatureList(.Data=new_x, .Dir=x@.Dir,
                           .ACCN=x@.ACCN, .DEF=x@.DEF)
          })

##' @export
setMethod("strand", "gbFeatureList",
          function(x, drop=TRUE)  {
            if (drop) {
              if (max(unlist(lapply(ans <- lapply(x, strand), length))) > 1L) {
                ans
              } else {
                unlist(ans)
              }
            } else {
              ans <- data.frame(do.call(rbind, lapply(x, strand)))
              names(ans) <- "strand"
              ans   
            }
          })

##' @export
setMethod("strand<-", "gbFeatureList",
          function(x, value){
            if (length(value) != length(x)) {
              value <- c(rep(value, length(x)%/%length(value)),
                         value[seq_len(length(x)%%length(value))])
            }
            new_x <- Map(function(Feature, val) { 
              strand(Feature) <- val
              Feature }, Feature=x, val=value)
            
            .gbFeatureList(.Data=new_x, .Dir=x@.Dir,
                           .ACCN=x@.ACCN, .DEF=x@.DEF)
          })

##' @export
setMethod("width", "gbFeatureList",
          function(x, drop=TRUE)  {
            if (drop) {
              if (max(unlist(lapply(ans <- lapply(x, width), length))) > 1L) {
                ans
              } else {
                unlist(ans)
              }
            } else {
              ans <- data.frame(do.call(rbind, lapply(x, width)))
              names(ans) <- "width"
              ans   
            }
          } )

##' @export
setMethod("range", "gbFeatureList",
          function(x) do.call( rbind, lapply(x, range) ))

##' @export
setMethod("getLocation", 
          #### getLocation-method ####
          signature(x="gbFeatureList"),
          function (x, attributes=TRUE, simplify=TRUE) {
            
            ans <- lapply(x, function(f) 
              cbind(index=f@.ID, key=f@key, range(f@location)))
            
            if (simplify) {
              ans <- do.call(rbind, ans)
              if (attributes) {
                return( structure(ans, accession=x@.ACCN,
                                  definition=x@.DEF,
                                  database=x@.Dir) )
              }
              else {
                return( ans )
              }
            } 
            else {
              if (attributes) {
                return( structure(ans, accession=x@.ACCN,
                                  definition=x@.DEF, 
                                  database=x@.Dir) )
              }
              else {
                return( ans )
              }
            } 
          })

##' @export
setMethod("getIndex",
          #### getIndex-method ####
          signature(x="gbFeatureList"),
          function (x, attributes=TRUE, simplify=TRUE) {
            
            ans <- lapply(x, function(f) f@.ID)
            
            if (simplify) {
              ans <- unlist(ans)
            }
            
            if (attributes) {
              ans <- structure(ans, accession=x@.ACCN,
                               definition=x@.DEF,
                               database=x@.Dir)
            }
            return( ans )
          })

##' @export
setMethod("getKey",
          #### getKey-method ####
          signature(x="gbFeatureList"),
          function (x, attributes=TRUE, simplify=TRUE) {
            
            ans <- lapply(x, function(f) f@key)
            
            if (simplify) {
              ans <- unlist(ans)
            }
            
            if (attributes) {
              ans <- structure(ans,
                               id=vapply(x, function(f) f@.ID, numeric(1)),
                               accession=x@.ACCN,
                               definition=x@.DEF,
                               database=x@.Dir)
            }
            
            return( ans )
          })


##' @export
setMethod("getQualifier",
          #### getQualifier-method ####
          signature(x="gbFeatureList"),
          function (x, which="", attributes=TRUE,
                    simplify=TRUE, fixed=FALSE) {
            
            ans <- lapply(x, getQualifier, which=which,
                          attributes=FALSE, fixed=fixed)
            
            if (simplify) {
              ans <- .simplify(x=ans, unlist=FALSE)
            }
            
            if (is.data.frame(ans) && attributes) {
              ans <- structure(
                data.frame(stringsAsFactors=FALSE,
                           cbind(getIndex(x, attributes=FALSE),
                               getKey(x, attributes=FALSE),
                                 structure(ans, names=NULL))),
                names=c("index", "key", names(ans)),
                accession=x@.ACCN,
                definition=x@.DEF,
                database=x@.Dir)
              
            } else if (is.list(ans) && attributes) {
              ans <- structure(
                ans,
                names=paste0(getKey(x, attributes=FALSE), ".", getIndex(x, attributes=FALSE)),
                accession=x@.ACCN,
                definition=x@.DEF,
                database=x@.Dir)
            }
            
            return( ans )
          })


##' @export
setMethod("dbXref",
          signature(x="gbFeatureList"),
          #### dbXref-method ####
          function (x, db=NULL, na.rm=TRUE, simplify=TRUE, ...) {     
            ans <- lapply(x, dbXref, db=db)
            names(ans) <- lapply(x, function(f) paste0(f@key, ".", f@.ID))
            
            if (simplify)
              ans <- unlist(ans)
            
            if (na.rm)
              ans <- ans[!is.na(ans)]
            
            if (length(ans) == 0) {
              return( NA_character_ )
            }

            return( ans )
          })

##' @export
setMethod("getSequence",
          #### getSequence-method ####
          signature("gbFeatureList"),
          function (x, db=NULL, na.) {
            stopifnot(hasValidDb(x))
            db <- initGB(x@.Dir, verbose=FALSE)
            .seq_access(s=dbFetch(db, "sequence"), x, type=dbFetch(db, "type"))
          })
          
##' @export
setMethod("hasKey",
          #### hasKey-method ####
          signature("gbFeatureList"), 
          function (x, key) {
            vapply(x, hasKey, key, FUN.VALUE=logical(1))
          })

##' @export
setMethod("hasQualifier",
          #### hasQualifier-method ####
          signature("gbFeatureList"), 
          function (x, qualifier) {
            vapply(x, hasQualifier, qualifier, FUN.VALUE=logical(1))
          })


# Subsetting ----------------------------------------------------------


##' @export
setMethod("[",
          signature(x = "gbFeatureList", i = "character", j = "missing", drop = "missing"),
          function (x, i, j, ..., drop = TRUE) {
            idx <- vapply(x@.Data, function(f) f@key, character(1L)) == i
            gbFeatureList(db_dir=x@.Dir, accession=x@.ACCN,
                          definition=x@.DEF, features=x@.Data[idx])
          })

##' @export
setMethod("[",
          signature(x = "gbFeatureList", i = "numeric", j = "missing", drop = "missing"),
          function (x, i, j, ..., drop = TRUE) {
            gbFeatureList(db_dir=x@.Dir, accession=x@.ACCN,
                          definition=x@.DEF, features=x@.Data[i])
          })

##' @export
setMethod("[",
          signature(x = "gbFeatureList", i = "logical", j = "missing", drop = "missing"),
          function (x, i, j, ..., drop = TRUE) {
            gbFeatureList(db_dir=x@.Dir, accession=x@.ACCN,
                          definition=x@.DEF, features=x@.Data[i])
          })

##' @export
setMethod("[",
          signature(x = "gbFeatureList", i = "missing", j = "missing", drop = "missing"),
          function (x, i, j, ..., drop = TRUE) {
            return(x)
          })

##' @export
setMethod("select",
          #### select-method ####
          signature(x="gbFeatureList"), 
          function (x, subset = "", select = "") {
            ans <- .select(x=x, which=subset)
            ans <- .retrieve(x=ans, which=select)
            ans
          })

    
# View ----------------------------------------------------------------


##' View method
##' 
##' @usage view(x, n, ...)
##' 
##' @param x x
##' @param n n
##' @param ... Additional parameters.
##' 
##' @export
##' @docType methods
##' @rdname view-method
setGeneric( "view", function(x, n, ...)
             standardGeneric("view") )

##' @export
setMethod("view",
          #### view-method ####
          signature(x="gbFeatureList"), 
          function (x, n)  {
            for (i in x[seq(if (missing(n)) length(x) else n)]){
              show(i)
              cat("\n")
            }
          })


# shift features ------------------------------------------------------


.shift_features <- function (x,
                             shift=0L,
                             split=FALSE,
                             order=FALSE,
                             update_db=FALSE) {
  
  if (is(x, "gbRecord")) {
    len <- x$length
    features <-x$features
  } else if (is(x, "gbFeatureList")) {
    if (!any(hasKey(x, "source"))) {
      stop("No source key in this gbFeatureList")
    }
    len <- end(x["source"])
    features <- x
  }
  
  update_split <- function(x, split_matrix) {
    if (!is.na(x@location@compound)) {
      stop("Cannot split a compound location")
    }
    x@location@.Data <- split_matrix
    x@location@compound <- "join"
    x@location@partial <- matrix(c(FALSE, TRUE, TRUE, FALSE), ncol=2)
    x@location@accession <- rep(x@location@accession, 2)
    x@location@remote <- rep(x@location@remote, 2)
    x@location@closed <- matrix(rep(x@location@closed, 2), ncol=2)
    x
  }
  
  src <- features[1]
  f <- features[-1]
  
  start_pos <- start(f)
  end_pos <- end(f)
  
  new_start <- Map("+", start_pos, shift)
  new_end <- Map("+", end_pos, shift)
  
  exceeds_len_start <- which(mapply(function (x) any(x > len), new_start) |
                             mapply(function (x) any(x < 0L), new_start)) 
  exceeds_len_end <-  which(mapply(function (x) any(x > len), new_end) |
                            mapply(function (x) any(x < 0L), new_end))
  
  if (length(exceeds_len_start) > 0L || length(exceeds_len_end) > 0L) {
    
    start_end <- intersect(exceeds_len_start, exceeds_len_end)
    
    get_len <- function (x, len) ifelse(x > len, x - len, ifelse(x < 0L, len + x, x))
    new_start[start_end] <- Map(get_len, new_start[start_end], len)
    new_end[start_end] <- Map(get_len, new_end[start_end], len)
    
    end_only <- setdiff(exceeds_len_end, exceeds_len_start)
    
    if (length(end_only) > 0L) {
      if (split) {
        ss <- mapply("-", new_start[end_only], len)
        se <- mapply("-", new_end[end_only], len)
        ## Split Matrix
        sm <- Map(function(ss, se) matrix(c(len + ss, 1, len, se), ncol=2), ss, se)
        f[end_only] <- Map(update_split, x=f[end_only], split_matrix=sm)
        new_start[end_only] <- Map(function(x) x[,1], sm)
        new_end[end_only] <- Map(function(x) x[,2], sm)
      } else {
        stop("This shiftwidth would split feature(s) ", paste(end_only, collapse=", "))
      }
    }
    
  }
  
  start(f) <- new_start
  end(f) <- new_end
  
  if (order) {
    f <- f[order(mapply("[", new_start, 1))]
  }
  
  f <- .gbFeatureList(.Data=c(src, f), .Dir=src@.Dir,
                      .ACCN=src@.ACCN, .DEF=src@.DEF)
  
  if (update_db) {
    db <- initGB(src@.Dir, verbose=FALSE)
    dbInsert(db, key="features", value=f)
    
    seq <- dbFetch(db, "sequence")
    
    if (shift > 0) {
      shift_point <- seq@ranges@width - shift + 1L
    } else {
      shift_point <- 0L - shift + 1L
    }
      
    new_seq <- xscat(subseq(seq, start = shift_point), 
                     subseq(seq, start = 1L, end = shift_point - 1L))
    names(new_seq) <- names(seq)
    dbInsert(db, key="sequence", value=new_seq)
  }
  
  return( f )
}


.revcomp_features <- function (x, order=FALSE, update_db=FALSE) {
  
  if (is(x, "gbRecord")) {
    max_len <- x$length
    f <-x$features
  } else if (is(x, "gbFeatureList")) {
    
    if (length(x["source"]) == 0) {
      stop("No source field in gbFeatureList")
    }
    
    max_len <- end(x["source"])
    f <- x
  }
  
  new_end <- max_len - start(f) + 1
  new_start <- max_len - end(f) + 1
  new_strand <- strand(f)*-1
  
  start(f) <- new_start
  end(f) <- new_end
  strand(f) <- new_strand
  
  if (order) {
    f <- f[order(new_start)]
  }
  
  f <- .gbFeatureList(.Data=f, .Dir=f@.Dir, .ACCN=f@.ACCN, .DEF=f@.DEF)
  
  if (update_db) {
    db <- initGB(f@.Dir)
    dbInsert(db, key="features", value=f)
    
    seq <- dbFetch(db, "sequence")
    new_seq <- reverseComplement(seq)
    dbInsert(db, key="sequence", value=new_seq)
  }
  
  return( f )
}


#' @export
setMethod("shift", "gbFeatureList",
          function(x, shift=0L, split=FALSE, order=FALSE, update_db=FALSE) {
            .shift_features(x=x, shift=shift, split=split, order=order,
                            update_db=update_db)
          })


##' reverse complement features in a GenBank record
##'
##' @usage revcomp(x, update_db=TRUE)
##'
##' @param x A gbFeatureList or gbRecord object
##' (gbFeatureLists must include a 'source' field).
##' @param order Should the resulting gbFeatureList be reordered.
##' @param update_db Should the sequence and the new feature locations be
##' updated in the underlying filehash database?
##' 
##' @return A gbFeatureList object
##'
##' @docType methods
##' @export
setGeneric("revcomp", function(x, order=FALSE, update_db=FALSE, ...) 
  standardGeneric("revcomp") )


#' @export
setMethod("revcomp", "gbFeatureList",
          function(x, order=FALSE, update_db=FALSE) {
            .revcomp_features(x=x, order=order, update_db=update_db)
          })
