#### gbFeatureList objects
setOldClass("list")

setClass("gbFeatureList", 
         representation(.Dir="character",
                        .ACCN="character",
                        .DEF="character"),
         contains="list")

#### Constructor
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
  
  fl <- new("gbFeatureList", .Data=features, .Dir=as.character(db_dir),
            .ACCN=as.character(accession), .DEF=as.character(definition))
  fl
}

#### Accessors
setMethod("index", "gbFeatureList",
          function (object) {
            idx <- vapply(object, index, FUN.VALUE=0)
            attr(idx, "accession") <- object@.ACCN
            attr(idx, "definition") <- object@.DEF
            attr(idx, "db_dir") <- object@.Dir
            idx
          })

setMethod("key", "gbFeatureList",
          function (object) {
            vapply(object, key, FUN.VALUE="")
          })

setMethod("start", "gbFeatureList", 
          function (object) {
            pos <- lapply(object, .location, "start")
            if (any(vapply(pos, length, FUN.VALUE=0) > 1)) {
              message("Features have multiple start postitions. Cannot return vector")
              return(pos)
            } else {
              return(unlist(pos))
            }
          })

setMethod("end", "gbFeatureList", 
          function (object) {
            pos <- lapply(object, .location, "end")
            if (any(vapply(pos, length, FUN.VALUE=0) > 1)) {
              message("Features have multiple end postitions. Cannot return vector")
              return(pos)
            } else {
              return(unlist(pos))
            }
          })

setMethod("strand", "gbFeatureList", 
          function (object) {
            strand <- vapply(object, .location, "strand", FUN.VALUE=0)
            strand
          })

setMethod("qualifiers", "gbFeatureList", 
          function (object) {
            lapply(object, qualifiers)
          })

setMethod("locusTag", "gbFeatureList", 
          function (object) {
            vapply(object, locusTag, FUN.VALUE="")
          })

setMethod("gene", "gbFeatureList", 
          function (object) {
            vapply(object, gene, FUN.VALUE="")
          })

setMethod("product", "gbFeatureList", 
          function (object) {
            vapply(product, gene, FUN.VALUE="")
          })

setMethod("note", "gbFeatureList", 
          function (object) {
            vapply(product, note, FUN.VALUE="")
          })

setMethod("proteinId", "gbFeatureList", 
          function (object) {
            vapply(product, proteinId, FUN.VALUE="")
          })

setMethod("dbxref", "gbFeatureList", 
          function (object) {
            lapply(product, dbxref)
          })

setMethod("translation", "gbFeatureList", 
          function (object) {
            
            translation2 <- function (object) {
              idx <- charmatch("translation", names(object@qualifiers))
              if (is.na(idx))
                return(NA_character_)
              else
                return(object@qualifiers[[idx]])
            }
            
            aa <- AAStringSet(vapply(object, translation2, FUN.VALUE=""))
            names(aa) <- locusTag(object)
            aa
          })

setMethod("is.pseudo", "gbFeatureList", 
          function (object) {
            vapply(object, is.pseudo, FUN.VALUE=FALSE)
          })

#### show method
setMethod("show", signature="gbFeatureList", 
          function (object)  {
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


#### Subsetting
setMethod("[", signature(x="gbFeatureList", i="character", j="missing", drop="missing"),
          function(x, i, j) {
            idx <- vapply(x@.Data, function(f) f@key, FUN.VALUE="") == i
            gbFeatureList(x@.Dir, x@.ACCN, x@.DEF, x@.Data[idx])
          })

setMethod("[", signature(x="gbFeatureList", i="numeric", j="missing", drop="missing"),
          function(x, i, j) {
            gbFeatureList(x@.Dir, x@.ACCN, x@.DEF, x@.Data[i])
          })

setMethod("[", signature(x="gbFeatureList", i="logical", j="missing", drop="missing"),
          function(x, i, j) {
            gbFeatureList(x@.Dir, x@.ACCN, x@.DEF, x@.Data[i])
          })

setMethod("[", signature(x="gbFeatureList", i="missing", j="missing", drop="missing"),
          function(x, i, j) {
            return(x)
          })

#### select method
setMethod("select", signature(db="gbFeatureList"), 
          definition=function (db, key, qualifier, location)  {
            .select(x=db, key, qualifier, location) 
          })
