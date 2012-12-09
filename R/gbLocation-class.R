#' @include utils.r
#' @include validate.r
#' @include all-generics.r
#' @importClassesFrom intervals Intervals_full
NULL


# gbLocation-class ----------------------------------------------------


#' gbLocation
#' 
#' \dQuote{gbLocation} is an S4 class that provides a container for
#' GenBank feature locations.
#' It extends \code{\linkS4class{Intervals_full}} and provides 5
#' additional Slots.
#' 
#' @slot strand A single integer -1 or 1 (for minus or plus strand,
#'   respectively)
#' @slot compound A character code specifying how multiple segments
#'   are joined. One of \sQuote{join} or \sQuote{order}.
#' @slot partial A logical matrix specifying whether residues are
#'   missing from the 5' and 3' ends respectively.
#' @slot accession
#' @slot remote
#'
#' @details
#' For more information see the 
#' \href{ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt}{GenBank Release Note}
#'
#' @rdname gbLocation
#' @export
#' @classHierarchy
#' @classMethods
setClass("gbLocation",
         representation(strand = "integer",
                        compound = "character",
                        partial = "matrix",
                        accession = "character",
                        remote = "logical"),
         prototype(type = "Z",
                   strand = NA_integer_,
                   compound = NA_character_,
                   partial = matrix( FALSE, 0, 2 ),
                   accession = NA_character_,
                   remote = FALSE ),
         contains = "Intervals_full")


#' @keywords internal
#' @autoImports
setValidity2("gbLocation",
            function (object) {  
              if (is_empty(object@strand) || !all(object@strand %in% c(1L, -1L, NA_integer_)))
                return("The 'strand' slot can contain either -1, 1, or NA")         
              if (length(object@compound) > 1L || is_empty(object@compound) ||
                   object@compound %ni% c("join", "order", NA_character_))
                return("The 'compound' slot can contain either 'join', 'order', or NA")
              
              TRUE
            })


#' @keywords internal
#' @autoImports
setMethod("initialize", "gbLocation",
          function (.Object, .Data, strand, compound, partial, remote, ...)  {
            if (missing(.Data)) {
              callNextMethod(.Object, ...)
            } else {
              if (missing(strand)) {
                strand <- NA_integer_
              } 
              if (missing(partial)) {
                partial <- matrix(FALSE, nrow(.Data), 2)
              }
              if (missing(compound)) {
                compound <- NA_character_
              }
              if (missing(remote)) {
                remote <- rep(FALSE, nrow(.Data))
              } else {
                remote <- recycle(remote, nrow(.Data))
              } 
              callNextMethod(.Object, .Data=.Data, strand=strand,
                             compound=compound, partial=partial,
                             remote=remote, ...)
            }
          })


# Getter-methods ---------------------------------------------------------


setMethod("start", "gbLocation",
          function (x, join = FALSE, drop = TRUE) {
            if (join)
              min(x@.Data[, 1, drop = drop])
            else
              x@.Data[, 1, drop = drop]
          })


setMethod("end", "gbLocation",
          function (x, join = FALSE, drop = TRUE) {
            if (join)
              max(x@.Data[, 2, drop = drop])
            else
              x@.Data[, 2, drop = drop]
          })


setMethod("width", "gbLocation",
          function (x, join = FALSE) {
            if (join) 
              max(x@.Data[, 2]) - min(x@.Data[, 1]) + 1
            else
              x@.Data[, 2] - x@.Data[, 1] + 1
          })


setMethod("strand", "gbLocation",
          function (x, join = FALSE) {
            if (join || nrow(x) == 1L)
              unique(x@strand)
            else
              x@strand     
          })


setMethod("range", "gbLocation",
          function (x, join = FALSE) {
            start <- as.integer(start(x, join = join))
            width <- as.integer(end(x, join = join)) - start + 1L
            strand <- strand(x, join = join)
            new('gbRange', start, width, strand)
          })


setMethod("partial", "gbLocation",
          function (x) x@partial)


setMethod("accession", "gbLocation",
          function (x) x@accession)


# Replace methods -----------------------------------------------------


setReplaceMethod("start", "gbLocation",
                 function (x, value) {
                   if (!is.numeric(value))
                     stop("replacement 'value' must be numeric")
                   if (length(value) != nrow(x)) {
                     stop(sprintf("This gbLocation contains %s start values", nrow(x)))
                   }
                   if (all(x@.Data[,1] == x@.Data[,2])) {
                     x@.Data[,1] <- value
                     x@.Data[,2] <- value
                   } else {
                     x@.Data[,1] <- value
                   }
                   validObject(x)
                   x
                 })


setReplaceMethod("end", "gbLocation",
                 function (x, value) {
                   if (!is.numeric(value))
                     stop("replacement 'value' must be numeric")
                   if (length(value) != nrow(x)) {
                     stop(sprintf("This gbLocation contains %s end values", nrow(x)))
                   }
                   if (all(x@.Data[,2] == x@.Data[,1])) {
                     x@.Data[,2] <- value
                     x@.Data[,1] <- value
                   } else {
                     x@.Data[,2] <- value
                   }
                   validObject(x)
                   x
                 })


setReplaceMethod("strand", "gbLocation",
                 function (x, value) {
                   if (length(value) > nrow(x)) {
                     value <- value[seq_len(nrow(x))]
                   }
                   if (length(value) < nrow(x)) {
                     value <- recycle(value, nrow(x))
                   }
                   if (is.character(value) && all(value %in% c("+", "-", NA_character_))) {
                     value <- vapply(value, switch, "+" = 1L, "-" = -1L,
                                     "NA" = NA_integer_, FUN.VALUE=integer(1))
                   }
                   x@strand <- as.integer(value)
                   validObject(x)
                   x
                 })

# Coerce-methods ------------------------------------------------------


#' @autoImports
#' @importFrom intervals closed
#' @importFrom intervals partial
setAs("gbLocation", "character",
      function (from) {
        if (nrow(from) == 0)
          return(character())
        else {
          clo <- from@closed
          par <- partial(from)
          str <- strand(from)
          cmp <- from@compound
          acc <- accession(from)
          rem <- from@remote
          
          span <- ifelse(clo[,1],
                         "..", 
                         ifelse(from[,2] == from[,1] + 1,
                                "^",
                                ".")
          )
          
          pos <- ifelse(from[,1] == from[,2],
                        from[,1], 
                        paste0(
                          ifelse( par[,1], "<", "" ),
                          from[,1],
                          span,
                          ifelse( par[,2], ">", "" ),
                          from[,2]
                        )
          )
          
          pos <- ifelse( rem,
                         paste0(acc, ":", pos),
                         pos)
          
          res <- 
            if (length(unique(str)) == 1) {
              paste0(
                ifelse( unique(str) == -1, "complement(", ""),
                ifelse( !is.na(cmp), paste0(cmp, "("), ""),
                paste0(pos, collapse=","),
                ifelse( !is.na(cmp), ")", ""),
                ifelse( unique(str) == -1, ")", "")
              )
            } else if (length(str) == nrow(from)) {
              paste0(
                ifelse( !is.na(cmp), paste0(cmp, "("), ""),
                paste0(
                  ifelse( str == -1L,
                          paste0("complement(", pos, ")"),
                          pos),
                  collapse = ","),
                ifelse( !is.na(cmp), ")", "")
              )  
            }
          
          res
        }
      })


setAs("character", "gbLocation",
      function (from) {
        l <- parse_gb_location(from)
        if (is_empty(l)) {
          err <- sprintf("The string %s cannot be parsed as a gbLocation.",
                         sQuote(from))
          stop(err)
        }
        l
      })


#' @export
as.gbLocation <- function (base_span) {
  as(as.character(base_span), "gbLocation")
}


# shift ---------------------------------------------------------------


setMethod("shift", "gbLocation",
          function (x, shift = 0L, ...) {
            if (!is.numeric(shift))
              stop("'shift' must be an integer")
            if (!is.integer(shift))
              shift <- as.integer(shift)
            if (length(shift) > 1L) {
              warning("'shift' must be a single integer. Only the first element is used")
              shift <- shift[[1L]]
            }
            
            x@.Data <- x@.Data + shift
            validObject(x)
            x
          })


# Show-method ---------------------------------------------------------


setMethod("show", "gbLocation",
          function (object) {
            res <- as(object, "character")
            cat(linebreak(res, FORCE=TRUE), "\n" )
          })

