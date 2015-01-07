#' @include utils.R
#' @include all-generics.R
#' @include gbHeader-class.R
NULL

#' Class \code{"gbLocation"}
#' 
#' \dQuote{gbLocation} is an S4 class that provides a container for
#' GenBank feature location descriptors.
#' 
#' @slot range An integer matrix indicating the base numbers delimiting a
#' sequence span.
#' @slot fuzzy A logical matrix indicating fuzzy start and/or end
#' (e.g. <1..200).
#' @slot strand An integer vector containing -1, 1, or NA.
#' @slot compound A character code specifying how multiple ranges
#' are joined. One of \sQuote{join}, \sQuote{order}, or \sQuote{bond}.
#' @slot accession A character vector; the accession number of the sequence
#' of the feature this location lives on.
#' @slot remote A logical vector
#' @slot type A character vector describing the type of the position. Normally
#' an "R" for \sQuote{range} (e.g., \code{1..200} or point position \code{200}),
#' a "B" for \sQuote{between bases} (e.g., \code{36^37}), or a "G" for gaps
#' (e.g., \code{gap()}, \code{gap(30)}, or \code{gap(unk30)}).
#'
#' @details
#' For more information see the 
#' \href{ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt}{GenBank Release Note}
#' 
#' @keywords classes internal
setClass("gbLocation",
         representation(range = "matrix",
                        fuzzy = "matrix",
                        strand = "integer",
                        compound = "character",
                        accession = "character",
                        remote = "logical",
                        type = "character"),
         prototype(range = matrix(0L, 0, 2),
                   fuzzy = matrix(FALSE, 0, 2),
                   strand = NA_integer_,
                   compound = NA_character_,
                   accession = NA_character_,
                   remote = FALSE,
                   type = 'R'))


#' @keywords internal
#' @importFrom S4Vectors setValidity2
setValidity2("gbLocation", function(object) {
  # check range matrix
  if (!is.integer(object@range) || dim(object@range)[2] != 2 )
    return( "The 'range' slot should be a two-column, integer matrix." )
  # Check for valid base ranges
  if (any(object@range[, 2] < object@range[, 1], na.rm = TRUE))
    return( "One or more ranges with second endpoint before first." )
  # check fuzzy matrix
  if (!is.logical(object@fuzzy) || dim(object@fuzzy)[2] != 2 )
    return( "The 'fuzzy' slot should be a two-column, logical matrix." )
  # check strand vector
  if (all_empty(object@strand) || !all(object@strand %in% c(1L, -1L, NA_integer_)))
    return("The 'strand' slot should only contain 1L, -1L, or NA")
  # check compound character
  if (length(object@compound) > 1L || all_empty(object@compound) ||
        !object@compound %in% c("join", "order", "bond", NA_character_))
    return("The 'compound' slot should contain either 'join', 'order', 'bond', 'gap', or NA")
  # For type 'B', check that nucleotides are adjoining
  if (any(object@type == "B") && any(object@range[,2] - object@range[,1][object@type == 'B'] != 1))
    return( "For span type 'B', start and end position must be adjacent" )
  
  TRUE
})


# Getter-methods ---------------------------------------------------------


#' @describeIn start
setMethod("start", "gbLocation", function(x, join = FALSE) {
  if (join) {
    min(x@range[, 1, drop = TRUE])
  } else {
    x@range[, 1, drop = TRUE]
  }
})

#' @describeIn end
setMethod("end", "gbLocation", function(x, join = FALSE) {
  if (join) {
    max(x@range[, 2, drop = TRUE])
  } else {
    x@range[, 2, drop = TRUE]
  }
})

#' @describeIn span
setMethod("span", "gbLocation", function(x, join = FALSE) {
  if (join) {
    max(x@range[, 2]) - min(x@range[, 1]) + 1L
  } else {
    x@range[, 2] - x@range[, 1] + 1L
  }
})

#' @describeIn span
setMethod("joint_range", "gbLocation", function(x) {
  range(x@range)
})

#' @describeIn strand
setMethod("strand", "gbLocation", function(x, join = FALSE) {
  if (join || dim(x@range)[1] == 1L) {
    unique(x@strand)
  } else {
    x@strand
  }
})

#' @describeIn fuzzy
setMethod("fuzzy", "gbLocation", function(x) x@fuzzy)

#' @rdname accessors
setMethod("getAccession", "gbLocation", function(x) x@accession)


# Replace methods -----------------------------------------------------

.gbLocation_replace_start <- function(x, check = TRUE, value) {
  nrow <- dim(x@range)[1]
  if (!is.numeric(value))
    stop("replacement 'value' must be numeric")
  if (length(value) != nrow)
    stop("This gbLocation contains ", nrow, " start values")
  
  if (all(x@range[, 1] == x@range[, 2])) {
    x@range[, 1] <- as.integer(value)
    x@range[, 2] <- as.integer(value)
  } else {
    x@range[, 1] <- as.integer(value)
  }
  if (check)
    validObject(x)
  x
}

#' @describeIn start<-
setReplaceMethod("start", "gbLocation", function(x, ..., value) 
  .gbLocation_replace_start(x, ..., value = value)
)

.gbLocation_replace_end <- function(x, check = TRUE, value) {
  nrow <- dim(x@range)[1]
  if (!is.numeric(value))
    stop("replacement 'value' must be numeric")
  if (length(value) != nrow)
    stop("This gbLocation contains ", nrow ," end values")
  
  if (all(x@range[, 2] == x@range[, 1])) {
    x@range[, 2] <- as.integer(value)
    x@range[, 1] <- as.integer(value)
  } else {
    x@range[, 2] <- as.integer(value)
  }
  if (check)
    validObject(x)
  x
}

#' @describeIn end
setReplaceMethod("end", "gbLocation", function(x, ..., value) 
  .gbLocation_replace_end(x, ..., value = value)
)

#' @describeIn strand
setReplaceMethod("strand", "gbLocation",
                 function(x, value) {
                   nrow <- dim(x@range)[1]
                   if (length(value) > nrow)
                     value <- value[seq_len(nrow)]
                   if (length(value) < nrow)
                     value <- recycle(value, nrow)
                   if (is.character(value))
                     value <- vapply(value, switch, '+' = 1L, '-' = -1L, NA_integer_,
                                     FUN.VALUE = integer(1))
                   x@strand <- as.integer(value)
                   x
                 })


# Coerce-methods ------------------------------------------------------


setAs("gbLocation", "character",
      function(from) {
        nrow <- dim(from@range)[1]
        if (nrow == 0)
          return(character())
        else {
          rng <- from@range
          fuz <- from@fuzzy
          str <- from@strand
          cmp <- from@compound
          acc <- from@accession
          rem <- from@remote
          typ <- from@type
          span <- vapply(typ, switch, "R" = "..", "B" = "^", "G" = "",
                         FUN.VALUE = "", USE.NAMES = FALSE)
          pos <- ifelse(rng[,1] == rng[,2],
                        paste0(
                          ifelse(fuz[,1],
                                 "<",
                                 ifelse(fuz[,2], ">", "")
                          ),
                          ifelse(typ == "G", "", rng[,1])
                        ),
                        paste0(
                          ifelse(fuz[,1],
                                 ifelse(typ == "G", "unk", "<"),
                                 ""),
                          ifelse(typ == "G", "", rng[,1]),
                          span,
                          ifelse(fuz[,2], ">", "" ),
                          rng[,2]
                        )
          )
          pos[pos == "<"] <- ""
          pos <- ifelse(rem,
                        paste0(acc, ":", pos),
                        pos)
          pos <- ifelse(typ == "G", paste0("gap(", pos, ")"), pos)
          res <- 
            if (length(unique(str)) == 1) {
              paste0(
                ifelse( unique(str) == -1, "complement(", ""),
                ifelse( !is.na(cmp), paste0(cmp, "("), ""),
                paste0(pos, collapse = ","),
                ifelse( !is.na(cmp), ")", ""),
                ifelse( unique(str) == -1, ")", "")
              )
            } else if (length(str) == nrow) {
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
      function(from) gbLocation(from))


#' Create a \code{gbLocation}.
#' 
#' Create a \code{gbLocation} object out of a character string.
#' 
#' @param base_span A character string representation of GenBank feature location
#' @return A \code{\linkS4class{gbLocation}} object.
#' @export
#' @examples
#' as.gbLocation("join(1..10,12..20)")
as.gbLocation <- function(base_span) {
  as(as.character(base_span), "gbLocation")
}


# shift ---------------------------------------------------------------


#' @describeIn shift
setMethod("shift", "gbLocation", function(x, shift = 0L, ...) {
  if (!is.numeric(shift))
    stop("'shift' must be an integer")
  if (!is.integer(shift))
    shift <- as.integer(shift)
  if (length(shift) > 1L) {
    warning("'shift' must be a single integer. Only the first element is used")
    shift <- shift[[1L]]
  }
  
  x@range <- x@range + shift
  validObject(x)
  x
})


# Show-method ---------------------------------------------------------


setMethod("show", "gbLocation", function(object) {
  res <- as(object, "character")
  cat(linebreak(res, FORCE = TRUE), "\n" )
})

