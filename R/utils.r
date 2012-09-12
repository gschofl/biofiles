##' @importFrom Biostrings read.DNAStringSet
##' @importFrom Biostrings read.RNAStringSet
##' @importFrom Biostrings read.AAStringSet
##' @importFrom Biostrings DNAStringSet
##' @importFrom Biostrings RNAStringSet
##' @importFrom Biostrings AAStringSet
##' @importFrom Biostrings reverseComplement
##' @importFrom Biostrings xscat
##' @importFrom Biostrings subseq
##' @importFrom Biostrings toString 
##' @importFrom plyr rbind.fill
##' @importFrom parallel mcmapply
##' @importFrom parallel mclapply
##' @importFrom parallel detectCores
##' @import intervals
##' @import IRanges
##' @import stringr
NULL


##' Format paragraphs
##' 
##' Similar to \code{\link{strwrap}} but returns a single string with
##' linefeeds inserted
##' 
##' @param s a character vector or a list of character vectors
##' @param width a positive integer giving the column for inserting
##' linefeeds
##' @param indent an integer giving the indentation of the first line of
##' the paragraph; negative values of \code{indent} are allowed and reduce
##' the width for the first line by that value.
##' @param offset a non-negative integer giving the indentation of all
##' but the first line
##' @param split regular expression used for splitting. Defaults to
##' a whitespace character.
##' @param FORCE if \code{TRUE} words are force split if the available width
##' is too small.
##' @param FULL_FORCE Always split at the specified position.
##' 
##' @return a character vector
##' @keywords internal
linebreak <- function (s, width=getOption("width") - 2, indent=0, offset=0,
                       split=" ", FORCE=FALSE, FULL_FORCE=FALSE) {
  
  if (!is.character(s)) 
    s <- as.character(s)
  
  if (length(s) == 0L)
    return("")
  
  # set indent string to "" if a negative value is given
  # this lets us shrink the available width for the first line by that value
  indent_string <- blanks(ifelse(indent < 0, 0, indent))
  offset_string <- paste0("\n", blanks(offset))
  
  s <- mapply(function (s, width, offset, indent, indent_string, split, FORCE, FULL_FORCE) {
    # remove leading and trailing blanks
    # convert newlines, tabs, spaces to " "
    # find first position where 'split' applies
    if (!FULL_FORCE) {
      s <- gsub("[[:space:]]+", " ", gsub("^[[:blank:]]+|[[:blank:]]+$", "", s), perl=TRUE)
    }
    fws <- regexpr(split, s, perl=TRUE)
    if (offset + indent + nchar(s) > width) {
      # if not everything fits on one line
      if (FULL_FORCE ||
        (fws == -1 || fws >= (width - offset - indent)) && FORCE) {
        # if no whitespace or first word too long and force break
        # cut through the middle of a word
        pat1 <- paste0("^.{", width - offset - indent, "}(?=.+)")
        pat2 <- paste0("(?<=^.{", width - offset - indent, "}).+")
        leading_string <- regmatches(s, regexpr(pat1, s, perl=TRUE))
        trailing_string <- regmatches(s, regexpr(pat2, s, perl=TRUE)) 
        s <- paste0(indent_string, leading_string, offset_string,
                    linebreak(s=trailing_string, width=width, indent=0,
                              offset=offset, split=split, FORCE=FORCE,
                              FULL_FORCE=FULL_FORCE))
        
      } else if ((fws == -1 || fws >= (width - offset + indent)) && !FORCE) {
        # if no whitespace or first word too long and NO force break
        # stop right here
        stop("Can't break in the middle of a word. Use the force!")
        
      } else {
        # break the line
        s_split <- unlist(strsplit(s, split))
        s_cum <- cumsum(nchar(s_split) + 1)
        leading_string <- 
          paste0(paste0(s_split[s_cum < width - offset - indent], collapse=split),
                 ifelse(split == " ", "", split))
        trailing_string <- 
          paste0(s_split[s_cum >= width - offset - indent], collapse=split)
        s <- paste0(indent_string, leading_string, offset_string,
                    linebreak(s=trailing_string, width=width, indent=0,
                              offset=offset, split=split, FORCE=FORCE, FULL_FORCE=FULL_FORCE))
      }
      
    } else
      # if everything fits on one line go with the string
      s
  }, s, width, offset, abs(indent), indent_string, split, FORCE, FULL_FORCE,
              SIMPLIFY=FALSE, USE.NAMES=FALSE)
  unlist(s)
}


mergeLines <- function (lines) {
  if (length(lines) == 1) {
    gsub("^\\s+|\\s+$", "", lines)
  } else {
    paste0(gsub("^\\s+|\\s+$", "", lines), collapse=" ")
  }
}


##' create blank strings with a given number of characters
##' @seealso Examples for \code{\link{regmatches}}
##' @keywords internal
blanks <- function(n) {
  vapply(Map(rep.int, rep.int(" ", length(n)), n, USE.NAMES=FALSE),
         paste0, collapse="", character(1))
}


## chain functions
`%@%` <- function(x, f) {
  eval.parent(as.call(append(as.list(substitute(f)), list(x), 1)))
}


# negate %in%
`%ni%` <- Negate(`%in%`)


#' Flatten (Nested) Lists.
#'
#' Flatten \code{lists} according to specifications made via
#' \code{start_after} and/or \code{stop_at}. When keeping 
#' the defaults, the function will traverse \code{src} to retrieve the
#' values at the respective bottom layers/bottom elements. These values are
#' arranged in a named \code{list} where the respective names can be
#' interpreted as the the paths to the retrieved values.   
#'
#' @param x An arbitrarily deeply nested \code{list}
#' @param start_after An \code{integer} specifying the layer after which to 
#' start the flattening. \code{NULL} means to start at the very top.
#' @param stop_at An \code{integer} specifying the layer at which to stop
#' the flattening. \code{NULL} means there is not stop criterion.
#' @param delim_path A \code{character} specifying how the names
#' of the resulting flattened list should be pasted.
#' @param ... Further args.
#' @return A named \code{list} that features the desired degree of flattening.
#' @keywords internal
#' @author Janko Thyson \email{janko.thyson.rstuff@@googlemail.com}
#' @examples
#'  ##
flatten <- function (x, 
                     start_after=NULL, 
                     stop_at=NULL, 
                     delim_path=".",
                     do_warn=TRUE,
                     ... )
{
  # VALIDATE
  if (!is.list(x)) {
    stop("'src' must be a list.")
  }
  if (!is.null(start_after) && !is.null(stop_at)) {
    if (start_after == 1 && stop_at == 1)
      stop(sprintf("Invalid specification:\nstart_after: %s\nstop_at: %s\n",
                   start_after, stop_at))
  }
  
  # INNER FUNCTIONS
  .startAfterInner <- function(envir, nms, out.1, ...) {
    idx_diff <- diff(c(envir$start_after, length(envir$counter)))
    
    # UPDATE IF DEGREE OF NESTEDNESS EXCEEDS START CRITERION
    if (idx_diff > 0) {
      idx_cutoff <-
        seq(from=(length(envir$counter) - idx_diff + 1), to=length(envir$counter))
      
      idx_left <- envir$counter[-idx_cutoff]
      nms.1 <- nms[idx_cutoff]
      names(out.1) <- paste(nms.1, collapse=envir$delim_path)
      # UPDATE SRC
      idx_append <- sapply(envir$history, function (x_hist) {
        all(idx_left == x_hist)        
      })
      
      if (any(idx_append)) {                                          
        envir$src[[idx_left]] <- append(envir$src[[idx_left]], values=out.1)                    
      } else {
        envir$src[[idx_left]] <- out.1
        # UPDATE HISTORY
        envir$history <- c(envir$history, list(idx_left))
      }
      envir$out <- envir$src          
    } else if (idx_diff < 0) {
      envir$out <- envir$src
    }
    
    # RESET
    envir$nms <- envir$nms[-length(envir$nms)]
    envir$counter <- envir$counter[-length(envir$counter)]
    
    TRUE
  }
  
  .updateOutInner <- function (envir, out.1, ...) {
    
    # UPDATE OUT
    envir$out <- c(get("out", envir = envir), out.1)
    
    # RESET
    envir$nms       <- envir$nms[-length(envir$nms)]
    envir$counter   <- envir$counter[-length(envir$counter)]
    
    TRUE
  }
  
  .flattenInner <- function(x, envir, ...) {
    if ( is(x, "list") && length(x) != 0 ) {
      
      # UPDATE
      envir$counter_history <- c(envir$counter_history, list(envir$counter))
      
      # EXIT IF DEGREE EXCEEDS CUTOFF
      if (!is.null(envir$stop_at)) {
        if (length(envir$counter) > envir$stop_at) { 
          nms <- get("nms", envir=envir)
          out.1 <- list(x)
          names(out.1) <- paste(nms, collapse=envir$delim_path)
          
          # DECISION ON FLATTENING
          if (!is.null(envir$start_after)) {
            .startAfterInner(envir=envir, nms=nms, out.1=out.1)
            return(NULL)
          } else {
            .updateOutInner(envir=envir, out.1=out.1)
            return(NULL)
          }
        }
      }
      
      # LOOP OVER ELEMENTS
      for (i in seq_along(x)) {
        # UPDATE COUNTER
        envir$counter <- c(envir$counter, i)
        # UPDATE NAMES
        list_names <- if (is.null(names(x[i]))) paste0("X", i) else names(x[i])
        assign("nms", c(get("nms", envir=envir), list_names), envir=envir)
        # RECURSIVE FLATTENING
        .flattenInner(x=x[[i]], envir) # call  recursively
        # RESET COUNTER
        if (i == length(x)) {
          envir$nms <- envir$nms[-length(envir$nms)]
          envir$counter <- envir$counter[-length(envir$counter)]
        }
      }
    } else {
      
      nms <- get("nms", envir=envir)
      out.1 <- list(x)
      names(out.1) <- paste(nms, collapse=envir$delim_path)
      
      # DECISION ON FLATTENING
      if (!is.null(envir$start_after))
        .startAfterInner(envir=envir, nms=nms, out.1=out.1)
      else
        .updateOutInner(envir=envir, out.1=out.1)
    }
    
    TRUE
  }
  
  out <- list()
  # ENVIR
  envir <- new.env()
  envir$counter <- NULL
  envir$counter_history <- NULL
  envir$delim_path <- delim_path
  envir$do_warn <- do_warn
  envir$do_block_warning <- FALSE
  envir$history <- NULL
  envir$nms <- NULL
  envir$out <- list()
  envir$src <- x
  envir$start_after <- start_after
  
  if (!is.null(stop_at)) {
    stop_at_0 <- stop_at
    if (stop_at == 1) {
      return(src)
    } else {
      stop_at <- stop_at - 1
    }
  }
  
  envir$stop_at <- stop_at
  
  .flattenInner(x, envir)
  
  if (envir$do_warn) {
    max_length <- max(sapply(envir$counter_history, length))
    
    if (!is.null(start_after)) {            
      if (start_after > max_length) {                        
        warning(paste("Argument 'start_after=", start_after, 
                      "' exceeds maximum degree of sublayer nestedness (=", 
                      max_length, ").", sep=""))
      }
    }
    if (!is.null(stop_at)) {
      if (stop_at_0 > max_length){
        warning(paste("Argument 'stop_at=", stop_at_0, 
                      "' exceeds maximum degree of sublayer nestedness (=", 
                      max_length, ").", sep=""))    
      }
    }
  }
  out <- envir$out
  out    
}


getCompounds <- function (x) {
  x <- x[which(is.compound(x))]
  if (length(x) == 0) return(NA_real_) 
  cL <- vapply(x, function (f) nrow(f@location@.Data), numeric(1))
  cL
}


expandIds <- function (x) {
  cmp_pos <- Position(is.compound, x)
  if (is.na(cmp_pos)) {
    return(list(ids=getIndex(x, FALSE), keys=getKey(x, FALSE)))
  }
  xHead <- x[1:(cmp_pos - 1)]
  xCmp <- x[cmp_pos]
  L <- length(x)
  xTail <- if (cmp_pos + 1 <= L) {
    x[(cmp_pos +  1):L]
  } else {
    .gbFeatureList()
  }
  id1 <- list(ids=getIndex(xHead, FALSE), keys=getKey(xHead, FALSE))
  nC <- getCompounds(xCmp)
  id2 <- list(ids=rep(getIndex(xCmp, FALSE), nC),
              keys=rep(getKey(xCmp, FALSE), nC))
  id3 <- list(ids=getIndex(xTail, FALSE), keys=getKey(xTail, FALSE))
  Map(function(a, b, c) c(a, b, c), a=id1, b=id2, c=Recall(xTail))
}

## Access qualifiers from gbFeature or gbFeatureList objects
.qualAccess <- function (x, qual = "", fixed = FALSE) {
  
  .access <- function (q) {
    q <- q@qualifiers
    if (fixed) qual <- paste0("\\b", qual, "\\b") 
    n <- length(q)
    
    if (n == 0) {
      return( structure(rep(NA_character_, length(qual)),
                        names=gsub("\\b", "", qual, fixed=TRUE)) )
    }
    
    idx <- matrix(
      vapply(qual, grepl, names(q),
             USE.NAMES=FALSE, FUN.VALUE=logical(n)),
      nrow=n)
    n_col <- dim(idx)[2]
    if (n_col == 1L) {
      if (any(idx))
        q[idx]
      else
        structure(NA_character_, names=gsub("\\b", "", qual, fixed=TRUE))
    } else {
      ans <- lapply(seq.int(n_col), function (i) q[ idx[ ,i]])
      if (any(na_idx <- !apply(idx, 2, any))) {
        for (na in which(na_idx)) {
          ans[[na]] <- structure(NA_character_, names=gsub("\\b", "", qual, fixed=TRUE)[na])
        }
      }
      unlist(ans)
    }
  }
  
  if (is(x, "gbFeature")) {
    .access(x)
  } else if (is(x, "gbFeatureList")) {
    lapply(x, .access)
  }
}


.simplify <- function (x, unlist = TRUE) {
  if (length(len <- unique(unlist(lapply(x, length)))) > 1L) {
    return(x)
  }
  if (len == 1L && unlist) {
    unlist(x, recursive=FALSE)
  } else if (len >= 1L) {
    n <- length(x)
    r <- as.vector(unlist(x, recursive=FALSE))
    if (prod(d <- c(len, n)) == length(r)) {
      return(data.frame(stringsAsFactors=FALSE,
                        matrix(r, nrow=n, byrow=TRUE,
                               dimnames=if (!(is.null(nm <- names(x[[1L]]))))
                                 list(NULL, nm))))
    } else {
      x
    }
  } else {
    x
  }
}



.seqAccess <- function (s, x, type) {
  ## merge Sequences
  mergeSeq <- function (s, x, type) {
    if (length(start(x)) == 1L) {
      seq <- subseq(s, start=start(x), end=end(x))
    } else {
      seq <- do.call(xscat, Map(subseq, s, start=start(x), end=end(x)))
    }
    seq <- switch(type,
                  DNA=DNAStringSet(seq),
                  AA=AAStringSet(seq),
                  RNA=RNAStringSet(seq))
    seq@ranges@NAMES <- sprintf("%s.%s.%s", x@.ACCN, x@key, x@.ID)
    seq
  }
  
  if (is(x, "gbFeatureList")) {
    ## initiate empty XStringSet
    seq <- switch(type, 
                  DNA=DNAStringSet(),
                  AA=AAStringSet(),
                  RNA=RNAStringSet())
    for (i in seq_along(x)) {
      seq[i] <- mergeSeq(s, x[[i]], type)             
    }
  } else if (is(x, "gbFeature")) {
    seq <- mergeSeq(s, x, type)
  }
  
  seq@metadata <- list(definition=x@.DEF, database=x@.Dir)
  seq
}

# --R-- vim:ft=r:sw=2:sts=2:ts=4:tw=76:
#       vim:fdm=marker:fmr={{{,}}}:fdl=0

