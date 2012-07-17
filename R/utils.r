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
##' @import intervals
##' @import IRanges
##' @import stringr
NULL


##' Join lines of text recursively
##' 
##' @param lines A character vector
##' @param extract_pat A regular expression extracting parts of the
##' character strings to join (defaults to \code{.*})
##' @param break_pat A regular expression that sets a break condition.
##' The line where this pattern matches is the last to be joined. 
##' If not set, all lines will be joined
##' @param sep if \code{FALSE} the joined lines are concatenated without
##' intervening spaces 
##' 
##' @return A list with two elements. The first contains the joined
##' character vector, the second the number of lines joined
##' 
##' @keywords internal
joinLines <- function (lines, extract_pat = ".*", break_pat = NULL, sep = TRUE) {
  i  <-  0
  list(eval(function (lines, extract_pat, break_pat) {
    l <- regmatches(lines[1], regexpr(extract_pat, lines[1], perl=T))
    i <<- i + 1
    
    if (length(l) == 0) {
      # jump out if we reach the last element of the character vector
      i <<- i - 1
      return(l) }
    else if (!is.null(break_pat) && grepl(break_pat, l, perl=TRUE))
      # or if a break pattern is set jump out when the break condition is met 
      return(l)
    
    if (sep)
      l <- paste(l, Recall(lines[-1], extract_pat, break_pat))
    else
      l <- paste0(l, Recall(lines[-1], extract_pat, break_pat))
  }) (lines, extract_pat=extract_pat, break_pat=break_pat), i)
}

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

##' Extract matched group(s) from a string.
##'
##' @param pattern character string containing a regular expression
##' @param str character vector where matches are sought
##' @param capture if \code{TRUE} capture groups are returned in addition
##' to the complete match
##' @param perl if \code{TRUE} perl-compatible regexps are used.
##' @param global if \code{TRUE} \code{gregexpr} is used for matching
##' otherwise \code{regexpr}.
##' @param ignore.case case sensitive matching
##' @return a list containing a \code{match} and a \code{capture} component
##' @keywords character
##' @keywords internal
##' @examples
##' ##
strmatch <- function (pattern, str, capture = TRUE, perl = TRUE, 
                      global = TRUE, ignore.case = FALSE) {
  
  if (!is.atomic(str))
    stop("String must be an atomic vector", call. = FALSE)
  
  if (!is.character(str)) 
    string <- as.character(str)
  
  if (!is.character(pattern)) 
    stop("Pattern must be a character vector", call. = FALSE)
  
  if (global)
    m <- gregexpr(pattern, str, perl=perl, ignore.case=ignore.case)
  else
    m <- regexpr(pattern, str, perl=perl, ignore.case=ignore.case)
  
  .matcher <- function (str, m) {
    Map( function (str, start, len) substring(str, start, start + len - 1L), 
         str, m, lapply(m, attr, "match.length"), USE.NAMES=FALSE)
  }
  
  match <- if (capture) {
    .capture.matcher <- function (str, m) {
      cap <- Map( function (str, start, len) {
        mapply( function (str, start, len) {
          substr(str, start, start + len - 1L) 
        }, str, start, len, USE.NAMES=FALSE)
      }, str, lapply(m, attr, "capture.start"),
                  lapply(m, attr, "capture.length"), USE.NAMES=FALSE)
      
      cap_names <- lapply(m, attr, "capture.names")
      if (all(nchar(cap_names) > 0)) {
        if (!all(mapply(function (c, n) length(c) == length(n), cap, cap_names)))
          warning("Mismatch between number of captures and capture names", call.=TRUE)
        
        cap <- mapply( function (val, name) `names<-`(val, name),
                       cap, cap_names, USE.NAMES=FALSE)
      }
      
      cap
    }
    
    list(match=.matcher(str, m),
         capture=if (!is.null(attributes(m[[1]])$capture.start))
           .capture.matcher(str, m) else NULL)
  } else {
    match <- .matcher(str, m)
  }
  match
}

.typeToFeature <- function (s) {
  s <- gsub("gene_component_region", "misc_feature", s)
  s
}

# substitute gff-specific attribute tags for GenBank qualifiers
.tagToQualifier <- function (s) {
  s <- sub("^ID$", "locus_tag", s)
  s <- sub("^Name$", "gene", s)
  s <- sub("^Dbxref$", "db_xref", s)
  s <- tolower(s)
  s <- sub("^", "/", s)
  s
}

## TAB - %09, NL - %0A, CR - %0D, 
## %3B (semicolon), %3D (equals),
## %25 (percent), %26 (ampersand)
## %2C (comma)
.unescape <- function (s) {
  s <- gsub("%3B", ";", s)
  s <- gsub("%3D", "=", s)
  s <- gsub("%25", "%", s)
  s <- gsub("%26", "&", s)
  s <- gsub("%2C", ",", s)
  s <- gsub("^", "\"", s)
  s <- gsub("$", "\"", s)
  s
}

## introduce line breaks after 80 cols
## the negative indent accounts for the length of the qualifier tag
## and the following "="
.cleanQualifiers <- function (q, v) {
  
  v <- linebreak(s=.unescape(v), width=79, offset=21,
                 indent=-(nchar(q)+1), FORCE=TRUE)
  
  l_pos <- which(pmatch(q, "/locus_tag", dup=TRUE, nomatch=0) == 1)
  p_pos <- which(pmatch(q, "/parent", dup=TRUE, nomatch=0) == 1)
  g_pos <- which(pmatch(q, "/gene", dup=TRUE, nomatch=0) == 1)
  l_val <- v[l_pos]
  p_val <- v[p_pos]
  g_val <- v[g_pos]
  lp <- c(l_val, p_val)
  lpg <- c(l_val, p_val, g_val)
  
  if (length(g_val) != 0 && 
    !any(grepl(pattern=g_val, x=c(l_val, p_val)))) {
    # if gene exists and is different from locus tag
    # retain gene and choose locus tag among parent and locus tag
    # and throw out parent
    v[l_pos] <- lp[which(nchar(lp) == min(nchar(lp)))]
    v <- v[-c(l_pos[-1], p_pos)]
    q <- q[-c(l_pos[-1], p_pos)]
  } else {
    # otherwise choose among locus tag, parent, and gene; throw
    # out parrent and gene
    v[l_pos] <- lpg[which(nchar(lpg) == min(nchar(lpg)))]
    v <- v[if (length(s <- c(l_pos[-1], p_pos, g_pos)) > 0) -s else TRUE]
    q <- q[if (length(s <- c(l_pos[-1], p_pos, g_pos)) > 0) -s else TRUE]
  }
  if (length(v) > 0) 
    return(list(q=q, v=v))
  else
    invisible(NULL)
}

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

