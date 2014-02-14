#' @import methods
#' @import Rcpp
#' @importFrom assertthat assert_that "on_failure<-" is.string is.scalar
#' @importFrom stats setNames
NULL

is.empty <- function(x) {
  is.null(x) || length(x) == 0L || (length(x) == 1L && !nzchar(x))
}
on_failure(is.empty) <- function(call, env) {
  paste0(deparse(call$x), " is not empty.")
}

are_empty <- function(x) {
  if (is.recursive(x) || length(x) > 1) {
    vapply(x, function(x) is.null(x) || length(x) == 0L, FALSE, USE.NAMES = FALSE) | !nzchar(x)
  } else {
    is.empty(x)
  }
}

all_empty <- function(x) all(are_empty(x))
on_failure(all_empty) <- function(call, env) {
  paste0("Not all elements in ", deparse(call$x), " are empty.")
}

is_in <- function(x, table) {
  assert_that(is.scalar(x))
  x %in% table
}
on_failure(is_in) <- function(call, env) {
  paste0(sQuote(deparse(call$x)), " is not an element of ",
         paste0(sQuote(eval(call$table, env)), collapse = ", "))
}

"%is_in%" <- is_in

"%ni%" <- Negate(`%in%`)

"%||%" <- function(a, b) {
  if (is.empty(a)) force(b) else a
}

"%|NA|%" <- function(a, b) {
  if (is.na(a)) force(b) else a
}

"%|AA|%" <- function(a, b) {
  if (a == 'AA') force(b) else a
}

dots <- function(...) {
  eval(substitute(alist(...)))
}

compact <- function(x) {
  x[!vapply(x, is.empty, FALSE, USE.NAMES = FALSE)]
}

compactChar <- function(x) {
  x[vapply(x, nzchar, FALSE, USE.NAMES = FALSE)]
}

compactNA <- function(x) {
  x[!is.na(x)]
}

merge_dups <- function(x) {
  if (all_empty(x)) {
    return(NULL)
  }
  x_names <- names(x)
  a <- x[!duplicated(x_names)]
  b <- x[duplicated(x_names)]
  modify_list(a, b, "merge")
}

modify_list <- function(a, b, mode = c("replace", "merge")) {
  assert_that(is.list(a), is.list(b))
  mode <- match.arg(mode)
  a_names <- names(a)
  for (v in names(b)) {
    a[[v]] <- if (v %in% a_names && is.list(a[[v]]) && is.list(b[[v]])) {
      modify_list(a[[v]], b[[v]])
    } else {
      switch(mode,
             replace = b[[v]],
             merge = unique(c(a[[v]], b[[v]])))
    }
  }
  a
}

#' @importFrom reutils make_flattener
flatten1 <- make_flattener(flatten.at = 1)
flatten2 <- make_flattener(flatten.at = 2)

re <- function(x) {
  assert_that(is.string(x))
  structure(x, class = "regexp")
}

"%~%" <- function(x, pattern) {
  grepl(pattern, x, fixed = FALSE)
}

## divide the data in vector x into groups, where the start of
## each group is defined by an element of index i
##
## e.g. x = (a, b, c, d, e, f, g)
##      i = (1, 3, 6)
##    res = (a,b), (c,d,e), (f,g)  
ixsplit <- function(x, i, include_i = TRUE, collapse_x = FALSE, ...) {
  l <- length(x)
  if (max(i) > l) {
    stop("index '", max(i), "' out of range")
  }
  j <- c(i[-1] - 1, l)
  if (!include_i) {
    i <- i + 1
  }
  if (any(i > j)) {
    stop("start point large than end point")
  }
  FUN <- if (collapse_x) {
    function(i) collapse(x[i], ...)
  } else {
    function(i) x[i]
  }
  lapply(.mapply(seq.int, list(i, j), NULL), FUN)
}

Call <- function(fn, ...) {
  fn <- match.fun(fn)
  fn(...)
}

Partial <- function(fn, ..., .env = parent.frame()) {
  fn <- match.fun(fn)
  fcall <- substitute(fn(...))
  if (!is.primitive(fn))
    fcall <- match.call(fn, fcall)  
  fcall[[length(fcall) + 1]] <- quote(...)
  args <- list("..." = quote(expr = ))
  eval(call("function", as.pairlist(args), fcall), .env)
}

Compose <- function(...) {
  fns <- lapply(list(...), match.fun)
  len <- length(fns)
  function(...) {
    res <- Call(fns[[len]], ...)
    for (fn in rev(fns[-len]))
      res <- fn(res)
    res
  }
}

usplit <- Compose("unlist", "strsplit")

uusplit <- Compose("unique", "unlist", "strsplit")

dup <- function(x, n) {
  if (any(n < 0)) n[n < 0] <- 0
  vapply(.mapply(rep.int, list(rep.int(x, length(n)), n), NULL),
         paste0, collapse = "", FUN.VALUE = "")
}

blanks <- Partial(dup, x = " ")

trim <- function(x, trim = '\\s+') {
  gsub(paste0("^", trim, "|", trim, "$"), '', x)
}

wrap <- function(x, wrap = '"') {
  sprintf('%s%s%s', wrap, x, wrap)
}

## UnitTests: inst/tests/test-utils.r
collapse <- function(x, sep = ' ') {
  if (is.list(x)) {
    vapply(x, collapse, sep = sep, FUN.VALUE = '')
  } else {
    paste0(trim(x), collapse = sep)
  }
}

ellipsize <- function(obj, width = getOption("width"), ellipsis = " ...") {
  str <- encodeString(obj)
  ifelse(nchar(str) > width - 1,
         paste0(substring(str, 1, width - nchar(ellipsis) - 1), ellipsis),
         str)
}

#' Pad a string
#' 
#' @param x Input character vector.
#' @param n Pad \code{x} to this (minimum) width.
#' @param where Side where the padding is added.
#' @param pad Padding character.
#' @return A character vector.
#' @keywords internal
pad <- function (x, n = 10, where = 'left', pad = ' ') {
  x <- as.character(x)
  where <- match.arg(where, c("left", "right", "both"))
  needed <- pmax(0, n - nchar(x))
  left <- switch(where, left = needed, right = 0, both = floor(needed/2))
  right <- switch(where, left = 0, right = needed, both = ceiling(needed/2))
  lengths <- unique(c(left, right))
  padding <- dup(pad, lengths)
  paste0(padding[match(left, lengths)], x, padding[match(right, lengths)])
}

count_re <- function(x, re) {
  vapply(gregexpr(re, x), function(x) sum(x > 0L), 0, USE.NAMES = FALSE)
}

#' Test if an external executable is available
#' 
#' Uses \code{\link{Sys.which}} internally, so it should work
#' on Windows and Unix.alikes.
#' 
#' @param cmd The exececutable to test for.
#' @param msg Additional message if the test fails.
#' @keywords internal
has_command <- function(cmd, msg = "") {
  assert_that(is.string(cmd))
  unname(Sys.which(cmd) != "")
}
on_failure(has_command) <- function(call, env) {
  paste0("Dependency ", sQuote(eval(call$cmd, env)), " is not installed\n",
         eval(call$msg, env))
}

## UnitTests: inst/tests/test-utils.r
##
#' Format paragraphs
#' 
#' Similar to \code{\link{strwrap}} but returns a single string with
#' linefeeds inserted
#' 
#' @param s a character vector.
#' @param width a positive integer giving the column for inserting
#' linefeeds
#' @param indent an integer giving the indentation of the first line of
#' the paragraph; negative values of \code{indent} are allowed and reduce
#' the width for the first line by that value.
#' @param offset a non-negative integer giving the indentation of all
#' but the first line.
#' @param split regular expression used for splitting. Defaults to whitespace.
#' @param FORCE Words are force-split if the available width is too small.
#' @param FULL_FORCE Lines are split exactly at the specified width
#' irrespective of whether there is whitespace or not.
#' 
#' @return a character vector
#' @keywords internal
linebreak <- function(s, width = getOption("width") - 2,
                      indent = 0, offset = 0, split = ' ',
                      FORCE = FALSE, FULL_FORCE = FALSE) {
  assert_that(offset >= 0)
  first_pass <- TRUE
  s <- as.character(s)
  if (length(s) == 0) return("")
  
  (function(s) {
    # remove leading and trailing blanks
    # convert newlines, tabs, spaces to " "
    # find first position where 'split' applies
    indent_string <- dup(' ', indent)
    indent <- abs(indent)
    offset_string <- paste0("\n", dup(' ', offset))
    if (!FULL_FORCE) {
      s <- gsub("\\s+", " ", trim(s), perl = TRUE)
    }
    fws <- regexpr(split, s, perl = TRUE)
    if (first_pass) {
      string_width <- indent + nchar(s)
      .offset <- 0
    } else {
      string_width <- offset + nchar(s)
      .offset <- offset
    }
    if (string_width > width) {
      # if not everything fits on one line
      if (FULL_FORCE || (FORCE && (fws == -1 || fws >= width - .offset - indent))) {
        # if no whitespace or first word too long and force break cut through the
        # middle of a word
        pat1 <- paste0("^.{", width - .offset - indent, "}(?=.+)")
        pat2 <- paste0("(?<=^.{", width - .offset - indent, "}).+")
        leading_string <- regmatches(s, regexpr(pat1, s, perl = TRUE))
        trailing_string <- regmatches(s, regexpr(pat2, s, perl = TRUE)) 
      } else if (!FORCE && (fws == -1 || fws >= (width - .offset + indent))) {
        # if no whitespace or first word too long and NO force break stop right here
        stop("Can't break in the middle of a word. Use the force!")
      } else {
        # break the line
        s_split <- usplit(s, split)
        s_cum   <- cumsum(nchar(s_split) + nchar(split))
        leading_string <- 
          paste0(s_split[s_cum < width - .offset - indent + 1],
                 ifelse(split == " ", "", split), collapse = split)
        trailing_string <- 
          paste0(s_split[s_cum >= width - .offset - indent + 1], collapse = split)
      }
      first_pass <<- FALSE
      indent <<- 0
      s <- paste0(indent_string, leading_string, offset_string, Recall(trailing_string))
    } else {
      # if everything fits on one line go with the string + indent
      paste0(indent_string, s)
    }
  })(s)
}

strsplitN <- function(x, split, n, from = "start", collapse = split, ...) {
  from <- match.arg(from, c("start", "end"))
  xs <- strsplit(x, split, ...)
  end <- vapply(xs, length, 0L)
  if (from == "end") {
    end <- end + 1L
    n <- lapply(end, `-`, n)
    n <- .mapply(`[<-`, list(x = n, i = lapply(n, `<`, 0), value = 0L), NULL)
  } else {
    n <- lapply(rep(0, length(xs)), `+`, n)
    n <- .mapply(`[<-`, list(x = n, i = Map(`>`, n, end), value = end), NULL)
  }  
  n <- lapply(n, Compose("sort", "unique"))
  unlist(.mapply(function(x, n) paste0(x[n], collapse = collapse), list(x = xs, n = n), NULL))
}

replace_ext <- function(file, replacement) {
  split_file <- usplit(file, '\\.')
  split_file[length(split_file)] <- replacement
  paste0(split_file, collapse = ".")
}

recycle <- function(val, len) {
  lv <- length(val)
  if (len > lv) {
    val <- c(rep(val, len%/%lv), val[seq_len(len%%lv)])
  }
  val
}

is_compound <- function(x) {
  if (is(x, "gbFeatureTable")) {
    return(vapply(x, function(f) !is.na(f@location@compound), FALSE))
  } else if (is(x, "gbFeature")) {
    return(!is.na(x@location@compound))
  } else if (is(x, "gbLocation")) {
    return(!is.na(x@compound))
  }
}

get_compounds <- function(x) {
  x <- x[which(is_compound(x))]
  if (length(x) == 0) return(NA_real_) 
  cL <- vapply(x, function(f) nrow(f@location@range), 0L)
  cL
}

.access <- function(x, which, dbxrefs, use.names = TRUE) {
  q <- x@qualifiers
  n <- length(q)
  els <- c(which[which != 'db_xref' & which != '\\bdb_xref\\b'], dbxrefs)
  if (n == 0) {
    ans <- rep(NA_character_, length(els))
    if (use.names) {
      return(setNames(ans, trim(els, "\\\\b")))
    } else {
      return(ans)
    }
  }
  if (length(which) == 1) {
    idx <- grepl(which, names(q))
    if (any(idx)) {
      ans <- q[idx]
    } else {
      els <- els %||% which
      ans <- setNames(rep(NA_character_, length(els)), trim(which, "\\\\b"))
    }
  } else {
    ans <- lapply(lapply(which, grepl, names(q)), function(i) q[i])
    na  <- which(vapply(ans, length, 0) == 0)
    if (length(na) > 0) {
      for (i in na) {
        ans[[i]] <- setNames(NA_character_, nm = trim(which[i], "\\\\b"))
      }
    }
  }
  ans <- unlist(ans)
  ## with db_xref remove the db identifier from the id number and
  ## attach it to the names
  if (any(dbx_idx <- names(ans) == 'db_xref')) {
    dbx <- strsplit(ans[dbx_idx], ':')
    dbx_dbnm <- vapply(dbx, `[`, 1L, FUN.VALUE = "", USE.NAMES = FALSE)
    dbx_dbid <- vapply(dbx, `[`, 2L, FUN.VALUE = "", USE.NAMES = FALSE)
    if (length(dbxrefs) > 0) {
      ans <- c(ans[!dbx_idx], setNames(dbx_dbid[match(dbxrefs, dbx_dbnm)],
                                       paste0('db_xref.', dbxrefs)))
    } else {
      nm <- if (all(is.na(dbx_dbnm))) "" else paste0('db_xref.', dbx_dbnm)
      ans <- c(ans[!dbx_idx], setNames(dbx_dbid, nm))
    }
  }
  return(if (use.names) ans else unname(ans))
}

#' @return A named character vector or a list of named character vectors
#' of the qualifiers specified in \code{which}
.qual_access <- function(x, which = "", fixed = FALSE, use.names = TRUE) {
  dbxrefs <- NULL
  dbx <- grepl('db_xref[:.].+', which)
  if (any(dbx)) {
    dbxrefs <- strsplitN(which[dbx], '[:.]', 2L)
    which <- c(which[!dbx], 'db_xref')
  }
  if (fixed) {
    which <- wrap(which, "\\b")
  }
  if (is(x, "gbFeature")) {
    .access(x, which, dbxrefs, use.names)
  } else if (is(x, "gbFeatureTable")) {
    lapply(x, .access, which, dbxrefs, use.names)
  }
}

.simplify <- function(x, unlist = TRUE) {
  # case 1) x is atomic OR a list of single, equally named, elements 
  #         if unlist == TRUE        => return as is
  #         if unlist == FALSE       => return data.frame
  # case 2) x is a list of equally long (equally named) elements
  #                                  => return data.frame
  # case 3) x is a list of unequally long, unnamed elements
  #                                  => return as is
  # case 4) x is a list of unequally named elements
  #                                  => return expanded data.frame
  len <- unique(unlist(lapply(x, length)))
  nm  <- compactChar(unique(unlist(lapply(x, names))))
  ## case 1
  if (is.atomic(x) || len == 1L && length(nm) <= 1L) {
    if (unlist) {
      unlist(x, recursive = FALSE, use.names = FALSE)
    } else {
      nm <- nm %||% 'X1'
      x <- unlist(x, recursive = FALSE, use.names = FALSE)
      data.frame(stringsAsFactors = FALSE,
                 matrix(x, ncol = 1, dimnames = list(NULL, nm)))
    }
  }
  ## case 2
  else if (length(len) == 1L &&
             (is.null(nm) || length(nm) == len)) {
    r <- as.vector(unlist(x, recursive = FALSE, use.names = FALSE))
    data.frame(stringsAsFactors = FALSE,
               matrix(r, ncol = len, byrow = TRUE,
                      dimnames = list(NULL, nm)))
  }
  # case 3
  else if (is.null(nm)) {
    x
  }
  # case 4
  else if ((lnm <- length(nm)) != min(len)) {
    m <- matrix(rep(NA, length(x) * lnm), ncol = lnm, dimnames = list(NULL, nm))
    for (i in seq_along(x)) {
      for (n in nm) {
        m[i, n] <- unname(x[[i]][n])
      }
    }
    data.frame(stringsAsFactors = FALSE, m)
  }
}


.seq_access <- function(x) {
  seq <- .sequence(x)
  if (length(seq) == 0) {
    return(seq)
  }
  SEQFUN <- match.fun(class(seq))
  if (is(x, "gbFeature")) {
    seq <- merge_seq(seq, x, SEQFUN)
  } else if (is(x, "gbFeatureTable")) {
    seq <- Reduce(append, lapply(x, merge_seq, seq = seq, SEQFUN = SEQFUN))
  }
  seq
}

# merge Sequences
#' @importFrom Biostrings xscat
merge_seq <- function(seq, x, SEQFUN) {
  if (length(start(x)) == 1L) {
    outseq <- XVector::subseq(x = seq, start = start(x), end = end(x))
  } else {
    outseq <- do.call(xscat, Map(subseq, x = seq, start = start(x), end = end(x)))
  }
  outseq <- SEQFUN(outseq)
  outseq@ranges@NAMES <- .defline(x)
  outseq
}

parse_dbsource <- function(dbsource) {
  if (is.na(dbsource)) {
    '|gb|'
  }
  else {
    db <- strsplitN(dbsource, ": | ", 1L)
    db <- switch(db, accession = 'gb', REFSEQ = 'ref', db)
    paste0('|', db, '|')
  }
}
