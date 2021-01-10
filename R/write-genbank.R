#' @export
#' @rdname write.GenBank-methods
setMethod("write.GenBank", "gbRecord", 
          function (x, file, header = TRUE, sequence = TRUE, append = FALSE) {
            .writeGenBank(x = x, file = file, header = header, sequence = sequence,
                          append = append) 
          })

#' @export
#' @rdname write.GenBank-methods
setMethod("write.GenBank", "gbFeatureTable", 
          function (x, file, header = TRUE, sequence = TRUE, append = FALSE) {
            x <- as(x, "gbRecord")
            .writeGenBank(x = x, file = file, header = header, sequence = sequence,
                          append = append) 
          })

.writeGenBank <- function(x, file, header = TRUE, sequence = TRUE, append = FALSE) {
  op <- options(useFancyQuotes = FALSE)
  on.exit(options(op))
  
  if (header) {
    header(x)$write(file = file, append = append, sep = "")
  }
  
  cat("FEATURES:            Location/Qualifiers:\n", file = file, append = TRUE)
  f <- lapply(.features(x), show_gbFeature, write_to_file = TRUE)
  cat(paste0(f, collapse = "\n"), file = file, append = TRUE)

  if (sequence) {
    .writeSequence(x, file)
  }
  invisible()
}

.writeSequence <- function (x, file = "out.gbk") {
  if (length(seq <- getSequence(x)) > 0L) {
    lineno <- seq(from = 1, to = seq@ranges@width, by = 60)
    lines <- seq_along(lineno)
    n_lines <- length(lines)
    s <- character(n_lines)
    for (i in lines) {
      seqw <- ifelse(i <  n_lines, i*60, seq@ranges@width)
      seqs <- XVector::toString(XVector::subseq(seq, 1 + (i - 1)*60, seqw))
      nnn <- seq(1, nnncc <- nchar(seqs), by = 10) 
      s[i] <- paste0(substring(seqs, nnn, c(nnn[-1]-1, nnncc)), collapse = " ")
    }
    s <- sprintf("%+9s %s", lineno, s)
    cat("\nORIGIN", file = file, sep = "\n", append = TRUE)
    cat(s, file = file, sep = "\n", append = TRUE)
    cat("//", file = file, append = TRUE)
  } else {
    cat("\n//", file = file, append = TRUE)
  }
  invisible()
}

#' @export
#' @rdname saveRecord-methods
setMethod("saveRecord", "gbRecord", function(x, file = NULL, dir = ".", ...) {
  if (!is.character(file)) {
    fname <- paste0(getAccession(x), '.rds')
    file <- normalizePath(file.path(dir, fname), mustWork = FALSE)
  } else {
    file <- normalizePath(file.path(dir, file), mustWork = FALSE)
  }
  saveRDS(object = x, file = file, ...)
  return(invisible())
})

#' @export
#' @rdname saveRecord-methods
setMethod("saveRecord", "gbRecordList", function(x, file = NULL, dir = ".", ...) {
  if(!is.character(file)) {
    fname <- paste0(ellipsize(collapse(getAccession(x), '_'), width = 60, ellipsis = '__'), '.rds')
    file <- normalizePath(file.path(dir, fname), mustWork = FALSE)
  } else {
    file <- normalizePath(file.path(dir, file), mustWork = FALSE)
  }
  saveRDS(object = x, file = file, ...)
  return(invisible())
})


#' @export
#' @rdname saveRecord-methods
loadRecord <- function(file, ...) {
  if (missing(file)) {
    stop("No filename provided")
  }
  readRDS(file = file, ...)
}




