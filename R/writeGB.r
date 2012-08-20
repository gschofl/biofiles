##' General function for writing out GenBank flat files
##'
##' For a description of the GenBank format see
##' \url{http://www.ncbi.nlm.nih.gov/collab/FT/}
##'
##' @usage writeGB(db, outfile = "out.gbk")
##'
##' @param db A \code{gbRecord} object.
##' @param outfile Output file.
##' 
##' @export
writeGB <- function(db, outfile = "out.gbk") {
  
  if (file.exists(outfile)) {
    unlink(outfile)
  }
  
  ## write header
  h <- .writeHeader(db, outfile)
  
  ## write features
  op <- options("useFancyQuotes")
  options(useFancyQuotes=FALSE)
  cat("Writing features\n")
  f <- unlist(lapply(db$features, .writeFeature))
  cat(paste(f, collapse="\n"), file=outfile, append=TRUE)
  options(op)
  
  ## write origin
  cat("Writing sequence\n")
  s <- .writeSequence(db, outfile)
  
  invisible(list(header=h, features=f, sequence=s))
  
}


.writeHeader <- function (db, outfile = "out.gbk") {
  
  type <- if (db$type == "DNA") "bp" else "  "
  loc_line <- sprintf("%-12s%-17s %+10s %s    %-6s  %-8s %s %s",
                      "LOCUS", db$locus, db$length, type, db$type, db$topology,
                      db$division, toupper(format(db$date, "%d-%b-%Y")))
  def_line <- sprintf("%-12s%s", "DEFINITION", 
                      linebreak(db$definition, width=79, offset=12))
  acc_line <- sprintf("%-12s%s", "ACCESSION", db$accession)
  ver_line <- sprintf("%-12s%-12s%s%s", "VERSION", db$version, "GI:", db$GI)
  
  dbl_line <- if (!is.null(db$dblink)) {
    sprintf("%-12s%s%s", "DBLINK", "Project: ", db$dblink)
  } else {
    character()
  }
  
  key_line <- sprintf("%-12s%s", "KEYWORDS",
                      linebreak(db$keywords, width=79, offset=12))
  src_line <- sprintf("%-12s%s", "SOURCE",
                      linebreak(db$source, width=79, offset=12))
  org_line <- sprintf("%-12s%s", "  ORGANISM",
                      paste(db$organism,
                            linebreak(db$lineage, width=79, indent=12, offset=12),
                            sep="\n"))
  ref_line <- sprintf("%-12s%-3s(bases %s to %s)\n%-12s%s\n%-12s%s\n%-12s%s",
                      "REFERENCE", 1, 1, db$length,
                      "  AUTHORS", "authors",
                      "  TITLE", "title",
                      "  JOURNAL", "journal")
  
  com_line <- if (!is.null(db$comment)) {
    sprintf("%-12s%s", "COMMENTS", linebreak(db$comment, width=79, offset=12))
  } else {
    character()
  }
  
  f_line <- sprintf("%-21s%s", "FEATURES", "Location/Qualifiers\n")
  
  header <- paste(loc_line, def_line, acc_line, ver_line, dbl_line,
                  key_line, src_line, org_line, ref_line, com_line,
                  f_line, sep="\n")
  header <- gsub("\n{2,}", "\n", header)
  cat(header, file=outfile)
  
  invisible(header)
}


.writeFeature <- function (f) {
  
  
  loc_line <- sprintf("%s%-16s%s",
                      blanks(5),
                      f@key,
                      linebreak(as(f@location, "character"),
                                width=79, offset=21, indent=0, split=","))
  qua <- names(f@qualifiers)
  val <- linebreak(dQuote(f@qualifiers), width=79, offset=21, FORCE=TRUE,
                   indent=-(nchar(qua) + 2))
  qua_line <- sprintf("%+22s%s=%s", "/", qua, val)
  
  feature <- paste0(loc_line, "\n", paste0(qua_line, collapse="\n"))
  feature 
}


.writeSequence <- function (db, outfile = "out.gbk") {
  
  if (!is.null(db$sequence)) {
    sequence <- db$sequence
    lineno <- seq(from=1, to=sequence@ranges@width, by=60)
    lines <- seq_along(lineno)
    n_lines <- length(lines)
    s <- character(n_lines)
    
    for (i in lines) {
      seqw <- ifelse(i <  n_lines, i*60, sequence@ranges@width)
      seq <- toString(subseq(sequence, 1 + (i - 1)*60, seqw))
      s[i] <- paste0(strsplit(seq, "(?<=.{10})(?=.)", perl=TRUE)[[1]], collapse=" ")     
    }
    
    s <- sprintf("%+9s %s", lineno, s)
    
    cat("\nORIGIN", file=outfile, sep="\n", append=TRUE)
    cat(s, file=outfile, sep="\n", append=TRUE)
    cat("//", file=outfile, append=TRUE)
  } else {
    cat("//", file=outfile, append=TRUE)
  }
  
  invisible(s)
}

