#' @autoImports
setMethod("write.GenBank", "gbRecord", 
          function (x, file, header = TRUE, append = FALSE) {
            # write header
            if (header) {
              h <- .writeHeader(x, file)
            } else {
              h <- ""
            }
            # write features
            op <- options("useFancyQuotes")
            options(useFancyQuotes=FALSE)
            #cat("Writing features\n")
            f <- unlist(lapply(getFeatures(x), .writeFeature))
            cat(paste0(f, collapse="\n"), file=file, append=TRUE)
            options(op)
            # write origin
            #cat("Writing sequence\n")
            s <- .writeSequence(x, file)
            
            invisible(list(header=h, features=f, sequence=s)) 
          })


#' @autoImports
.writeHeader <- function (x, outfile = "out.gbk") {
  
  type <- if (x@type == "DNA") "bp" else "  "
  loc_line <- sprintf("%-12s%-17s %+10s %s    %-6s  %-8s %s %s",
                      "LOCUS", getLocus(x), getLength(x), type, x@moltype, x@topology,
                      x@division, base::toupper(format(x@date, "%d-%b-%Y")))
  def_line <- sprintf("%-12s%s", "DEFINITION", 
                      linebreak(definition(x), width=79, offset=12))
  acc_line <- sprintf("%-12s%s", "ACCESSION", getAccession(x))
  ver_line <- sprintf("%-12s%-12s%s%s", "VERSION", x@version, "GI:", x@GI)
  
  dbl_line <- if (!is.null(x@dblink)) {
    sprintf("%-12s%s%s", "DBLINK", "Project: ", x@dblink)
  } else {
    character()
  }
  
  key_line <- sprintf("%-12s%s", "KEYWORDS",
                      linebreak(x@keywords, width=79, offset=12))
  src_line <- sprintf("%-12s%s", "SOURCE",
                      linebreak(x@source, width=79, offset=12))
  org_line <- sprintf("%-12s%s", "  ORGANISM",
                      paste0(x@organism,
                             linebreak(x@lineage, width=79, indent=12, offset=12),
                             sep="\n"))
  ref_line <- sprintf("%-12s%-3s(bases %s to %s)\n%-12s%s\n%-12s%s\n%-12s%s",
                      "REFERENCE", 1, 1, getLength(x),
                      "  AUTHORS", "authors",
                      "  TITLE", "title",
                      "  JOURNAL", "journal")
  
  com_line <- if (!is.null(x@comment)) {
    sprintf("%-12s%s", "COMMENTS", linebreak(x@comment, width=79, offset=12))
  } else {
    character()
  }
  
  f_line <- sprintf("%-21s%s", "FEATURES", "Location/Qualifiers\n")
  
  header <- paste0(loc_line, def_line, acc_line, ver_line, dbl_line,
                   key_line, src_line, org_line, ref_line, com_line,
                   f_line, sep="\n")
  header <- base::gsub("\n{2,}", "\n", header)
  cat(header, file=outfile)
  
  invisible(header)
}


#' @autoImports
.writeFeature <- function (f) {
  loc_line <- sprintf("%s%-16s%s",
                      blanks(5),
                      f@key,
                      linebreak(as(f@location, "character"),
                                width=79, offset=21, indent=0, split=","))
  qua <- names(f@qualifiers)
  val <- linebreak(dQuote(f@qualifiers), width=79, offset=21, FORCE=TRUE,
                   indent=-(base::nchar(qua) + 2))
  qua_line <- sprintf("%+22s%s=%s", "/", qua, val)
  
  feature <- paste0(loc_line, "\n", paste0(qua_line, collapse="\n"))
  feature 
}


#' @autoImports
.writeSequence <- function (x, outfile = "out.gbk") {
  
  if (exists("sequence", envir=x@seqinfo)) {
    seq <- getSequence(x)
    lineno <- seq(from=1, to=seq@ranges@width, by=60)
    lines <- seq_along(lineno)
    n_lines <- length(lines)
    s <- character(n_lines)
    
    for (i in lines) {
      seqw <- base::ifelse(i <  n_lines, i*60, seq@ranges@width)
      seqs <- XVector::toString(XVector::subseq(seq, 1 + (i - 1)*60, seqw))
      s[i] <- paste0(strsplit(seqs, "(?<=.{10})(?=.)", perl=TRUE)[[1]], collapse=" ")     
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

