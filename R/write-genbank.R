
setMethod("write.GenBank", "gbRecord", 
          function (x, file, header = TRUE, sequence = TRUE, append = FALSE) {
            # write header
            h <- if (header) .writeHeader(x, file) else ""
            # write features
            op <- options(useFancyQuotes=FALSE)
            #cat("Writing features\n")
            f <- unlist(lapply(getFeatures(x), .writeFeature))
            cat(paste0(f, collapse="\n"), file=file, append=TRUE)
            options(op)
            # write origin
            s <- if (sequence) .writeSequence(x, file) else "//"
            invisible(list(header=h, features=f, sequence=s)) 
          })



.writeHeader <- function (x, outfile = "out.gbk") {
  type <- if (getMoltype(x) == "DNA") "bp" else "  "
  loc_line <- sprintf("%-12s%-17s %+10s %s    %-6s  %-8s %s %s",
                      "LOCUS", getLocus(x), getLength(x), type, getMoltype(x),
                      getTopology(x), getDivision(x), base::toupper(format(getDate(x)[2], "%d-%b-%Y")))
  def_line <- sprintf("%-12s%s", "DEFINITION", linebreak(getDefinition(x), width=79, offset=12))
  acc_line <- sprintf("%-12s%s", "ACCESSION", getAccession(x))
  ver_line <- sprintf("%-12s%-12s%s%s", "VERSION", getVersion(x), "GI:", getGeneID(x))
  dbl_line <- if (!is.null(dblink <- getDBLink(x))) {
    sprintf("%-12s%s%s", "DBLINK", "Project: ", dblink)
  } else {
    character(0)
  }
  key_line <- sprintf("%-12s%s", "KEYWORDS", linebreak(getKeywords(x), width=79, offset=12))
  src_line <- sprintf("%-12s%s", "SOURCE", linebreak(getSource(x), width=79, offset=12))
  org_line <- sprintf("%-12s%s", "  ORGANISM", paste0(getOrganism(x),
                             linebreak(getTaxonomy(x), width=79, indent=12, offset=12),
                             sep="\n"))
  ref_line <- sprintf("%-12s%-3s(bases %s to %s)\n%-12s%s\n%-12s%s\n%-12s%s",
                      "REFERENCE", 1, 1, getLength(x),
                      "  AUTHORS", "authors",
                      "  TITLE", "title",
                      "  JOURNAL", "journal")
  com_line <- if (!is.null(comment <- getComment(x))) {
    sprintf("%-12s%s", "COMMENTS", linebreak(comment, width=79, offset=12))
  } else {
    character(0)
  }
  f_line <- sprintf("%-21s%s", "FEATURES", "Location/Qualifiers\n")
  header <- paste0(loc_line, def_line, acc_line, ver_line, dbl_line,
                   key_line, src_line, org_line, ref_line, com_line,
                   f_line, collapse="\n")
  header <- base::gsub("\n{2,}", "\n", header)
  cat(header, file=outfile)
  invisible(header)
}



.writeFeature <- function (f) {
  loc_line <- sprintf("%s%-16s%s",
                      dup(' ', 5),
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



.writeSequence <- function (x, outfile = "out.gbk") {
  
  if (exists("sequence", envir=x@seqinfo)) {
    seq <- getSequence(x)
    lineno <- seq(from=1, to=seq@ranges@width, by=60)
    lines <- seq_along(lineno)
    n_lines <- length(lines)
    s <- character(n_lines)
    
    for (i in lines) {
      seqw <- ifelse(i <  n_lines, i*60, seq@ranges@width)
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

