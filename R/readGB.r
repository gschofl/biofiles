#' General function for importing GenBank flat files
#'
#' For a description of the GenBank format see
#' \url{http://www.ncbi.nlm.nih.gov/collab/FT/}
#'
#' @usage readGB(gb, with_sequence=TRUE, force=FALSE)
#'
#' @param gb Path to a GenBank flat file or an 
#' \code{\link[rentrez]{efetch-class}} object containing GenBank record(s).
#' @param with_sequence If \code{TRUE}, sequence information
#' will be included.
#' @param force If \code{TRUE} existing database directories are
#' overwritten without prompting.
#' 
#' @return A (list of) \code{\link{gbRecord-class}} object(s).
#' 
#' @export
readGB <- function (gb, with_sequence = TRUE, force = FALSE) {
  # we store efetch data in temporary files
  if (is(gb, "efetch")) {
    ## we can parse rettype = gbwithparts, gb, gp and retmode =  text
    if (!grepl("^gb|^gp", gb@type) || gb@mode != "text")
      stop("Must use efetch with rettype='gbwithparts','gb', or 'gp' and retmode='text'")
    
    split_gb <- strsplit(gb@content, "\n\n")[[1L]]
    n <- length(split_gb)
    db_path <- replicate(n, tempfile(fileext=".db"))
    
    for (i in seq_len(n)) {
      gb_data <- strsplit(split_gb[i], "\n")[[1L]]
      cat(gettextf("Importing into %s\n", dQuote(basename(db_path[i]))))
      .parseGB(gb_data, db_path[i], with_sequence=with_sequence, force=force)
    }
    
  } else if (!isS4(gb) && file.exists(gb)) {
    con <- file(gb, open="rt")
    on.exit(close(con))
    db_path <- paste0(gb, ".db")
    .parseGB(gb_data=readLines(con), db_path, with_sequence=with_sequence,
             force=force)
    
  } else {
    stop("'gb' must be a valid GenBank flat file or an 'efetch' object containing GenBank records")
  }
  return(invisible(db_path))
}


# --R-- vim:ft=r:sw=2:sts=2:ts=4:tw=76:
#       vim:fdm=marker:fmr={{{,}}}:fdl=0
