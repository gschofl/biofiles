.GBFIELDS <- c("@G@I", "accession", "comment", "date", "dblink",
               "dbsource", "definition", "division", "features",
               "keywords", "length", "lineage", "locus", "organism",
               "references", "sequence", "source", "topology",
               "type", "version")

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.biofiles <- list(
    ## switch off tests that don't work with R CMD check
    biofiles.test.parser = FALSE
  )
  toset <- !(names(op.biofiles) %in% names(op))
  if (any(toset)) {
    options(op.biofiles[toset])
  }
  invisible()
}
