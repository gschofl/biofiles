#' @importFrom Biostrings DNAStringSet
#' @importFrom IRanges new2
#' @importFrom parallel mclapply mcmapply detectCores
#' @importFrom foreach foreach registerDoSEQ "%dopar%"
#' @importFrom iterators iter
NULL

## declare "rec" global so that "codetools" don't complain
## see "http://r.789695.n4.nabble.com/R-CMD-check-and-foreach-code-td4687660.html"
globalVariables("rec")

parse_record <- function(rcd) {
  ## check that the gbk/embl record is not empty
  if (length(rcd) == 0L) {
    stop("This \"gbRecord\" is empty.")
  }
  ## check if gbk/embl record contains multiple entries
  end_of_record <- grep('^//$', rcd)
  n_records <- length(end_of_record)
  if (n_records > 1L) {
    start_of_record <- c(1, end_of_record[-n_records] + 1)
    irec <- iter(ixsplit(rcd, start_of_record))
  } else {
    registerDoSEQ()
    irec <- iter(list(rcd))
  }
  rcd_list <- if (is.gbk(rcd)) {
    foreach(rec = irec, .inorder = FALSE) %dopar% parse_gbk_record(rec)
  } else if (is.embl(rcd)) {
    foreach(rec = irec, .inorder = FALSE) %dopar% parse_embl_record(rec)
  } else {
    stop("Unknown format <", rcd[1], ">")
  }
  if (length(rcd_list) == 1L) {
    rcd_list[[1L]]
  } else {
    gbRecordList(rcd_list)
  }
}
