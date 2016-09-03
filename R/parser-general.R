#' @importFrom foreach foreach registerDoSEQ "%dopar%"
#' @importFrom iterators iter
NULL

## declare "rec" global so that "codetools" don't complain
## see "http://r.789695.n4.nabble.com/R-CMD-check-and-foreach-code-td4687660.html"
globalVariables("rec")

parse_record <- function(rcd, progress = FALSE) {
  ## check that the gbk/embl record is not empty
  if (length(rcd) == 0L) {
    stop("This \"gbRecord\" is empty.")
  }
  ## check if gbk/embl record contains multiple entries
  end_of_record <- grep('^//$', rcd)
  n_records <- length(end_of_record)
  pb <- make_progress_bar(n_records)
  if (n_records > 1L) {
    start_of_record <- c(1, end_of_record[-n_records] + 1)
    irec <- iter(ixsplit(rcd, start_of_record))
    #rec <- nextElem(irec)
  } else {
    registerDoSEQ()
    irec <- iter(list(rcd))
  }
  rcd_list <- if (is.gbk(rcd)) {
    foreach(rec = irec, .inorder = FALSE) %dopar% {
      if (progress) pb()
      gbk_record(rec)
    }
  } else if (is.embl(rcd)) {
    foreach(rec = irec, .inorder = FALSE) %dopar% {
      if (progress) pb()
      embl_record(rec)
    }
  } else {
    stop("Unknown format <", rcd[1], ">")
  }
  if (length(rcd_list) == 1L) {
    rcd_list[[1L]]
  } else {
    gbRecordList(rcd_list)
  }
}

make_progress_bar <- function(n) {
  now <- function() {
    proc.time()[[3]]
  }
  
  tick <- function() {
    if (i < n) {
      i <<- i + 1L
    }
  }
  
  repstr <- function(x, i) {
    paste(rep.int(x, i), collapse = "")
  }
  
  width <- function() {
    getOption("width") - nchar("|[100%] - 100h to go.") - 2
  }
  
  cat_line <- function(...) {
    msg <- paste0(..., collapse = "")
    gap <- max(c(0, getOption("width") - nchar(msg, "width")))
    cat("\r", msg, rep.int(" ", gap), sep = "")
    flush.console()
  }
  
  show_time <- function(x) {
    if (x < 60) {
      paste0(round(x), "s")
    } else if (x < 60*60) {
      paste0(round(x/60), "m")
    } else {
      paste0(round(x/(60*60)), "h")
    }
  }
  
  i <- 0L
  init_time <- now()
  
  function() {
    tick()
    avg <- (now() - init_time) / i
    time_left <- (n - i)*avg
    nbars <- trunc(i/n*width())
    cat_line("|", repstr("=", nbars), repstr(" ", width() - nbars), "| [",
             format(round(i/n*100), width = 3), "%] - ", show_time(time_left),
             " to go.")
  }
}

## characters permitted to occur in feature table component names:
## [A-Z], [a-z], [0-9], [_'*-] --> [[:alnum]_'*\\-]
##
## Qualifiers:
## A qualifier takes the form of a slash (/) followed by its name,
## if applicable, an equals sign ( = ) and its value.
## 
## Qualifier value formats:
## Free text; Controlled vocabulary or enumerated values; Citation or 
## reference numbers
##
## Free text: "text"
##   enclosed in double quotation marks; composed of printable
##   characters (ASCII values 32-126 decimal); internal quotation markes are
##   ""escaped"" by placing a second quotation mark immediately before.
##
## Controlled vocabulary/Enumerators:
##   /anticodon = (pos:<base_range>,aa:<amino_acid>)
##   /codon_start = 1, 2, or 3
##   /direction = left, right, or both
##   /estimated_length = unknown or <integer>
##   /calculated_mol_wt = <integer>
##   /mod_base = name of modified base
##   /number = unquoted text
##   /rpt_type = <repeat_type>
##   /rpt_unit_range = <base_range>
##   /tag_peptide = <base_range>
##   /transl_except = (pos:location,aa:<amino_acid>)
##   /transl_table = <integer>
##  
## Citations:
##   /citation = [number] e.g. /citation = [3]
##   /compare = [accession-number.sequence-version] e.g. /compare = AJ634337.1
#' @keywords internal
parse_features <- function(x, seqinfo) {
  feature_start <- which(substr(x, 6, 6) != " ")
  fl <- ixsplit(x, feature_start)
  mc_cores <- floor(parallel::detectCores()*0.75)
  id <- seq_along(feature_start)
  ftbl <- parallel::mcmapply(
    gbFeature, feature = fl, id = id,
    MoreArgs = list(accession = getAccession(seqinfo)[1]),
    SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.cores = mc_cores
  )                   
  S4Vectors::new2('gbFeatureTable', .Data = ftbl, .id = id, .seqinfo = seqinfo, check = FALSE) 
}

#' @keywords internal
parse_sequence <- function(seq, acc, seqtype, src) {
  # read.BStringSet() does not support connections and
  # currently only accepts fasta format. So we write out gb_sequence as
  # a temporary fasta file and read it back in as an AAStringSet or
  # DNAStringSet (mRNA etc seems to be encoded with Ts rather then Us,
  # so we use DNAStringSets for RNA)
  if (is.null(seq)) {
    return(Biostrings::BStringSet())
  } else {
    tmp <- tempfile()
    on.exit(unlink(tmp))
    writeLines(text = join_seq(seq, acc, src), tmp)
    origin <- switch(seqtype,
                     AA = Biostrings::readAAStringSet(tmp, format = "fasta"),
                     Biostrings::readDNAStringSet(tmp, format = "fasta"))
    origin
  }
}

join_seq <- function(seq, acc, src = c("gbk", "embl")) {
  src <- match.arg(src, c("gbk", "embl"))
  .substr <- switch(src,
                    gbk  = Partial('substr', start = 11, stop = 75),
                    embl = Partial('substr', start = 6, stop = 70)
  )
  mc_cores <- floor(parallel::detectCores()*0.75)
  s <- unlist(parallel::mclapply(seq, function(x) {
    paste0(strsplit(.substr(x), ' ')[[1L]], collapse = '')
  }, mc.cores = mc_cores))
  s <- c(paste0(">", acc), s)
  s
}

