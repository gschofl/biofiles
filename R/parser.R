#' @importFrom Biostrings DNAStringSet
#' @importFrom IRanges new2
#' @importFrom parallel mclapply mcmapply detectCores
NULL

parse_gb_record <- function(gb_data, with_sequence = TRUE) {
  # get a vector with the positions of the main GenBank fields
  gbf <- grep("^[A-Z//]+", gb_data)
  gbf_names <- strsplitN(gb_data[gbf], split=" +", 1)
  names(gbf) <- gbf_names
  gb_contig <- gb_sequence <-  NULL
  
  # Check the presence of mandatory fields
  essential <- c("LOCUS", "DEFINITION", "ACCESSION", "VERSION", "FEATURES")
  if (any(is.na(charmatch(essential, gbf_names)))) {
    stop("Some fields seem to be missing from the GenBank file")
  }
  
  ## HEADER
  seqenv <- seqinfo(gbHeader(gb_data[seq.int(gbf["FEATURES"]) - 1]), NULL)
    
  ## ORIGIN
  if (!is.na(origin <- gbf["ORIGIN"])) {
    seq_idx <- seq.int(origin + 1, gbf["//"] - 1)
    if (length(seq_idx) > 1L && seq_idx[2] < seq_idx[1]) {
      # happens if "//" is next after "ORIGIN", i.e. no  sequence is present
      gb_sequence <- NULL
    } else {
      gb_sequence <- gb_data[seq_idx]
    }
  } else if (!is.na(contig <- gbf["CONTIG"])) {
    gb_contig <- gbLocation(strsplitN(gb_data[contig], "CONTIG", 2, fixed=TRUE))
  }
  seqenv$sequence <- gbSequence(gb_sequence, getAccession(seqenv), getMoltype(seqenv))
  
  ## FEATURES
  if (match("FEATURES", gbf_names) == length(gbf)) {
    gb_features <- gb_data[seq.int(gbf["FEATURES"] + 1, length(gb_data) - 1)]
  } else {
    gb_features <- gb_data[seq.int(gbf["FEATURES"] + 1, gbf[match("FEATURES", gbf_names) + 1] - 1)]
  }
  if (length(gb_features) < 2L) {
    stop("No features in the GenBank file")
  }
  features <- gbFeatures(gb_features, seqinfo=seqenv)
  
  new_gbRecord(seqinfo=seqenv, features=features, contig=gb_contig)
}

## characters permitted to occur in feature table component names:
## [A-Z], [a-z], [0-9], [_'*-] --> [[:alnum]_'*\\-]
##
## Qualifiers:
## A qualifier takes the form of a slash (/) followed by its name,
## if applicable, an equals sign (=) and its value.
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
##   /anticodon=(pos:<base_range>,aa:<amino_acid>)
##   /codon_start=1, 2, or 3
##   /direction=left, right, or both
##   /estimated_length=unknown or <integer>
##   /calculated_mol_wt=<integer>
##   /mod_base=name of modified base
##   /number=unquoted text
##   /rpt_type=<repeat_type>
##   /rpt_unit_range=<base_range>
##   /tag_peptide=<base_range>
##   /transl_except=(pos:location,aa:<amino_acid>)
##   /transl_table=<integer>
##  
## Citations:
##   /citation=[number] e.g. /citation=[3]
##   /compare=[accession-number.sequence-version] e.g. /compare=AJ634337.1
##
gbFeatures <- function(gb_features, seqinfo) {
  # where do all the features start
  feature_start <- which(substr(gb_features, 6, 6) != " ")
  # where do all the features end
  feature_end <- c(feature_start[-1] - 1, length(gb_features))
  # indeces for all features
  feature_idx <- .mapply(seq.int, list(from=feature_start, to=feature_end), NULL)
  accession <- getAccession(seqinfo)
  
#   ftbl <- mapply(function(idx, n) {
#     gbFeature(gb_features[idx], seqinfo, accession, n)
#   }, idx=feature_idx, n=seq_along(feature_start),
#                    SIMPLIFY=FALSE, USE.NAMES=FALSE)
  
#   p <- MulticoreParam(workers=floor(detectCores()*0.75), verbose=TRUE)
#   ftbl <- bpmapply(function(idx, n) {
#     gbFeature(feature=gb_features[idx], seqinfo=seqinfo, accession=accession, id=n)
#   },
#   idx = feature_idx, n = seq_along(feature_start),
#   SIMPLIFY=FALSE, USE.NAMES=FALSE, PBPARAM=p)
  
  mc_cores <- floor(detectCores()*0.75)
  ftbl <- mcmapply(function(idx, n) {
    gbFeature(gb_features[idx], accession, n)
  },
  idx = feature_idx, n = seq_along(feature_start),
  SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=mc_cores, mc.silent=FALSE)
                   
  IRanges::new2('gbFeatureList', .Data=ftbl, .seqinfo=seqinfo, check=FALSE) 
}


join_seq <- function(seq, accession_no) {
  mc_cores <- floor(detectCores()*0.75)
  s <- unlist(mclapply(seq, function(x) {
    paste0(strsplit(substr(x, 11, 75), " ")[[1L]], collapse="")
  }, mc.cores=mc_cores))
  s <- c(paste0(">", accession_no), s)
  s
}

#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings readAAStringSet
#' @importFrom Biostrings BStringSet
gbSequence <- function(gb_sequence, accession_no, seq_type) {
  # read.BStringSet() does not support connections and
  # currently only accepts fasta format. So we write out gb_sequence as
  # a temporary fasta file and read it back in as an AAStringSet or
  # DNAStringSet (mRNA etc seems to be encoded with Ts rather then Us,
  # so we use DNAStringSets for RNA)
  if (is.null(gb_sequence)) {
    return(BStringSet())
  } else {
    tmp <- tempfile()
    on.exit(unlink(tmp))
    writeLines(text=join_seq(gb_sequence, accession_no), tmp)
    origin <- switch(seq_type,
                     AA=readAAStringSet(tmp, format="fasta"),
                     readDNAStringSet(tmp, format="fasta"))
    origin
  }
}

