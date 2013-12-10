#' @importFrom Biostrings DNAStringSet
#' @importFrom IRanges new2
#' @importFrom parallel mclapply mcmapply detectCores
#' @importFrom foreach foreach registerDoSEQ "%dopar%"
#' @importFrom iterators iter
NULL

.mandatory <- c("LOCUS", "DEFINITION", "ACCESSION", "VERSION", "FEATURES", "//")

parse_gb_record <- function(gb_record) {
  ## check that the gb_record is not empty
  if (length(gb_record) == 0L) {
    stop("This \"gbRecord\" is empty.")
  }
  ## check if gb_record contains multiple entries
  end_of_record <- grep('^//$', gb_record)
  n_records <- length(end_of_record)
  if (n_records > 1L) {
    start_of_record <- c(1, end_of_record[-n_records] + 1)
    irec <- iter(ixsplit(gb_record, start_of_record))
  } else {
    registerDoSEQ()
    irec <- iter(list(gb_record))
  }
  gbk_list <- foreach(rec = irec, .inorder=FALSE) %dopar% {
    # get a vector with the positions of the main GenBank fields
    rec_idx <- grep("^[A-Z//]+", rec)
    rec_kwd <- strsplitN(rec[rec_idx], " +", 1L)
    gb_contig <- gb_sequence <-  NULL
    # Check the presence of mandatory fields
    if (any(is.na(charmatch(.mandatory, rec_kwd)))) {
      stop("mandatory fields are missing from the GenBank file")
    }
    ## get positions of features, origin, contig and end_of_record
    ftb_idx <- rec_idx[rec_kwd == "FEATURES"]
    ori_idx <- rec_idx[rec_kwd == "ORIGIN"]
    ctg_idx <- rec_idx[rec_kwd == "CONTIG"]
    end_idx <- rec_idx[rec_kwd == "//"]
    ftb_end_idx <- rec_idx[which(rec_kwd == "FEATURES") + 1] - 1
    
    ## HEADER
    seqenv <- seqinfo(gbHeader(gb_header=rec[seq.int(ftb_idx - 1)]), NULL)
    ## SEQUENCE
    if (length(ori_idx) > 0L) {
      # if "//" is right after "ORIGIN" there is no sequence
      # and gb_sequence stays set to NULL
      if (end_idx - ori_idx > 1L) {
        gb_sequence <- rec[seq.int(ori_idx + 1, end_idx - 1)]
      }
      ## CONTIG
    } else if (length(ctg_idx) > 0L) {
      contig_line <- strsplitN(collapse(rec[seq.int(ctg_idx, end_idx-1)], ''),
                               'CONTIG', 2L, fixed=TRUE)
      gb_contig <- gbLocation(contig_line)
    }
    seqenv$sequence <- gbSequence(gb_sequence, getAccession(seqenv), getMoltype(seqenv))
    ## FEATURES
    gb_features <- rec[seq.int(ftb_idx + 1, ftb_end_idx)]
    gb_features <- gbFeatures(gb_features, seqinfo=seqenv)
    new_gbRecord(seqinfo=seqenv, features=gb_features, contig=gb_contig) 
  }
   
  if (length(gbk_list) == 1L) {
    gbk_list[[1L]]
  } else {
    gbRecordList(gbk_list)
  }
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
  feature_start <- which(substr(gb_features, 6, 6) != " ")
  fl <- ixsplit(gb_features, feature_start)
  mc_cores <- floor(detectCores()*0.75)
  ftbl <- mcmapply(gbFeature, feature=fl, id=seq_along(feature_start),
                   MoreArgs=list(accession=getAccession(seqinfo)),
                   SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=mc_cores)                   
  IRanges::new2('gbFeatureList', .Data=ftbl, .seqinfo=seqinfo, check=FALSE) 
}


join_seq <- function(seq, accession_no) {
  mc_cores <- floor(detectCores()*0.75)
  s <- unlist(mclapply(seq, function(x) {
    collapse(strsplit(substr(x, 11, 75), " ")[[1L]], '')
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

