#' @importFrom Biostrings DNAStringSet
NULL
#' @importFrom IRanges new2
NULL
#' @importFrom parallel mclapply mcmapply detectCores
NULL

.parseGB <- function (gb_data, with_sequence = TRUE) {
  # get a vector with the positions of the main GenBank fields
  gbf <- grep("^[A-Z//]+", gb_data)
  gbf_names <- strsplitN(gb_data[gbf], split=" +", 1)
  names(gbf) <- gbf_names
  gb_contig <- gb_sequence <-  NULL
  
  # Check the presence of a number of the absolutely essential data fields
  essential <- c("DEFINITION", "ACCESSION", "FEATURES")
  if (any(is.na(charmatch(essential, gbf_names))))
    stop("Some fields seem to be missing from the GenBank file")
  
  # Split the GenBank file into HEADER, FEATURES, ORIGIN
  gb_header <- gb_data[seq.int(gbf["FEATURES"]) - 1]
  
  if (match("FEATURES", gbf_names) == length(gbf)) {
    gb_features <- gb_data[seq.int(gbf["FEATURES"] + 1, length(gb_data) - 1)]
  } else {
    gb_features <- gb_data[seq.int(gbf["FEATURES"] + 1, gbf[match("FEATURES", gbf_names) + 1] - 1)]
  }
  
  if (length(gb_features) < 2L) 
    stop("No features in the GenBank file")
  
  if (!is.na(origin <- gbf["ORIGIN"])) {
    seq_idx <- seq.int(origin + 1, gbf["//"] - 1)
    if (length(seq_idx) > 1L && seq_idx[2] < seq_idx[1]) {
      # happens if "//" is next after "ORIGIN", i.e. no  sequence is present
      gb_sequence <- NULL
    } else {
      gb_sequence <- gb_data[seq_idx]
    }
  } else if (!is.na(contig <- gbf["CONTIG"])) {
    gb_contig <- gbLocation(strsplitN(gb_data[contig], "CONTIG", 2))
  }

  # generate the seqenv environment that will hold the Seqinfo and Sequence
  # information. Set "seqinfo" and "sequence" to NULL
  seqenv <- new.env(parent=emptyenv())
  
  ## parse HEADER, FEATURES, and ORIGIN and construct 'gbData' object
  header <- .parseGbHeader(gb_header=gb_header, gb_fields=gbf)
  seqenv[["seqinfo"]] <- header[["seqinfo"]]
  
  sequence <- .parseGbSequence(gb_sequence=gb_sequence,
                               accession_no=seqnames(header[["seqinfo"]]),
                               seq_type=header[["moltype"]])
  seqenv[["sequence"]] <- sequence
  
  features <- .parseGbFeatures(gb_features=gb_features, seqenv=seqenv)
  list(header=header, features=features, contig=gb_contig, seqenv=seqenv)
}


.parseGbHeader <- function (gb_header, gb_fields) {
  #### LOCUS
  locus_line <- strsplit(gb_header[gb_fields[names(gb_fields) == "LOCUS"]], "\\s+")[[1]]
  # if split by whitespace the LOCUS field seems to be made up of up to 8
  # different elements, 6 of which are data: LOCUS (1), locus name (2)
  # sequence length (3), bp or aa (4), molecule type (5, not given if protein,
  # GenBank division (6), topology (7, optional?), and modification date (8)
  if (length(locus_line) == 8 || 
     (length(locus_line) == 7 && locus_line[4] == "aa"))
  {
    locus <- locus_line[2]
    length <- as.integer(locus_line[3])
    if (locus_line[4] == 'aa') {
      # these are GenPept files; they don't have a 'molecule type' but we set it 'AA'
      moltype <- 'AA'
      topology <- locus_line[5]
      division <- locus_line[6]
      update_date <- as.POSIXlt(locus_line[7], format="%d-%b-%Y")
    } else if (locus_line[4] == 'bp') {
      moltype <- locus_line[5]
      topology <- locus_line[6]
      division <- locus_line[7]
      update_date <- as.POSIXlt(locus_line[8], format="%d-%b-%Y")
    }
  } 
  else
  {
    # some GenBank files just don't seem to have topology or division
    locus <- locus_line[2]
    length <- as.integer(locus_line[3])
    if (locus_line[4] == 'aa') {
      moltype <- 'AA'
      topology <- locus_line[5]
      division <- NA_character_
      update_date <- as.POSIXlt(locus_line[6], format="%d-%b-%Y")
    } else if (locus_line[4] == 'bp') {
      moltype <- locus_line[5]
      topology <- NA_character_
      division <- locus_line[6]
      update_date <- as.POSIXlt(locus_line[7], format="%d-%b-%Y")
    }
  }
  
  #### DEFINITION
  def_idx <- which(names(gb_fields) == "DEFINITION")
  def_line <- gb_header[seq.int(gb_fields[def_idx], gb_fields[def_idx + 1] - 1)]
  definition <- paste(gsub(" +", " ", sub("DEFINITION  ", "", def_line)), collapse=" ")
  
  #### ACCESSION
  acc_line <- gb_header[gb_fields[names(gb_fields) == "ACCESSION"]]
  accession <- strsplit(acc_line, " +")[[1]][2]
  if (is.na(accession)) accession <- locus
  
  #### VERSION and GI
  ver_line <- gb_header[gb_fields[names(gb_fields) == "VERSION"]]
  version <- strsplit(ver_line, " +")[[1]][2]
  seqid <- paste0('gi|', strsplit(ver_line, "GI:")[[1]][2])
  
  #### DBLINK (not seen everywhere)
  if (length(db_line <- gb_header[gb_fields[names(gb_fields) == "DBLINK"]]) > 0L) {
    dblink <- strsplit(db_line, "Project: ")[[1]][2]
  } else {
    dblink <- NA_character_
  }
  
  #### DBSOURCE (only in GenPept files)
  # sometimes more than one line, but always followed by "KEYWORDS"?
  if (length(dbsrc_idx <- which(names(gb_fields) == "DBSOURCE")) > 0L) {
    dbs_lines <- 
      gb_header[seq.int(gb_fields[dbsrc_idx], gb_fields[dbsrc_idx + 1] - 1)]
    dbsource <- paste(gsub("^ +", "", sub("DBSOURCE", "", dbs_lines)), collapse="\n")
  } else {
    dbsource <- NA_character_
  }
  
  #### KEYWORDS
  key_line <- gb_header[gb_fields[names(gb_fields) == "KEYWORDS"]]
  keywords <- sub("KEYWORDS    ", "", key_line)
  
  #### SOURCE with ORGANISM and the complete lineage
  src_idx <- which(names(gb_fields) == 'SOURCE')
  source_lines <- 
    gb_header[seq.int(gb_fields[src_idx], gb_fields[src_idx + 1] - 1)]                  
  source <- sub("SOURCE      ", "", source_lines[1L])
  organism <- sub("  ORGANISM  ", "", source_lines[2L])
  taxonomy <- paste(gsub("^ +", "", source_lines[-c(1L,2L)]), collapse=" ")
  
  #### REFERENCES (not seen everywhere)
  if (length(ref_idx <- which(names(gb_fields) == "REFERENCE")) > 0L) {
    
    ref_lines <- gb_header[seq.int(gb_fields[ref_idx[1]],
                                   gb_fields[ref_idx[length(ref_idx)] + 1] - 1)]
    references <- .parseGbReferences(ref_lines)
  } else {
    references <- "Not available"
  }
  
  #### COMMENT (not always there)
  if (length(gb_fields[names(gb_fields) == "COMMENT"]) > 0L) {
    com_lines <- 
      gb_header[seq.int(gb_fields[names(gb_fields) == "COMMENT"],
                        length(gb_header))]
    comment <- paste(gsub("^ +", "", sub("COMMENT", "", com_lines)), collapse="\n")
  } else {
    comment <- NA_character_
  }

  #### Seqinfo
  seqinfo <- Seqinfo(seqnames=accession, seqlengths=length,
                     isCircular=topology == 'circular',
                     genome=definition)

  # References are assigned to the 'gb_reference' class
  list(seqinfo=seqinfo, locus=locus, moltype=moltype, topology=topology,
       division=division, update_date=update_date, create_date=as.POSIXlt(NA),
       version=version, seqid=seqid, dblink=dblink, dbsource=dbsource,
       keywords=keywords, source=source, organism=organism, taxonomy=taxonomy,
       references=references, comment=comment)
} 

.parseGbReferences <- function (ref_lines) {
  return("Not implemented yet")
  ref_idx <- grep("^REFERENCE", ref_lines)
  idx <- Map(seq, ref_idx, c(ref_idx[-1] - 1, length(ref_lines)))
  ref <- ref_lines[idx[[1]]]
  auth_line <- grep("^  AUTHORS", ref, value=TRUE)
  .parseAuthors <- function (auth_line) {}
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
.parseGbFeatures <- function (gb_features, seqenv) {
  # where do all the features start
  feature_start <- which(substr(gb_features, 6, 6) != " ")
  # where do all the features end
  feature_end <- c(feature_start[-1] - 1, length(gb_features))
  # indeces for all features
  feature_idx <- mapply(seq.int, feature_start, feature_end,
                        SIMPLIFY=FALSE, USE.NAMES=FALSE)
  ftbl <- mcmapply(function (idx, n) {
    gbFeature(gb_features[idx], seqenv, n)
  }, idx=feature_idx, n=seq_along(feature_start),
                  SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=detectCores())
  
  IRanges::new2('gbFeatureList', .Data=ftbl, .seqinfo=seqenv, check=FALSE) 
}

.joinSeq <- function (seq, accession_no) {
  mc_cores <- detectCores()
  s <- unlist(mclapply(seq, function(x) {
    paste0(strsplit(substr(x, 11, 75), " ")[[1L]], collapse="")
  }, mc.cores=mc_cores))
  s <- c(paste0(">", accession_no), s)
  s
}

#' @autoImports
.parseGbSequence <- function (gb_sequence, accession_no, seq_type) {
  # read.BStringSet() does not support connections and
  # currently only accepts fasta format. So we write out gb_sequence as
  # a temporary fasta file and read it back in as an AAStringSet or
  # DNAStringSet (mRNA etc seems to be encoded with Ts rather then Us,
  # so we use DNAStringSets for RNA)
  if (is.null(gb_sequence)) {
    return( BStringSet() )
    } else {
    tmp <- tempfile()
    on.exit(unlink(tmp))
    writeLines(text=.joinSeq( gb_sequence, accession_no ), tmp)
    origin <- switch(seq_type,
                     AA=readAAStringSet(tmp, format="fasta"),
                     readDNAStringSet(tmp, format="fasta"))
    origin
  }
}

