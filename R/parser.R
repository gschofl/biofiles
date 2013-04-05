.parseGB <- function (gb_data, db_path, with_sequence = TRUE, force = FALSE) {
  # get a vector with the positions of the main GenBank fields
  gbf <- grep("^[A-Z//]+", gb_data)
  gbf_names <- strsplitN(gb_data[gbf], split=" +", 1)
  names(gbf) <- gbf_names
  
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
  
  seq_idx <- seq.int(gbf["ORIGIN"] + 1, gbf["//"] - 1)
  if (length(seq_idx) > 1L && seq_idx[2] < seq_idx[1]) {
    # happens if "//" is next after "ORIGIN", i.e. no  sequence is present
    gb_sequence <- NULL
  } else {
    gb_sequence <- gb_data[seq_idx]
  }
  
  # prepare setting up a filehash db
  if (file.exists(db_path)) {
    if (!force && readline("Target directory exists. Do you want to overwrite (y/n): ") != "y")
      stop(gettextf("Aborted creation of db directory '%s'", db_path))
    else
      unlink(db_path, recursive=TRUE)
  }
  
  # set up filehash db
  db <- init_db(db_path, create = TRUE)
  
  ## parse HEADER, FEATURES, and ORIGIN and construct 'gbData' object
  header <- .parseGbHeader(gb_header, gb_fields=gbf)
  seqinfo <- IRanges::new2("gbInfo", db=db,
                           seqnames=header[["accession"]],
                           seqlengths=header[["length"]],
                           is_circular=header[["topology"]] == "circular", 
                           genome=header[["definition"]],
                           check=FALSE)
  features <- .parseGbFeatures(seqinfo=seqinfo, gb_features=gb_features)
  
  sequence <- .parseGbSequence(gb_sequence=gb_sequence,
                               accession_no=header[["accession"]],
                               seq_type=header[["type"]])
  
  list(db=db, header=header, features=features, sequence=sequence)
}


.parseGbHeader <- function (gb_header, gb_fields) {
  #### LOCUS
  locus_line <- strsplit(gb_header[gb_fields[names(gb_fields) == "LOCUS"]], "\\s+")[[1]]
  # if split by whitespace the LOCUS field seems to be made up of up to 8
  # different elements, 6 of which are data: LOCUS (1), locus name (2)
  # sequence length (3), bp or aa (4), molecule type (5, not given if protein,
  # GenBank division (6), topology (7, optional?), and modification date (8)
  if (length(locus_line) == 8 || 
     (length(locus_line) == 7 && locus_line[4] == "aa")) {
    locus <- locus_line[2]
    length <- setNames(as.integer(locus_line[3]), locus_line[4])
    
    if (names(length) == 'aa') {
      # these are GenPept files; they don't have a 'molecule type' but we set it 'AA'
      type <- 'AA'
      topology <- locus_line[5]
      division <- locus_line[6]
      date <- locus_line[7]
    } else if (names(length) == 'bp') {
      type <- locus_line[5]
      topology <- locus_line[6]
      division <- locus_line[7]
      date <- as.POSIXlt(locus_line[8], format="%d-%b-%Y")
    }
  } else {
    # some GenBank files just don't seem to have topology or division
    locus <- locus_line[2]
    length <- setNames(as.integer(locus_line[3]), locus_line[4])
    if (names(length) == 'aa') {
      type <- 'AA'
      topology <- locus_line[5]
      division <- NULL
      date <- as.POSIXlt(locus_line[6], format="%d-%b-%Y")
    } else if (names(length) == 'bp') {
      type <- locus_line[5]
      topology <- NULL
      division <- locus_line[6]
      date <- as.POSIXlt(locus_line[7], format="%d-%b-%Y")
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
  GI <- strsplit(ver_line, "GI:")[[1]][2]
  
  #### DBLINK (not seen everywhere)
  if (length(db_line <- gb_header[gb_fields[names(gb_fields) == "DBLINK"]]) > 0L) {
    dblink <- strsplit(db_line, "Project: ")[[1]][2]
  } else {
    dblink <- NULL
  }
  
  #### DBSOURCE (only in GenPept files)
  # sometimes more than one line, but always followed by "KEYWORDS"?
  if (length(dbsrc_idx <- which(names(gb_fields) == "DBSOURCE")) > 0L) {
    dbs_lines <- 
      gb_header[seq.int(gb_fields[dbsrc_idx], gb_fields[dbsrc_idx + 1] - 1)]
    dbsource <- paste(gsub("^ +", "", sub("DBSOURCE", "", dbs_lines)), collapse="\n")
  } else {
    dbsource <- NULL
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
  lineage <- paste(gsub("^ +", "", source_lines[-c(1L,2L)]), collapse=" ")
  
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
    comment <- NULL
  }

  # References are assigned to the 'gb_reference' class
  list(locus=locus, length=length, type=type, topology=topology,
       division=division, date=date, definition=definition,
       accession=accession, version=version, GI=GI, dblink=dblink,
       dbsource=dbsource, keywords=keywords, source=source,
       organism=organism, lineage=lineage, references=references,
       comment=comment)
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
#' @autoImports
.parseGbFeatures <- function (seqinfo, gb_features) {
  # where do all the features start
  feature_start <- which(substr(gb_features, 6, 6) != " ")
  # where do all the features end
  feature_end <- c(feature_start[-1] - 1, length(gb_features))
  # indeces for all features
  feature_idx <- mapply(seq.int, feature_start, feature_end,
                        SIMPLIFY=FALSE, USE.NAMES=FALSE)
  
  message("Parsing features into ", dQuote(seqinfo@db@name))
  
#   idx <- feature_idx[[1]]
#   n <- 1
  
  ftr <- mcmapply(function (idx, n) {
    parse_feature_table(id=n, lines=gb_features[idx], seqinfo=seqinfo)
  }, idx=feature_idx, n=seq_along(feature_start),
     SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=detectCores())
  
  IRanges::new2('gbFeatureList', .Data=ftr, .Info=seqinfo, check=FALSE) 
}


#' @autoImports
.parseGbSequence <- function (gb_sequence, accession_no, seq_type) {
  # read.BStringSet() does not support connections and
  # currently only accepts fasta format. So we write out gb_sequence as
  # a temporary fasta file and read it back in as an AAStringSet or
  # DNAStringSet (mRNA etc seems to be encoded with Ts rather then Us,
  # so we use DNAStringSets for RNA)
  if (is.null(gb_sequence)) {
    return(NULL)
  } else {
    tmp <- tempfile()
    on.exit(unlink(tmp))
    writeLines(text=.joinSeq(gb_sequence, accession_no), tmp)
    origin <- switch(seq_type,
                     AA=readAAStringSet(tmp, format="fasta"),
                     readDNAStringSet(tmp, format="fasta"))
    origin
  }
}


#' @autoImports
.joinSeq <- function (seq, accession_no) {
  mc_cores <- detectCores()
  s <- unlist(mclapply(seq, function(x) {
    paste0(strsplit(substr(x, 11, 75), " ")[[1L]], collapse="")
  }, mc.cores=mc_cores))
  s <- c(paste0(">", accession_no), s)
  s
}

