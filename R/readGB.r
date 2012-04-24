##' General function for importing GenBank flat files
##'
##' For a description of the GenBank format see
##' \url{http://www.ncbi.nlm.nih.gov/collab/FT/}
##'
##' @usage readGB(gb, with_sequence=TRUE, force=FALSE)
##'
##' @param gb Path to a GenBank flat file or an 
##' \code{\link[rentrez]{efetch-class}} object containing GenBank record(s).
##' @param with_sequence If \code{TRUE}, sequence information
##' will be included.
##' @param force If \code{TRUE} existing database directories are
##' overwritten without prompting.
##' 
##' @return A (list of) \code{\link{gbRecord-class}} object(s).
##' 
##' @export
readGB <- function (gb,
                    with_sequence=TRUE,
                    force=FALSE)
{
  # we store efetch data in temporary files
  if (is(gb, "efetch")) {
    ## we can parse rettype = gbwithparts, gb, gp and retmode =  text
    if (!grepl("^gb|^gp", gb@type) || gb@mode != "text")
      stop("Must use efetch with rettype='gbwithparts','gb', or 'gp' and retmode='text'")
    
    split_gb <- strsplit(gb@data, split="\n\n")[[1L]]
    n <- length(split_gb)
    db_path <- replicate(n, tempfile(fileext=".db"))
    
    for (i in seq_len(n)) {
      gb_data <- strsplit(split_gb[i], "\n")[[1L]]
      cat(gettextf("Importing into %s\n", dQuote(basename(db_path[i]))))
      .parseGB(gb_data, db_path[i], with_sequence=with_sequence, force=force)
    }
  }
  else if (!isS4(gb) && file.exists(gb)) {
    con <- file(gb, open="rt")
    on.exit(close(con))
    db_path <- paste0(gb, ".db")
    
    cat(gettextf("Importing into %s\n", dQuote(basename(db_path))))
    .parseGB(readLines(con), db_path,
             with_sequence=with_sequence,
             force=force)
  }
  else {
    stop("'gb' must be a valid GenBank flat file or an 'efetch' object containing GenBank records")
  }
  return(invisible(db_path))
}

# Helper functions ----------------------------------------------------

.parseGB <- function (gb_data,
                      db_path,
                      with_sequence=TRUE,
                      force=FALSE)
{
  # get a vector with the positions of the main fields
  gb_fields <- grep("^[A-Z//]+", gb_data)
  names(gb_fields) <- regmatches(gb_data[gb_fields], regexpr("^.[^ ]+", gb_data[gb_fields]))
  
  # Check the presence of a number of the absolutely essential data fields
  essential_fields <- "DEFINITION|ACCESSION|FEATURES"
  if (sum(grepl(essential_fields, names(gb_fields))) < 3)
    stop("Some fields seem to be missing from the GenBank file")
  
  # Split the GenBank file into HEADER, FEATURES, ORIGIN
  gb_header <- gb_data[seq(gb_fields["FEATURES"])-1]
  
  if (which(names(gb_fields) == "FEATURES") == length(gb_fields)) {
    gb_features <- gb_data[seq(gb_fields["FEATURES"]+1, length(gb_data)-1)]
  } else {
    gb_features <- gb_data[seq(gb_fields["FEATURES"]+1, gb_fields[which(names(gb_fields) == "FEATURES")+1]-1)]
  }
  
  if (length(gb_features) < 2) 
    stop("No features in the GenBank file")
  
  seq_idx <- seq(gb_fields["ORIGIN"]+1, gb_fields["//"]-1)
  if (seq_idx[2] < seq_idx[1]) {
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
  
  ## parse HEADER, FEATURES, and ORIGIN and construct 'gbData' object
  header <- .parseGbHeader(gb_header, gb_fields)
  features <- .parseGbFeatures(db_dir=db_path, accession=header$accession, 
                               definition=header$definition, 
                               gb_features=gb_features)
  sequence <- .parseGbSequence(gb_sequence, header$accession, header$type)
  gbRecord(db_dir=db_path, header=header, features=features, sequence=sequence)
}

.parseGbHeader <- function (gb_header, gb_fields)
{
  #### LOCUS
  locus_line <- strsplit(gb_header[gb_fields[names(gb_fields) == "LOCUS"]], split=" +")[[1]]
  # if split by whitespace the LOCUS field seems to be made up of up to 8
  # different elements, 6 of which are data: LOCUS (1), locus name (2)
  # sequence length (3), bp or aa (4), molecule type (5, not given if protein,
  # GenBank division (6), topology (7, optional?), and modification date (8)
  if (length(locus_line) == 8 || 
    (length(locus_line) == 7 && locus_line[4] == "aa")) {
    locus <- locus_line[2]
    length <- as.numeric(locus_line[3])
    names(length) <- locus_line[4]
    if (identical(locus_line[4], "aa")) {
      # these are GenPept files; they don't have a 'moecule type' but we set it 'AA'
      type <- "AA"
      topology <- locus_line[5]
      division <- locus_line[6]
      date <- locus_line[7]
    } else {
      type <- locus_line[5]
      topology <- locus_line[6]
      division <- locus_line[7]
      date <- as.POSIXlt(locus_line[8], format="%d-%b-%Y")
    }
  }
  else {
    # some GenBank files just don't seem to have topology
    warning("Sequence topology might be missing from the locus definition.\nBetter check the result.")
    locus <- locus_line[2]
    length <- as.numeric(locus_line[3])
    name(length) <- locus_line[4]
    type <- locus_line[5]
    topology <- NULL
    division <- locus_line[6]
    date <- locus_line[7]
  }
  #### DEFINITION
  def_line <- gb_header[seq.int(gb_fields[names(gb_fields) == "DEFINITION"],
                                gb_fields[names(gb_fields) == "ACCESSION"] - 1L)]
  definition <- paste(gsub(" +", " ", sub("DEFINITION  ", "", def_line)), collapse=" ")
  #### ACCESSION
  acc_line <- gb_header[gb_fields[names(gb_fields) == "ACCESSION"]]
  accession <- strsplit(acc_line, split=" +")[[1]][2]
  #### VERSION and GI
  ver_line <- gb_header[gb_fields[names(gb_fields) == "VERSION"]]
  version <- strsplit(ver_line, split=" +")[[1]][2]
  GI <- strsplit(ver_line, split="GI:")[[1]][2]
  #### DBLINK (not seen everywhere)
  if (length(db_line <- gb_header[gb_fields[names(gb_fields) == "DBLINK"]]) > 0L)
    dblink <- strsplit(db_line, "Project: ")[[1]][2]
  else
    dblink <- NULL
  #### DBSOURCE (only in GenPept files)
  # sometimes more than one line, but always followed by "KEYWORDS"?
  if (length(gb_fields[names(gb_fields) == "DBSOURCE"]) > 0L) {
    dbs_lines <- 
      gb_header[seq.int(gb_fields[names(gb_fields) == "DBSOURCE"],
                        gb_fields[names(gb_fields) == "KEYWORDS"] - 1L)]
    dbsource <- paste(gsub("^ +", "", sub("DBSOURCE", "", dbs_lines)), collapse="\n")
  }
  else
    dbsource <- NULL
  #### KEYWORDS
  key_line <- gb_header[gb_fields[names(gb_fields) == "KEYWORDS"]]
  keywords <- sub("KEYWORDS    ", "", key_line)
  #### SOURCE with ORGANISM and the complete lineage
  source_lines <- 
    gb_header[seq.int(gb_fields[names(gb_fields) == "SOURCE"], 
                      gb_fields[names(gb_fields) == "REFERENCE"][1L] - 1L)]
  source <- sub("SOURCE      ", "", source_lines[1L])
  organism <- sub("  ORGANISM  ", "", source_lines[2L])
  lineage <- paste(gsub("^ +", "", source_lines[-c(1L,2L)]), collapse=" ")
  #### COMMENT (not always there)
  if (length(gb_fields[names(gb_fields) == "COMMENT"]) > 0L) {
    com_lines <- 
      gb_header[seq.int(gb_fields[names(gb_fields) == "COMMENT"],
                        length(gb_header))]
    comment <- paste(gsub("^ +", "", sub("COMMENT", "", com_lines)), collapse="\n")
  }
  else
    comment <- NULL
  #### REFERENCE parsing references isn't implemented yet
  # References are assigned to the 'gb_reference' class
  return(list(locus=locus, length=length, type=type, topology=topology,
              division=division, date=date, definition=definition,
              accession=accession, version=version, GI=GI, dblink=dblink,
              dbsource=dbsource, keywords=keywords, source=source,
              organism=organism, lineage=lineage, references=list("Not implemented yet"),
              comment=comment))
} 

.parseGbFeatures <- function (db_dir, accession, definition, gb_features)
{
  # where do all the features start
  feature_start <- grep("^\\s{5}[[:alpha:]]+", gb_features)
  # where do all the features end
  feature_end <- c(feature_start[-1] - 1, length(gb_features))
  # indeces for all features
  feature_idx <- Map(seq.int, feature_start, feature_end)
  
#   f_list <- mapply( function (idx, n) {
  f_list <- mcmapply( function (idx, n) {
    .parseFeatureField(db_dir=db_dir, accession=accession,
                       definition=definition, id=n,
                       lines=gb_features[idx])
  }, feature_idx, seq_along(feature_start),
                      SIMPLIFY=FALSE,
                      USE.NAMES=FALSE,
                      mc.cores=detectCores())
    
  ans <- gbFeatureList(db_dir=db_dir, accession=accession,
                       definition=definition, features=f_list)
  ans
}

.parseFeatureField <- function (db_dir, accession, definition, id, lines,
                                key_pat="(?<=^\\s{5})\\S+") {
  key <- regmatches(lines[1], regexpr(key_pat, lines[1], perl=TRUE))
  ## concatenate feature locations if they span multiple lines
  ## loc[[2]] contains the number of lines concatenated (mostly 1 anyways)
  ## loc[[1]] contains the feature's actual base span
  loc <- .joinLocation(lines)
  qual <- .joinQualifiers(lines[-seq.int(loc[[2L]])])
  ans <- gbFeature(db_dir=db_dir, accession=accession, definition=definition,
                   id=id, key=key, location=loc[[1L]], qualifiers=qual)
  ans
}

.parseGbSequence <- function(gb_sequence, accession_no, seq_type) {
  # read.BStringSet() does not support connections and
  # currently only accepts fasta format. So we write out gb_sequence as
  # a temporary fasta file and read it back in as a DNAStringSet,
  # RNAStringSet, AAStringSet depending on seq_type (DNA, AA, and
  # everyting else is RNA)
  if (is.null(gb_sequence)) {
    return(NULL)
  } else {
    tmp <- tempfile()
    on.exit(unlink(tmp))
    writeLines(text=.joinSeq(gb_sequence, accession_no), con=tmp)
    origin <- switch(seq_type,
              DNA=read.DNAStringSet(tmp, format="fasta"),
              AA=read.AAStringSet(tmp, format="fasta"),
              read.RNAStringSet(tmp, format="fasta"))
    origin
  }
}

.joinLocation <- function (lines, loc_pat="\\b\\S+$") {
  if (grepl("\\b\\S+[^,]$", lines[1]))
    return(list(regmatches(lines[1],
                           regexpr(loc_pat, lines[1], perl=TRUE)), 1))
  else
    return(joinLines(lines, extract_pat=loc_pat,
                     break_pat="\\b\\S+[^,]$", sep=FALSE))
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
##   characters (ASCII values 32-126 decimal); interal quotation markes are
##   ""escaped"" by placing a second quotation mark immediately before.
##
## Controlled vocabulary/Enumerators:
##   /anticodon=(pos:<base_range>,aa:<amino_acid>)
##   /codon_start=1, 2, or 3
##   /direction=left, right, or both
##   /estimated_length=unknown or <integer>
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
.joinQualifiers <- function (lines) {
  i <- 0
  Q <- c()
  Q <- eval(function (lines, 
                      qual_pat="(?<=^\\s{21}/).[^=]+\\b", # everything between (/) and (=)
                      text_pat="(?<=\\=\")(.+)(?=\"$)", # everything between ("") after (=) before EOL
                      ctrl_voc="^ {21}/(anticodon|codon_start|direction|estimated_length|
                      mod_base|number|rpt_type|rpt_unit_range|tag_peptide|
                      transl_except|transl_table|citation|compare)\\=.+") {

    if (!grepl("^ {21}/", lines[1])) return(Q)
    if (grepl(text_pat, lines[1], perl=TRUE)) {
      # standard free text enclosed by double quotation marks
      v <- regmatches(lines[1], regexpr(text_pat, lines[1], perl=TRUE)) 
      names(v) <- regmatches(lines[1], regexpr(qual_pat, lines[1], perl=TRUE))
      i <<- i + 1
      Q <- c(v, Recall(lines[-1])) }
    else if (grepl(ctrl_voc, lines[1], perl=TRUE)) {
      # any one of the cases of controled vocabulary
      v <- regmatches(lines[1], regexpr("(?<=\\=)(.+$)", lines[1], perl=TRUE))
      names(v) <- sub(ctrl_voc, "\\1", lines[1])
      i <<- i + 1
      Q <- c(v, Recall(lines[-1]))}
    else if (!grepl("^ {21}/.+=", lines[1])) {
      # qualifiers without value
      v <- TRUE
      names(v) <- regmatches(lines[1], regexpr("(?<= {21}/).+", lines[1], perl=TRUE))
      i <<- i + 1
      Q <- c(v, Recall(lines[-1]))}
    else {
      # as a last resort join lines; if we are at the /translation qualifier
      # remove the whitespace between lines
      j_line  <- joinLines(lines, extract_pat="(?<=^\\s{21}).+", break_pat="\"$",
                           sep=!grepl("/translation", lines[1]))
      v <- regmatches(j_line[[1]], regexpr(text_pat, j_line[[1]], perl=TRUE))
      names(v) <- regmatches(j_line[[1]], regexpr("(?<=^/).[^=]+\\b", j_line[[1]], perl=TRUE))         
      i <<- i + j_line[[2]]
      Q <- c(v, Recall(lines[-seq.int(j_line[[2]])]))
    } 
  }) (lines)
}

.joinSeq <- function (seq, accession_no) {
  mc_cores <- detectCores()
  s <- unlist(mclapply(strsplit(substring(text=seq, first=11, last=75), " "), 
                     function (x) { 
                       paste(x, collapse="") 
                     }, mc.cores=mc_cores))
  # paste0 pastes more efficiently than paste(..., sep="")
  s <- c(paste0(">", accession_no), s)
  s
}

####
# --R-- vim:ft=r:sw=2:sts=2:ts=4:tw=76:
#       vim:fdm=marker:fmr={{{,}}}:fdl=0
