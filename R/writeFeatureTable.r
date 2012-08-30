
##' General function for writing out NCBI feature tables
##'
##' Feature tables are simple five-column tab-delimited tables specifying the
##' location and type of each feature. They can be used as input for tbl2asn
##' or Sequin to generate annotation.
##'
##' @param db A \code{gbRecord} object.
##' @param tablename (Optional) Optional table name to appear in the first line
##' of the feature table.
##' @param dbname Data base name associated with the CDS qualifier protein_id.
##' @param with_sequence Additionally autput fasta file 
##' @param outfile Output file.
##' 
##' @export
writeFeatureTable <- function(db, tablename="", dbname="",
                              with_sequence=TRUE, outfile="out.tbl") {
  
  if (file.exists(outfile)) {
    unlink(outfile)
  }
  
  fasta_outfile <- sub("tbl$", "fna", outfile)
  if (file.exists(fasta_outfile)) {
    unlink(fasta_info)
  }
  
  # write header
  header <- gsub("<\\s+|\\s+$", "",
                 sprintf(">Feature %s %s", db$accession, tablename))
  cat(paste(header, "\n"), file=outfile)
  
  # kick out source if present
  features <- dbFetch(db, "features")
  features <- features[getKey(features) != "source"]
  # get index of genes
  gene_idx <- getIndex(features[getKey(features) == "gene"])
  # write features
  f_table <- unlist(lapply(features, .getTableFeature, gene_idx, dbname))
  cat(paste(f_table, collapse="\n"), file=outfile, append=TRUE)
  
  if (with_sequence) {
    seq <- db$sequence
    names(seq) <- db$accession
    write.XStringSet(seq, filepath=fasta_outfile, format="fasta")
  }
  
  invisible(list(features=f_table))
}


.getTableFeature <- function (f, gene_idx, dbname) {
  
  getLoc <- function (l, p, strand) {
    if (strand == 1) {
      start <- as.character(l[,1])
      start_prefix <- ifelse(p[,1], "<", "")
      end <- as.character(l[,2])
      end_prefix <- ifelse(p[,2], ">", "")
    } else {
      start <- as.character(l[,2])
      start_prefix <- ifelse(p[,2], ">", "")
      end <- as.character(l[,1])
      end_prefix <- ifelse(p[,1], "<", "")
    }
    
    list(start_prefix=start_prefix, start=start,
         end_prefix=end_prefix, end=end)
  }  

  l <- f@location  
  loc <- getLoc(l@.Data[1,,drop=FALSE], l@partial[1,,drop=FALSE], l@strand)
  loc_line <- sprintf("%s%s\t%s%s\t%s\n", loc[[1]], loc[[2]], loc[[3]],
                      loc[[4]], f@key)
  loc_line2 <- if (nrow(l@.Data) > 1) {
    loc <- getLoc(l=l@.Data[-1,,drop=FALSE], p=l@partial[-1,,drop=FALSE],
                  strand=l@strand)
    sprintf("%s%s\t%s%s\n", loc[[1]], loc[[2]], loc[[3]],
            loc[[4]])
  } else { 
    "" 
  }
  loc_line <- sprintf("%s%s", loc_line, paste(loc_line2, collapse=""))
  qua <- names(f@qualifiers)
  val <- unname(f@qualifiers)

  ## add protein_id if the feature is a CDS
  if (f@key == "CDS" && "protein_id" %ni% qua) {
    locus_tag_idx <- match("locus_tag", qua)
    if (!is.na(locus_tag_idx)) {
      locus_tag <- val[locus_tag_idx]
    } else {
      previous_gene_idx <- max(gene_idx[gene_idx < f@.ID])
      previous_gene <- dbFetch(db, "features")[previous_gene_idx]
      locus_tag <- unname(select(previous_gene, cols="locus_tag"))
    }
    prot_id_line <- paste0("\t\t\tprotein_id\tgnl|", dbname, "|", locus_tag)
  } else {
    prot_id_line <- NULL
  }
  
  ## Qualifiers that can be used in CDS features
  #   product
  #   prot_desc
  #   function
  #   EC_number
  #   note
  #   experiment
  #   inference
  #   go_component
  #   go_process
  #   go_function
  #   db_xref
  #   pseudo
  #   exception
  #   transl_except
  #
  ## Therefor we remove
  #   locus_tag
  #   codon_start
  #   transl_table
  #   translation
  
  ## replace any "TRUE" feature values
  val[val %in% "TRUE"] <- ""
  
  ## remove locus_tag and gene from all but /gene qualifiers
  if (f@key != "gene") {
    tag_idx <- qua %in% c("locus_tag","gene")
    if (any(tag_idx)) {
      qua <- qua[!tag_idx]
      val <- val[!tag_idx] 
    }
  }
  
  ## remove codon_start, transl_table, translation
  qual_to_remove <- c("codon_start", "transl_table", "translation")
  qual_to_retain <- qua %ni% qual_to_remove
  qua <- qua[qual_to_retain]
  val <- val[qual_to_retain]
  
  qua_line <- c(sprintf("\t\t\t%s\t%s", qua, val), prot_id_line)
  feature <- paste0(loc_line, paste0(qua_line, collapse="\n"))
  feature 
}



