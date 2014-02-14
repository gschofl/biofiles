#' @importFrom Biostrings writeXStringSet
#' @export
#' @rdname write.FeatureTable-methods
setMethod("write.FeatureTable", "gbRecord", 
          function(x, file, tablename = "", dbname = "",
                   sequence = FALSE, append = FALSE) {
            .write_feature_table(x = x, file = file, tablename = tablename,
                                 dbname = dbname, sequence = sequence,
                                 append = append)
          })

#' @export
#' @rdname write.FeatureTable-methods
setMethod("write.FeatureTable", "gbFeatureTable", 
          function(x, file, tablename = "", dbname = "",
                   sequence = FALSE, append = FALSE) {
            x <- as(x, "gbRecord")
            .write_feature_table(x = x, file = file, tablename = tablename,
                                 dbname = dbname, sequence = sequence,
                                 append = append)
          })

.write_feature_table <- function(x, file, tablename = "", dbname = "",
                                 sequence = FALSE, append = FALSE) {
  # write header
  header <- trim(sprintf(">Feature %s %s", getAccession(x), tablename))
  cat(paste0(header, sep = "\n"), file = file)
  
  # kick out source if present
  f_list <- .features(x)
  f_list <- f_list[key(f_list) != "source"]
  # get index of genes
  gene_idx <- index(f_list["gene"])
  # write features
  f_table <- unlist(lapply(f_list, .get_table_feature, gene_idx, dbname, f_list))
  cat(paste0(f_table, collapse = "\n"), file = file, append = TRUE)
  
  if (sequence) {
    seq <- getSequence(x)
    names(seq) <- getAccession(x)
    writeXStringSet(seq, filepath = replace_ext(file, "fna"), format = "fasta")
  }
  
  invisible(NULL)
}

 
.get_table_feature <- function(f, gene_idx, dbname, f_list) {
  l <- location(f)  
  loc <- get_loc(l@range[1, , drop = FALSE], l@fuzzy[1, , drop = FALSE], l@strand)
  loc_line <- sprintf("%s%s\t%s%s\t%s\n", loc[[1]], loc[[2]], loc[[3]], loc[[4]], key(f))
  loc_line2 <- if (nrow(l@range) > 1) {
    loc <- get_loc(l@range[-1, , drop = FALSE], l@fuzzy[-1, , drop = FALSE], l@strand)
    sprintf("%s%s\t%s%s\n", loc[[1]], loc[[2]], loc[[3]], loc[[4]])
  } else { 
    "" 
  }
  loc_line <- sprintf("%s%s", loc_line, paste0(loc_line2, collapse = ""))
  qua <- names(f@qualifiers)
  val <- unname(f@qualifiers)

  ## add protein_id if the feature is a CDS
  if (key(f) == "CDS" && "protein_id" %ni% qua) {
    locus_tag_idx <- match("locus_tag", qua)
    if (!is.na(locus_tag_idx)) {
      locus_tag <- val[locus_tag_idx]
    } else {
      previous_gene_idx <- max(gene_idx[gene_idx < f@.id])
      previous_gene <- f_list[previous_gene_idx]
      locus_tag <- unname(qualif(previous_gene, "locus_tag"))
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
  if (key(f) != "gene") {
    tag_idx <- qua %in% c("locus_tag", "gene")
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
  feature <- paste0(loc_line, paste0(qua_line, collapse = "\n"), sep = " ")
  feature 
}

get_loc <- function(l, p, strand) {
  if (unique(strand) == 1L) {
    start <- as.character(l[, 1])
    start_prefix <- ifelse(p[, 1], "<", "")
    end <- as.character(l[, 2])
    end_prefix <- ifelse(p[, 2], ">", "")
  } else {
    start <- as.character(l[, 2])
    start_prefix <- ifelse(p[, 2], ">", "")
    end <- as.character(l[, 1])
    end_prefix <- ifelse(p[, 1], "<", "")
  }
  list(start_prefix = start_prefix, start = start,
       end_prefix = end_prefix, end = end)
} 
