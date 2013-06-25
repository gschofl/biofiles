setOldClass("environment")

#' gbHeader-class
#' 
#' @rdname gbHeader
#' @export
#' @classHierarchy
#' @classMethods
setClass("gbHeader",
         representation(
           locus = "character", length = "numeric", moltype = "character",
           topology = "character", division = "character", update_date = "POSIXlt",
           create_date = "POSIXlt", definition = "character", accession = "character",
           version = "character", seqid = "character", dblink = "character",
           dbsource = "character", keywords = "character", source = "character",
           organism = "character", taxonomy = "character", references = "character",
           comment = "character"),
         prototype(
           locus=NA_character_, length=NA_real_, moltype=NA_character_,
           topology=NA_character_, division=NA_character_, update_date=as.POSIXlt(NA),
           create_date=as.POSIXlt(NA), definition=NA_character_, accession=NA_character_,
           version=NA_character_, seqid=NA_character_, dblink=NA_character_,
           dbsource=NA_character_, keywords=NA_character_, source=NA_character_,
           organism=NA_character_, taxonomy=NA_character_, references=NA_character_,
           comment=NA_character_))


#' seqinfo-class
#' 
#' @rdname seqinfo
#' @export
#' @classHierarchy
#' @classMethods
setClass("seqinfo", contains="environment",
         prototype={
           .xData <- new.env(parent=emptyenv())
           .xData$header <- new('gbHeader')
           .xData$sequence <- new('DNAStringSet')
           .xData
         })


setMethod('.header', 'seqinfo', function (x) get('header', x))


setMethod('.sequence', 'seqinfo', function (x) {
  tryCatch(get("sequence", x), error = function (e) {
    warning("No sequence associated with this object", call.=FALSE)
    new('BStringSet')
  })
})

setMethod("getLocus", "seqinfo", function (x) .header(x)@locus)


setMethod("getLength", "seqinfo", function (x) .header(x)@length)


setMethod("getMoltype", "seqinfo", function (x) .header(x)@moltype)


setMethod("getTopology", "seqinfo", function (x) .header(x)@topology)


setMethod("getDivision", "seqinfo", function (x) .header(x)@division)


setMethod("getDate", "seqinfo", function (x) {
  c(create_date = .header(x)@create_date,
    update_date = .header(x)@update_date)
})


setMethod("getDefinition", "seqinfo", function (x) .header(x)@definition)


setMethod("getAccession", "seqinfo", function (x) .header(x)@accession)


setMethod("getVersion", "seqinfo", function (x) .header(x)@version)


setMethod("getGeneID", "seqinfo", 
          function (x, db = 'gi') {
            seqid <- .header(x)@seqid
            if (is.na(seqid)) {
              seqid
            }
            else {
              db.idx <- which(strsplitN(seqid, "|", 1, fixed = TRUE) == db)
              strsplitN(seqid, "|", 2, fixed = TRUE)[db.idx]
            }
          })


setMethod("getDBLink", "seqinfo", function (x) .header(x)@dblink)


setMethod("getDBSource", "seqinfo", function (x) .header(x)@dbsource)


setMethod("getSource", "seqinfo", function (x) .header(x)@source)


setMethod("getOrganism", "seqinfo", function (x) .header(x)@organism)


setMethod("getTaxonomy", "seqinfo", function (x) .header(x)@taxonomy)


setMethod("getReference", "seqinfo", function (x) .header(x)@references)


setMethod("getKeywords", "seqinfo", function (x) .header(x)@keywords)


setMethod("getComment", "seqinfo", function (x) .header(x)@comment)


setMethod("getSequence", "seqinfo", function (x) .sequence(x))


setMethod("show", "seqinfo", function (object) {
  if (is.na(getAccession(object))) {
    acc <- len <- def <- ""
  }
  else {
    acc <- getAccession(object)
    len <- paste0(getLength(object), " ", getMoltype(object))
    def <- getDefinition(object)
    acc <- pad(acc, nchar(acc) + 2, "right")
    len <- pad(len, nchar(len) + 2, "right")
    def <- ellipsize(def, width=getOption("width") - 
                     nchar(acc) - nchar(len) - 3)
  }
  cat(sprintf("%s%s%s", acc, len, def))
})






