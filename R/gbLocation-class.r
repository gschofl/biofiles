
# gbLocation-class ----------------------------------------------------

##' @importClassesFrom intervals Intervals_full
##' @importClassesFrom intervals Intervals_virtual
##' @import stringr
NULL

##' gbLocation class
##' 
##' gbLocation is a container for GenBank Feature Locations. It extends
##' \code{\link[Intervals]{Intervals_full}} and has xxx additional Slots:
##' \describe{
##'   \item{strand}{An integer code for minus (-1) or plus (1) strand}
##'   \item{compound}{A character code specifying how multiple segments
##'   are joined. One of 'join' or 'order'}
##'   \item{partial}{A logical matrix specifying whether residues are
##'   missing from the 5' and 3' ends respectively.}
##'   \item{accession}{}
##'   \item{remote}{}
##' }
##' 
##' For more information see the 
##' \href{ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt}{GenBank Release Note}
##'
##' @exportClass gbLocation
##' @name gbLocation-class
##' @rdname gbLocation-class
##' @aliases show,gbLocation-method
.gbLocation <- 
  #### gbLocation ####
  setClass("gbLocation",
           representation( strand = "integer",
                           compound =  "character",
                           partial = "matrix",
                           accession = "character",
                           remote = "logical"),
           prototype( type = "Z",
                      strand = NA_integer_,
                      compound = NA_character_,
                      partial = matrix( FALSE, 0, 2 ),
                      accession = NA_character_,
                      remote = FALSE ),
           contains = "Intervals_full",
           validity = function( object ) {
             if ( !all(object@strand %in% c(1L, -1L, NA_integer_)) )
               return( "The 'strand' slot should contain -1, 1, or NA" )
             if ( !all(object@compound %in% c("join","order",NA_character_)) )
               return( "The 'compound' slot should contain 'join', 'order', or NA" )
             
             return( TRUE )
           })

##' @keywords internal
setMethod("initialize",
          #### initialize-method ####
          signature(.Object = "gbLocation"),
          function (.Object, .Data, strand, compound, partial, remote, ...) 
          {
            if ( missing(.Data)) {
              callNextMethod(.Object, ...)
            } else {
              if ( !is.matrix( .Data ) )
                .Data <- matrix( .Data, ncol = 2 )
              
              if ( missing(strand) )
                strand <- NA_integer_
              if ( all(strand %in% c("+","-")) ) {
                strand <- if ( strand == "+") 1L else -1L
              } else if ( all(strand %in% c(1,-1)) ) {
                strand <- as.integer(strand)
              }
              
              if ( missing(partial) )
                partial <- matrix( FALSE, nrow(.Data), 2 )
              if ( is.vector(partial) ) {
                if ( length(partial) > 2 )
                  stop( "The 'partial' argument should be a matrix, or a vector of length 1 or 2." )
                partial <- matrix(
                  if ( nrow(.Data) == 0 ) logical() else partial,
                  nrow=nrow(.Data),
                  ncol=2, byrow=TRUE )
              }
              
              if ( missing(compound) && nrow(.Data) == 1L )
                compound <- NA_character_
              
              if ( missing(remote) )
                remote <- rep(FALSE, nrow(.Data))
              else if ( length(remote) != nrow(.Data) )
                remote <- c(rep(remote, nrow(.Data)%/%length(remote)),
                            remote[seq_len(nrow(.Data)%%length(remote))])  
              
              callNextMethod(.Object, .Data=.Data, strand=strand,
                             compound=compound, partial=partial,
                             remote=remote, ...)
            }
          })


setAs("gbLocation", "character",
      #### as.character-method ####
      function(from) {
        if (nrow(from) == 0)
          return(character())
        else {
          clo <- closed(from)
          par <- partial(from)
          str <- from@strand
          cmp <- from@compound
          acc <- from@accession
          rem <- from@remote
          
          span <- ifelse(clo[,1],
                         "..", 
                         ifelse(from[,2] == from[,1] + 1,
                                "^",
                                ".")
          )
          
          pos <- ifelse(from[,1] == from[,2],
                        from[,1], 
                        paste0(
                          ifelse( par[,1], "<", "" ),
                          from[,1],
                          span,
                          ifelse( par[,2], ">", "" ),
                          from[,2]
                        )
          )
          
          pos <- ifelse( rem,
                         paste0(acc, ":", pos),
                         pos)
          
          res <- 
            if (length(str) == 1) {
              paste0(
                ifelse( identical(str, -1L), "complement(", ""),
                ifelse( !is.na(cmp), paste0(cmp, "("), ""),
                paste0(pos, collapse=","),
                ifelse( !is.na(cmp), ")", ""),
                ifelse( identical(str, -1L), ")", "")
              )
            } else if (length(str) == nrow(from)) {
              paste0(
                ifelse( !is.na(cmp), paste0(cmp, "("), ""),
                paste0(
                  ifelse( str == -1L,
                          paste0("complement(", pos, ")"),
                          pos),
                  collapse = ","),
                ifelse( !is.na(cmp), ")", "")
              )  
            }
          
          return(res)
        }
      })

##' @keywords internal
setGeneric("partial",
           #### partial-generic ####
           function(x) {
             standardGeneric("partial")
           })

##' @keywords internal
setMethod("partial",
          #### partial-method ####
          signature("gbLocation"),
          function(x) {
            x@partial
          })

##' @keywords internal
setGeneric("partial<-",
           #### 'partial<-'-generic #### 
           function(x, value) {
             standardGeneric("partial<-")
           })

##' @keywords internal
setReplaceMethod("partial",
                 #### 'partial<-'-method ####
                 signature("gbLocation"),
                 function(x, value) {                   
                   error_msg <- "The 'value' argument should be a matrix, or a vector of length 1 or 2." 
                   if (is.vector(value)) {
                     if (length(value) > 2 )
                       stop(error_msg)
                     value <- matrix(
                       if (nrow(x) == 0) logical() else value,
                       nrow=nrow(x),
                       ncol=2, byrow = TRUE)
                   }
                   if ( !is.matrix( value ) || nrow(value) != nrow(x) || ncol(value) != 2 )
                     stop( error_msg )
                   x@partial <- value
                   return(x)
                 })

##' @export 
setMethod("show",
          #### show-method ####
          signature("gbLocation"),
          function( object ) {
            res <- as(object, "character")
            cat(linebreak(res, FORCE=TRUE), "\n" )
          })

## Extract location information from a genBank base span line
.getLocation <- function(gb_base_span)
{
  # transforms location information in the GenBank format (e.g. 1..23,
  # or complement(join(345..543,567..567)) into a named vector
  # start(1, start2), end(1, end2), length(1, length2), strand (1,-1), length
  strand <- ifelse(grepl("complement", gb_base_span), -1, 1)
  join <- ifelse(grepl("join", gb_base_span), 1, 0)
  order <- ifelse(grepl("order", gb_base_span), 1, 0)
  split_loc <- strsplit(unlist(strsplit(gsub("[^0-9\\.,]+", "", gb_base_span), ",")), "\\.\\.")
  start <- as.numeric(lapply(split_loc, "[", 1))
  end <- as.numeric(lapply(split_loc, "[", 2))
  end[is.na(end)] <- start[is.na(end)]
  length <- end - start + 1
  loc <- c(start=start, end=end, length=length, strand=strand, join=join, order=order)
  loc
}

# test cases
# gb_base_span  <- "340"
# gb_base_span <- "340..565"
# gb_base_span <- "<340..565"
# gb_base_span <- "566..>567"
# gb_base_span <- "102.110"
# gb_base_span <- "123^124"
# gb_base_span <- "J00194.1:100..202"
# gb_base_span <- "complement(565..>567)"
# gb_base_span <- "join(345..543,567..>590)"
# gb_base_span <- "order(<345..543,<567..>569,666..7000)"
# gb_base_span <- "join(complement(4918..5163),complement(2691..4571),7665..7899)"
# gb_base_span <- "complement(join(345..543,AL121804.2:567..>569,AL121804.2:<600..603))"

.getLocationS4 <- function(gb_base_span)
{                       
  # single location
  sil <- "^\\d+$"
  
  # between location
  bl <- "\\d+\\^\\d+"
  # within location
  wl <- "[<]?\\d+\\.[>]?\\d+"
  # paired location
  pl <- "[<]?\\d+\\.\\.[>]?\\d+"
  
  # simple location
  sl <- sprintf("([a-zA-z][a-zA-Z0-9]*(\\.[a-zA-Z0-9]+)?\\:)?(%s|%s|%s)", bl, wl, pl)
  # complemented simple location
  csl <- sprintf("complement\\(%s\\)", sl)
  # possibly complemented simplex location
  pcsl <- sprintf("(%s|%s)", sl, csl)
  
  
  # remote accession
  ra <- "([a-zA-Z][a-zA-Z0-9]*(\\.[a-zA-Z0-9]+)?)" 
  # compound location
  cl <- sprintf("(join|order)\\(%s(,%s)*\\)", pcsl, pcsl)
  # complemented compound location
  ccl <- sprintf("complement\\(%s\\)", cl)
  
  .parseSimpleSpan <- function (base_span) {  ## test for strand
    strand <- ifelse(grepl(csl, base_span), -1L, 1L)
    ## get span string
    span_str <- str_extract(base_span,  sl)
    ## get remote accession number
    accn <- str_extract(span_str, ra)
    remote <- ifelse(!is.na(accn), TRUE, FALSE)
    ## get closed and span
    span <- gsub(paste0(ra, "\\:"), "", span_str)
    closed <- ifelse(grepl(wl, span), FALSE, TRUE)
    span <- do.call(rbind, strsplit(span, "\\.\\.|\\.|\\^"))
    ## get partial
    partial <- matrix(grepl("^(<|>)", span), ncol=2)
    span <- matrix(as.integer(gsub("^(<|>)", "", span)), ncol=2)
    return(list(span=span, strand=strand, partial=partial, accn=accn,
                remote=remote, closed=closed))
  }
  
  # test for possibly complemented simple location first
  if (str_detect(gb_base_span, sprintf("^%s$", pcsl))) {
    l <- .parseSimpleSpan(gb_base_span)
    return(.gbLocation(.Data=l$span, strand=l$strand,
                       compound=NA_character_, partial=l$partial,
                       accession=l$accn, remote=l$remote,
                       closed=l$closed))
  }
  # test for possibly complemented compound
  else if ( str_detect(gb_base_span, cl) ) {
    ## test for complementary strand
    strand <- ifelse(grepl(ccl, gb_base_span), -1L, 1L)
    ## get compound
    cmpnd_str <- str_extract(gb_base_span, cl)
    compound <- str_extract(cmpnd_str, "(join|order)")
    ## get span strings
    span_str <- strsplit(str_extract(cmpnd_str, sprintf("%s(,%s)*", pcsl, pcsl)), ",")[[1L]]
    l <- lapply(span_str, .parseSimpleSpan)
    
    if (any(vapply(l, "[[", "strand", FUN.VALUE=integer(1)) < 1L)) {
      strand <- vapply(l, "[[", "strand", FUN.VALUE=integer(1))
    }
    
    return(.gbLocation(.Data=do.call(rbind, lapply(l, "[[", "span")), 
                       strand=strand, compound=compound,
                       partial=do.call(rbind, lapply(l, "[[", "partial")),
                       accession=vapply(l, "[[", "accn", FUN.VALUE=character(1)),
                       remote=vapply(l, "[[", "remote", FUN.VALUE=logical(1))))
  }
  else if ( str_detect(gb_base_span, sil) ) {
    return( .gbLocation(.Data=as.integer(str_extract( gb_base_span, sil ))) )
  }
}
