
# gbLocation-class ----------------------------------------------------

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
##' @importClassesFrom intervals Intervals_full
##'
##' @keywords internal
##' @name gbLocation-class
##' @rdname gbLocation-class
##' @aliases show,gbLocation-method
.gbLocation <- 
  #### gbLocation-class ####
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

##' @keywords internal
setMethod("coerce",
          #### coerce-method, character ####
          signature( from = "gbLocation", to = "character" ),
          function( from, to, strict ) {
            if (nrow(from) == 0)
              return(character())
            else {
              
              cl <- closed(from)
              par <- partial(from)
              str <- from@strand
              cmp <- from@compound
              acc <- from@accession
              rem <- from@remote
              
              span <- ifelse(cl[,1],
                             "..", 
                             ifelse(from[,2] == from[,1] + 1,
                                    "^",
                                    ".")
              )
              
              pos <- ifelse(from[,1] == from[,2],
                            from[,1], 
                            paste0(
                              ifelse( par[,1], "", "<" ),
                              from[,1],
                              span,
                              from[,2],
                              ifelse( par[,2], "", ">" ),
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
            cat("Feature location:\n")
            loc <- as(object, "character")
            cat(linebreak(loc, FORCE=TRUE), "\n" )
          })
