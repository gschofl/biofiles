#' Flatten (Nested) Lists or Environments.
#'
#' Flatten \code{lists} or \code{environments} according to specifications
#' made via arg \code{start.after} and/or arg \code{stop.at}. When keeping 
#' the defaults, the function will traverse arg \code{src} (if \code{src} is 
#' an \code{environment}, it is coerced to a \code{list} 
#' via \code{\link{envirAsList}} first) to retrieve the values at the 
#' respective bottom layers/bottom elements. These values are arranged in a 
#' named \code{list} where the respective names can be interpreted as the
#' the paths to the retrieved values. See examples.   
#'
#' @param src A named (arbitrarily deeply nested) \code{list} or an 
#' \code{environment} that should be flattened.
#' @param start.after An \code{integer} specifying the layer after which to 
#' start the flattening. \code{NULL} means to start at the very top. See
#' examples.
#' @param stop.at An \code{integer} specifying the layer at which to stop
#' the flattening. \code{NULL} means there is not stop criterion.
#' @param delim.path A \code{character} (length: 1) specifying how the names
#' of the resulting flattened list should be pasted.
#' @param .do.debug If \code{TRUE}, print information that might be helpful
#' for debugging.
#' @param ... Further args.
#' @return A named \code{list} that features the desired degree of flattening.
#' @export
#' @author Janko Thyson \email{janko.thyson.rstuff@@googlemail.com}
#' @seealso \code{\link{envirAsList}}
#' @example inst/examples/flatten.r
flatten <- function(
                    src, 
                    start.after=NULL, 
                    stop.at=NULL, 
                    delim.path="/",
                    do.warn=TRUE,
                    .do.debug=FALSE,
                    ...
                    ){
  #---------------------------------------------------------------------------
  # VALIDATE
  #---------------------------------------------------------------------------
  
  if(!is.list(src) & !is.environment(src)){
    stop("Arg 'src' must be a 'list' or an 'environment'.")
  }
  if(!is.null(start.after) & !is.null(stop.at)){
    if(start.after == 1& stop.at == 1){
      msg <- c(
        "Invalid specification:",
        paste("* start.after: ", start.after, sep=""),
        paste("* stop.at:     ", stop.at, sep="")
        )
      stop(cat(msg, sep="\n"))
    }
  }
  # /VALIDATE ----------
  
  #---------------------------------------------------------------------------
  # INNER FUNCTIONS
  #---------------------------------------------------------------------------
  
  .startAfterInner <- function(
                                envir,
                                nms,
                                out.1,
                                do.reset=FALSE,
                                ...
                                ){
    .do.debug <- envir$.do.debug
    idx.diff <- diff(c(envir$start.after, length(envir$counter)))
    if(.do.debug){
      cat(c("", "+++", ""), sep="\n")
      #                print("+++")
      cat("names:", sep="\n")
      print(names(out.1))
      cat("envir$counter:", sep="\n")
      print(envir$counter)
      cat("idx.diff:", sep="\n")
      print(idx.diff)            
    }
    # UPDATE IF DEGREE OF NESTEDNESS EXCEEDS START CRITERION
    if(idx.diff > 0){
      idx.cutoff      <- (
        length(envir$counter)-idx.diff+1):length(envir$counter
                                                 ) 
      idx.left        <- envir$counter[-idx.cutoff]
      nms.1           <- nms[idx.cutoff]
      names(out.1)    <- paste(nms.1, collapse="/") 
      # UPDATE SRC
      idx.append <- sapply(envir$history, function(x.hist){
        all(idx.left == x.hist)        
      })
      if(.do.debug){
        cat("idx.cutoff:", sep="\n")
        print(idx.cutoff)
        cat("idx.left:", sep="\n")
        print(idx.left)
        cat("idx.append:", sep="\n")
        print(idx.append)
        cat("names remaining:", sep="\n")
        print(names(out.1))
      }
      if(any(idx.append)){                                          
        envir$src[[idx.left]] <- append(envir$src[[idx.left]], 
                                        values=out.1)                    
      } else {
        envir$src[[idx.left]] <- out.1
        # UPDATE HISTORY
        envir$history <- c(envir$history, list(idx.left))
      }
      envir$out <- envir$src
      # /             
    } 
    if(idx.diff < 0){
      envir$out <- envir$src
      #            if(envir$do.warn & !envir$do.block.warning){
      #                warning(paste("Argument 'start.after=", envir$start.after, 
      #                    "' exceeds maximum degree of nestedness (=", 
      #                    envir$start.after + idx.diff, ").", sep=""))
      #                envir$do.block.warning <- TRUE
      #            }
    }
    # /
    # RESET
    if(do.reset){
      envir$nms       <- envir$nms[-length(envir$nms)]
      envir$counter   <- envir$counter[-length(envir$counter)]
    }
    # /
    return(TRUE)
  }
  
  .updateOutInner <- function(
                              envir,
                              out.1,
                              do.reset=FALSE,
                              ...
                              ){
    .do.debug <- envir$.do.debug
    # UPDATE OUT
    out.0           <- get("out", envir = envir)
    out             <- c(out.0, out.1)
    envir$out       <- out
    # /
    # RESET
    if(do.reset){
      envir$nms       <- envir$nms[-length(envir$nms)]
      envir$counter   <- envir$counter[-length(envir$counter)]
    }
    # /
    return(TRUE)
  }
  
  .flattenInner <- function(
                            x, 
                            envir, 
                            ...
                            ){
    .do.debug <- envir$.do.debug
    if( (class(x)=="list" & length(x) != 0) |
      (class(x) == "environment" & length(x) != 0)
        ){
      if(class(x) == "environment"){
        x <- as.list(x)
      }
      # UPDATE
      envir$counter.history <- c(envir$counter.history, list(envir$counter))
      # EXIT IF DEGREE EXCEEDS CUTOFF
      if(!is.null(envir$stop.at)){
        if(length(envir$counter) > envir$stop.at){ 
          # THIS
          nms             <- get("nms", envir=envir)
          if(.do.debug){
            cat("names:", sep="\n")
            print(paste(nms, collapse=envir$delim.path))
          }
          out.1           <- list(x)
          names(out.1)    <- paste(nms, collapse=envir$delim.path)
          # /
          # DECISION ON FLATTENING
          if(!is.null(envir$start.after)){
            .startAfterInner(envir=envir, nms=nms, out.1=out.1, 
                             do.reset=TRUE)
            return(NULL)
            #                    }
            # /
        } else {
          .updateOutInner(envir=envir, out.1=out.1, do.reset=TRUE)
          return(NULL)
        }
      }
    }
      # /
      # LOOP OVER ELEMENTS
      for(i in seq(along=x)){
        # UPDATE COUNTER
        envir$counter <- c(envir$counter, i)
        # UPDATE NAMES
        assign("nms", c(get("nms", envir=envir), names(x[i])), envir=envir)
        # RECURSIVE FLATTENING
        .flattenInner(x[[i]], envir) # call  recursively
        # RESET COUNTER
        if(i == length(x)){
          envir$nms       <- envir$nms[-length(envir$nms)]
          envir$counter   <- envir$counter[-length(envir$counter)]
        }
        # /
      }
      # /
  } else {
    # THIS
    nms             <- get("nms", envir=envir)
    if(.do.debug){
      cat("names:", sep="\n")
      print(paste(nms, collapse=envir$delim.path))
    }
    out.1           <- list(x)
    names(out.1)    <- paste(nms, collapse=envir$delim.path)
    # /
    # DECISION ON FLATTENING
    if(!is.null(envir$start.after)){
      .startAfterInner(envir=envir, nms=nms, out.1=out.1)
    } else {
      .updateOutInner(envir=envir, out.1=out.1)
    }
    if(.do.debug){
      cat("out.1:", sep="\n")
      print(out.1)
    }
    # RESET
    envir$nms       <- envir$nms[-length(envir$nms)]
    envir$counter   <- envir$counter[-length(envir$counter)]
    # /
  }
    return(TRUE)
}
  
  # /INNER FUNCTIONS ----------
  
  #---------------------------------------------------------------------------
  # ACTUAL PROCESSING
  #---------------------------------------------------------------------------
  
  # COERCE TO LIST
  if(class(src) == "environment"){
    src <- envirAsList(src=src)
  }
  # /
  # PRESERVE ORIGINAL (just in case)
  src.0               <- src
  out                 <- list()
  # ENVIR
  envir               <- new.env()
  envir$.do.debug     <- .do.debug
  envir$counter       <- NULL
  envir$counter.history <- NULL
  envir$delim.path    <- delim.path
  envir$do.warn       <- do.warn
  envir$do.block.warning    <- FALSE
  envir$history       <- NULL
  envir$nms           <- NULL
  envir$out           <- list()
  envir$src           <- src
  envir$start.after   <- start.after
  if(!is.null(stop.at)){
    stop.at.0 <- stop.at
    if(stop.at == 1){
      return(src)
    } else {
      stop.at <- stop.at - 1
    }
  }
  envir$stop.at       <- stop.at
  # /
  # APPLY INNER
  .flattenInner(src, envir)
  
  if(envir$do.warn){
    max.length <- max(sapply(envir$counter.history, function(x){
      length(x)        
    }))
    #        if(!envir$do.block.warning){
    if(!is.null(start.after)){            
      if(start.after > max.length){                        
        warning(paste("Argument 'start.after=", start.after, 
                      "' exceeds maximum degree of sublayer nestedness (=", 
                      max.length, ").", sep=""))
      }
    }
    if(!is.null(stop.at)){
      if(stop.at.0 > max.length){
        warning(paste("Argument 'stop.at=", stop.at.0, 
                      "' exceeds maximum degree of sublayer nestedness (=", 
                      max.length, ").", sep=""))    
      }
    }
    }
  
  out <- envir$out
  
  # /ACTUAL PROCESSING ----------
  
  return(out)    
  }


#' Coerce Environment to List (Recursively).
#'
#' Recursively coerces an \code{environment} to a \code{list}.  
#'
#' @param src A an \code{environment} that should be coerced.
#' @param ... Further args.
#' @return A named \code{list} that corresponds to the recursively coerced
#' initial \code{environment}.
#' @export
#' @author Janko Thyson \email{janko.thyson.rstuff@@googlemail.com}
#' @seealso \code{\link{flatten}}
#' @example inst/examples/envirAsList.r
envirAsList <- function(
                        src, 
                        ...
                        ){
  if(class(src) == "environment"){
    envir   <- new.env()  
    src     <- as.list(src)
    # LOOP OVER ELEMENTS
    out <- lapply(seq(along=src), function(x.src){
      envir$names <- c(envir$names, names(src[x.src]))
      # RECURSIVE FLATTENING
      out <- envirAsList(src[[x.src]])
      return(out)
    })
    names(out) <- envir$names
    # /
  } else {
    out <- src 
  }
  return(out)
}









