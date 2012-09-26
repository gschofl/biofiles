.revcomp_features <- function (x, order=FALSE, updateDb=FALSE) {
  
  if (is(x, "gbRecord")) {
    max_len <- x$length
    f <-x$features
  } else if (is(x, "gbFeatureList")) {
    
    if (length(x["source"]) == 0) {
      stop("No source field in gbFeatureList")
    }
    
    max_len <- end(x["source"])
    f <- x
  }
  
  new_end <- max_len - start(f) + 1
  new_start <- max_len - end(f) + 1
  new_strand <- strand(f)*-1
  
  start(f) <- new_start
  end(f) <- new_end
  strand(f) <- new_strand
  
  if (order) {
    f <- f[order(new_start)]
  }
  
  f <- .gbFeatureList(.Data=f, .Dir=f@.Dir, .ACCN=f@.ACCN, .DEF=f@.DEF)
  
  if (updateDb) {
    db <- init_db(f@.Dir)
    dbInsert(db, key="features", value=f)
    
    seq <- dbFetch(db, "sequence")
    new_seq <- reverseComplement(seq)
    dbInsert(db, key="sequence", value=new_seq)
  }
  
  return( f )
}

