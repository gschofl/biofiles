% Built by roxygen3 0.1
\name{end}
\alias{end}
\alias{end<-}
\title{Get or set the end of genomic features}
\usage{
  start(x, join=FALSE, ...)

  end(x, check = TRUE, value) <- value
}
\arguments{
  \item{x}{A \code{gbFeature} or \code{gbFeatureList} object.}

  \item{join}{Join compound genomic locations into a single range.}

  \item{...}{Further arguments passed to methods.}
}

\value{
An integer vector
}
\description{
Get or set the end of genomic features
}
\section{Method signatures for generic end}{
\itemize{
  \item \linkS4class{ANY}
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
  \item \linkS4class{gbLocation}
  \item \linkS4class{gbRecord}
  \item \linkS4class{gbRecordList}
}

}

\section{Method signatures for generic end<-}{
\itemize{
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
  \item \linkS4class{gbLocation}
}

}

\keyword{methods}
