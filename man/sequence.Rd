% Built by roxygen3 0.1
\name{getSequence}
\alias{getSequence}
\title{Retrieve sequence from a GenBank record or a Feature.}
\usage{
  getSequence(x, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbRecord}}, \code{\linkS4class{gbFeature}},
 or \code{\linkS4class{gbFeatureList}} instance.}

  \item{...}{Additional arguments passed to methods.}
}

\value{
An \code{\linkS4class{XStringSet}} object.
}
\description{
Retrieve sequence from a GenBank record or a Feature.
}
\section{Method signatures for generic getSequence}{
\itemize{
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
  \item \linkS4class{gbRecord}
  \item \linkS4class{gbRecordList}
}

}

\keyword{methods}
