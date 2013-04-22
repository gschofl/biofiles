% Built by roxygen3 0.1
\name{shift}
\alias{shift}
\title{Shift location of features in a GenBank record}
\usage{
  shift(x, shift=0L, split=FALSE, order=FALSE)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeatureList}} or
\code{\linkS4class{gbRecord}} instance (gbFeatureLists must 
be complete and include a 'source' field).}

  \item{shift}{Number of basepairs (or aa residues) to shift.}

  \item{split}{Split features that after the shift extends across the end of
the sequence.}

  \item{order}{Reorder features after the shift.}
}

\value{
A \code{\linkS4class{gbFeatureList}} object.
}
\description{
Shift location of features in a GenBank record
}
\section{Method signatures for generic shift}{
\itemize{
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
  \item \linkS4class{gbLocation}
  \item \linkS4class{gbRecord}
}

}

\keyword{methods}
