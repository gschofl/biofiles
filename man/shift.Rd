% Built by roxygen3 0.1
\name{shift}
\alias{shift}
\title{Shift location of features in a GenBank record}
\usage{
  shift(x, shift=0L, split=FALSE, order=FALSE,
    updateDb=FALSE)
}
\arguments{
  \item{x}{A gbLocation, gbFeature, gbFeatureList, or gbRecord instance
(gbFeatureLists must include a 'source' field).}

  \item{shift}{Number of basepairs (or aa residues) to shift.}

  \item{split}{(For gbFeatureList and gbRecord objects) Should a feature
that spans across the end of the sequence be split.}

  \item{order}{(For gbFeatureList and gbRecord objects) Should the
resulting gbFeatureList be reordered.}

  \item{updateDb}{Should filehash database be updated with new feature
locations.}
}

\value{
A \code{\link{gbLocation-class}}, \code{\link{gbFeature-class}},
or\code{\link{gbFeatureList-class}} object
}
\description{
Shift location of features in a GenBank record
}
\section{Method signatures for generic shift}{
\itemize{
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
  \item \linkS4class{gbLocation}
  \item \linkS4class{gbRange}
  \item \linkS4class{gbRecord}
}

}

\keyword{methods}
