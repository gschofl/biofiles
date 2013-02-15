% Built by roxygen3 0.1
\name{strand}
\alias{strand}
\alias{strand<-}
\title{Get or set the strand of genomic features}
\usage{
  strand(x, join=FALSE, ...)

  strand(x, ..., value) <- value
}
\arguments{
  \item{x}{A \code{gbFeature} or \code{gbFeatureList}object.}

  \item{join}{Join compound genomic locations into a single range.}

  \item{...}{Further arguments passed to methods.}
}

\value{
An integer vector of 1 (plus strand), -1 (minus strand), or
\code{NA}
}
\description{
Get or set the strand of genomic features
}
\section{Method signatures for generic strand}{
\itemize{
  \item \linkS4class{character}
  \item \linkS4class{DataTable}
  \item factor (S3)
  \item \linkS4class{GappedAlignmentPairs}
  \item \linkS4class{GappedAlignments}
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
  \item \linkS4class{gbLocation}
  \item \linkS4class{GRanges}
  \item \linkS4class{GRangesList}
  \item \linkS4class{integer}
  \item \linkS4class{logical}
  \item \linkS4class{missing}
  \item \linkS4class{Rle}
  \item \linkS4class{SummarizedExperiment}
}

}

\section{Method signatures for generic strand<-}{
\itemize{
  \item \linkS4class{DataTable}
  \item \linkS4class{GappedAlignments}
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
  \item \linkS4class{gbLocation}
  \item \linkS4class{GenomicRanges}
  \item \linkS4class{GRangesList}
  \item \linkS4class{SummarizedExperiment}
}

}

\keyword{methods}
