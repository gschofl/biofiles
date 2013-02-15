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
  \item \linkS4class{AlignedXStringSet0}
  \item \linkS4class{ANY}
  \item \linkS4class{CompressedIRangesList}
  \item \linkS4class{GappedAlignments}
  \item \linkS4class{GappedRanges}
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
  \item \linkS4class{gbLocation}
  \item \linkS4class{GenomicRanges}
  \item \linkS4class{GRangesList}
  \item \linkS4class{IntervalTree}
  \item \linkS4class{PartitioningByEnd}
  \item \linkS4class{PartitioningByWidth}
  \item \linkS4class{RangedData}
  \item \linkS4class{Ranges}
  \item \linkS4class{RangesList}
  \item \linkS4class{Rle}
  \item \linkS4class{SimpleViewsList}
  \item \linkS4class{SummarizedExperiment}
  \item \linkS4class{Views}
}

}

\section{Method signatures for generic end<-}{
\itemize{
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
  \item \linkS4class{gbLocation}
  \item \linkS4class{GenomicRanges}
  \item \linkS4class{GRangesList}
  \item \linkS4class{IRanges}
  \item \linkS4class{RangedData}
  \item \linkS4class{RangesList}
  \item \linkS4class{SummarizedExperiment}
  \item \linkS4class{Views}
}

}

\keyword{methods}
