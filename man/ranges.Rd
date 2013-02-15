% Built by roxygen3 0.1
\name{ranges}
\alias{ranges}
\alias{ranges<-}
\title{Get or set the range of genomic features}
\usage{
  ranges(x, join=FALSE, with_qual="none", without_qual="")

  ranges(x, ..., value) <- value
}
\arguments{
  \item{x}{A \code{gbFeature} or \code{gbFeatureList} object.}

  \item{join}{Join compound genomic locations into a single range.}

  \item{with_qual}{Include qualifiers as metadata columns. Can be "none",
"all", or a character vector of qualifier tags.}

  \item{without_qual}{Exclude specific qualifiers.}

  \item{...}{Further arguments passed to methods.}
}

\value{
A \code{\linkS4class{GRanges}} object.
}
\description{
Get or set the range of genomic features
}
\keyword{methods}
