% Built by roxygen3 0.1
\docType{methods}
\name{strand,gbFeature-method}
\alias{strand,gbFeature-method}
\title{Get the strand of genomic features}
\usage{
  \S4method{strand}{gbFeature}(x, ...)
}
\arguments{
  \item{x}{A \code{gbFeature} object.}

  \item{join}{Join compound genomic locations into a single range.}

  \item{...}{Further arguments passed to methods.}
}

\value{
An integer vector of 1 (plus strand), -1 (minus strand), or
\code{NA}
}
\description{
Get the strand of genomic features
}
\keyword{methods}
