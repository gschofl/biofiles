% Built by roxygen3 0.1
\docType{methods}
\name{location}
\alias{location}
\alias{location,gbFeature-method}
\alias{location,gbFeatureList-method}
\title{Get genomic locations of features}
\usage{
  location(x, attributes = FALSE, join = FALSE, ...)

  \S4method{location}{gbFeature}(x, attributes = FALSE,
    join = FALSE, ...)

  \S4method{location}{gbFeatureList}(x, attributes = FALSE,
    join = FALSE, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}} instance.}

  \item{attributes}{Include the \code{accession}, \code{definition},
\code{database} attributes of the feature.}

  \item{join}{Join compound genomic locations into a single range.}

  \item{x}{A \code{\linkS4class{gbFeatureList}} instance.}

  \item{attributes}{Include the \code{accession}, \code{definition},
\code{database} attributes of the feature.}

  \item{join}{Join compound genomic locations into a single range.}
}

\value{
A \code{\linkS4class{gbRange}} object including the feature key
and the feature index.

A \code{\linkS4class{gbRange}} object including feature keys
and feature indices.
}
\description{
Get genomic locations of features

Get genomic locations of features
}
\section{Method signatures for generic location}{
\itemize{
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
}

}

\keyword{methods}
