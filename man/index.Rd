\name{index}
\alias{index}
\alias{index-methods}
\title{Get feature indices from a GenBank record}
\usage{
  index(x, attributes = FALSE)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}} or
  \code{\linkS4class{gbFeatureList}} object.}

  \item{attributes}{Set the \code{accession},
  \code{definition}, \code{database} attributes.}

  \item{...}{Additional arguments passed to methods.}
}
\value{
  A numeric vector of feature indeces
}
\description{
  Get feature indices from a GenBank record
}

