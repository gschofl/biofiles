\name{location}
\alias{location}
\alias{location-methods}
\title{Get feature locations from a GenBank record}
\usage{
  location(x, attributes = TRUE, join = FALSE)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}} or
  \code{\linkS4class{gbFeatureList}} object.}

  \item{attributes}{Set the \code{accession},
  \code{definition}, \code{database} attributes.}

  \item{join}{combine compound locations}

  \item{...}{Additional arguments passed to methods.}
}
\value{
  A data frame
}
\description{
  Get feature locations from a GenBank record
}

