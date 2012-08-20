\name{getLocation}
\alias{getLocation}
\alias{getLocation-methods}
\title{Retrieve feature locations from a GenBank record}
\usage{
  getLocation(x, attributes = TRUE, join = FALSE)
}
\arguments{
  \item{x}{A \code{\link{gbFeature}} or
  \code{\link{gbFeatureList}} object}

  \item{attributes}{Set the \code{accession},
  \code{definition}, \code{database} attributes.}

  \item{join}{combine compound locations}
}
\value{
  A data frame
}
\description{
  Retrieve feature locations from a GenBank record
}

