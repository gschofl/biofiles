\docType{methods}
\name{revcomp}
\alias{revcomp}
\title{Reverse-complement features in a GenBank record}
\usage{
revcomp(x, order=TRUE)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeatureList}} or
  \code{\linkS4class{gbRecord}} object (gbFeatureLists must
  be complete and include a 'source' field).}

  \item{order}{Reorder features after reverse-complementing
  them.}
}
\description{
Reverse-complement features in a GenBank record
}

