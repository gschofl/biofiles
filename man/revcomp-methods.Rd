\docType{methods}
\name{revcomp}
\alias{revcomp}
\alias{revcomp,gbFeatureList-method}
\alias{revcomp,gbRecord-method}
\title{Reverse-complement features in a GenBank record}
\usage{
revcomp(x, order=TRUE)

\S4method{revcomp}{gbFeatureList}(x, order = TRUE)

\S4method{revcomp}{gbRecord}(x, order = TRUE)
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

