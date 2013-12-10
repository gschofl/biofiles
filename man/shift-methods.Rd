\docType{methods}
\name{shift}
\alias{shift}
\alias{shift,gbFeature-method}
\alias{shift,gbFeatureList-method}
\alias{shift,gbRecord-method}
\title{Shift location of features in a GenBank record}
\usage{
shift(x, shift=0L, split=FALSE, order=FALSE)

\S4method{shift}{gbFeature}(x, shift = 0L, ...)

\S4method{shift}{gbFeatureList}(x, shift = 0L, split = FALSE,
  order = FALSE)

\S4method{shift}{gbRecord}(x, shift, split = FALSE, order = FALSE)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeatureList}} or
  \code{\linkS4class{gbRecord}} instance (gbFeatureLists
  must be complete and include a 'source' field).}

  \item{shift}{Number of basepairs (or aa residues) to
  shift.}

  \item{split}{Split features that after the shift extends
  across the end of the sequence.}

  \item{order}{Reorder features after the shift.}
}
\value{
A \code{\linkS4class{gbFeatureList}} object.
}
\description{
Shift location of features in a GenBank record
}

