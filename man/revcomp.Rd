\name{revcomp}
\alias{revcomp}
\alias{revcomp-methods}
\title{reverse complement features in a GenBank record}
\usage{
  revcomp(x, order = FALSE, updateDb = TRUE, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeatureList}} or
  \code{\linkS4class{gbRecord}} instance
  (\code{gbFeatureList} instances must include a 'source'
  field).}

  \item{order}{Should the resulting gbFeatureList be
  reordered.}

  \item{update_db}{Should the sequence and the new feature
  locations be updated in the underlying filehash
  database?}

  \item{...}{Additional arguments passed to methods.}
}
\value{
  A \code{\linkS4class{gbFeatureList}} instance
}
\description{
  reverse complement features in a GenBank record
}

