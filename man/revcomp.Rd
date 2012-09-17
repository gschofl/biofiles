\name{revcomp}
\alias{revcomp}
\alias{revcomp-methods}
\title{reverse complement features in a GenBank record}
\usage{
  revcomp(x, updateDb=TRUE)
}
\arguments{
  \item{x}{A gbFeatureList or gbRecord object
  (gbFeatureLists must include a 'source' field).}

  \item{order}{Should the resulting gbFeatureList be
  reordered.}

  \item{update_db}{Should the sequence and the new feature
  locations be updated in the underlying filehash
  database?}
}
\value{
  A gbFeatureList object
}
\description{
  reverse complement features in a GenBank record
}

