\name{shift}
\alias{shift}
\alias{shift-methods}
\title{shift location of features in a GenBank record}
\usage{
  shift(x, shift=0L, split=FALSE, order=FALSE,
    update_db=FALSE)
}
\arguments{
  \item{x}{A gbLocation, gbFeature, gbFeatureList, or
  gbRecord object (gbFeatureLists must include a 'source'
  field).}

  \item{shift}{Number of basepairs (or aa residues) to
  shift.}

  \item{split}{(For gbFeatureList and gbRecord objects)
  Should a feature that spans across the end of the
  sequence be split.}

  \item{order}{(For gbFeatureList and gbRecord objects)
  Should the resulting gbFeatureList be reordered.}

  \item{update_db}{Should filehash database be updated with
  new feature locations.}
}
\value{
  A \code{\link{gbLocation-class}},
  \code{\link{gbFeature-class}},
  or\code{\link{gbFeatureList-class}} object
}
\description{
  shift location of features in a GenBank record
}

