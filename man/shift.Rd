\docType{methods}
\name{shift}
\alias{shift}
\title{shift location of features in a GenBank record}
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
  A gbLocation, gbFeature, or gbFeatureList object
}
\description{
  shift location of features in a GenBank record
}

