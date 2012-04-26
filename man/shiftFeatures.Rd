\docType{methods}
\name{shiftFeatures}
\alias{shiftFeatures}
\title{shift genomic location of features in a GenBank record}
\arguments{
  \item{x}{A gbRecord or complete gbFeatureList object
  (including the 'source' field)}

  \item{shift}{Number of basepairs (or aa residues) to
  shift}

  \item{update_db}{Update filehas database with new feature
  locations.}
}
\value{
  A gbFeatureList object
}
\description{
  shift genomic location of features in a GenBank record
}

