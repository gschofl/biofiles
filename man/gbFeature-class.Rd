\docType{methods}
\name{gbFeature-class}
\alias{[[,gbFeature-method}
\alias{$,gbFeature-method}
\alias{gbFeature-class}
\alias{getIndex}
\alias{getIndex,gbFeature-method}
\alias{getKey}
\alias{getKey,gbFeature-method}
\alias{getLocation}
\alias{getLocation,gbFeature-method}
\alias{getQualifier}
\alias{getQualifier,gbFeature-method}
\alias{getSequence}
\alias{getSequence,gbFeature-method}
\alias{hasKey}
\alias{hasKey,gbFeature-method}
\alias{hasQualifier}
\alias{hasQualifier,gbFeature-method}
\alias{show,gbFeature-method}
\title{gbFeature class}
\usage{
  .gbFeature(...)
}
\arguments{
  \item{data}{An instance of \code{\link{gbFeature-class}}
  or \code{\link{gbFeatureList-class}}}
}
\value{
  A numeric vector of feature indeces
}
\description{
  \sQuote{gbFeature} is an S4 class that extends the
  \code{\link{gbFeatureList-class}}. This class provides a
  container for feature data retrived from GenBank flat
  files.

  Retrieve index of a GenBank feature.
}
\details{
  \describe{ \item{.Dir}{The path to the database file
  containing the GenBank record the feature is part of.}
  \item{.ACCN}{Accession number of the GenBank record that
  the feature is part of.} \item{.DEF}{The definition line
  (brief description of the sequence) of the GenBank record
  the feature is part of.} \item{.ID}{Identifier
  (sequential index) of the feature in the GenBank record
  the feature is part of.} \item{key}{Feature key (e.g.
  Source, CDS, gene, etc.)} \item{location}{Named numeric
  vector. Name attributes: 'start<N>', 'end<N>',
  'length<N>', 'strand', 'join', 'order'.}
  \item{qualifiers}{Named character vector. Name attributes
  correspond to GenBank qualifier tags.} }
}
\seealso{
  Other accessors: \code{\link{getFeatures}}
}

