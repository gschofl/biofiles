\name{gbFeature-class}
\alias{[[,gbFeature-method}
\alias{$,gbFeature-method}
\alias{dbXref,gbFeature-method}
\alias{end,gbFeature-method}
\alias{gbFeature-class}
\alias{getIndex,gbFeature-method}
\alias{getKey,gbFeature-method}
\alias{getLocation,gbFeature-method}
\alias{getQualifier,gbFeature-method}
\alias{getSequence,gbFeature-method}
\alias{hasKey,gbFeature-method}
\alias{hasQualifier,gbFeature-method}
\alias{partial,gbFeature-method}
\alias{range,gbFeature-method}
\alias{show,gbFeature-method}
\alias{start,gbFeature-method}
\alias{strand,gbFeature-method}
\alias{width,gbFeature-method}
\title{gbFeature class}
\usage{
  .gbFeature(...)
}
\arguments{
  \item{...}{Slots of gbFeature}
}
\description{
  \sQuote{gbFeature} is an S4 class that extends the
  \code{\link{gbFeatureList-class}}. This class provides a
  container for feature data retrived from GenBank flat
  files.
}
\details{
  \code{gbFeature} provide the following slots:

  \describe{ \item{.Dir}{The path to the database file
  containing the GenBank record the feature is part of.}
  \item{.ACCN}{Accession number of the GenBank record that
  the feature is part of.} \item{.DEF}{The definition line
  (brief description of the sequence) of the GenBank record
  the feature is part of.} \item{.ID}{Identifier
  (sequential index) of the feature in the GenBank record
  the feature is part of.} \item{key}{Feature key (e.g.
  Source, CDS, gene, etc.)} \item{location}{An object of
  \code{\link{gbLocation-class}}} \item{qualifiers}{Named
  character vector. Name attributes correspond to GenBank
  qualifier tags.} }
}

