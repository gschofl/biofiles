\name{gbFeatureList-class}
\alias{[,gbFeatureList-method}
\alias{dbxref,gbFeatureList-method}
\alias{end,gbFeatureList-method}
\alias{end<-,gbFeatureList-method}
\alias{gbFeatureList-class}
\alias{hasKey,gbFeatureList-method}
\alias{hasQualif,gbFeatureList-method}
\alias{index,gbFeatureList-method}
\alias{key,gbFeatureList-method}
\alias{location,gbFeatureList-method}
\alias{qualif,gbFeatureList-method}
\alias{range,gbFeatureList-method}
\alias{select,select-method,gbFeatureList-method}
\alias{sequence,gbFeatureList-method}
\alias{show,gbFeatureList-method}
\alias{start,gbFeatureList-method}
\alias{start<-,gbFeatureList-method}
\alias{strand,gbFeatureList-method}
\alias{strand<-,gbFeatureList-method}
\alias{summary,gbFeatureList-method}
\alias{view,view-method,gbFeatureList-method}
\alias{width,gbFeatureList-method}
\title{gbFeatureList class}
\usage{
  .gbFeatureList(...)
}
\arguments{
  \item{...}{Slots of gbFeature}
}
\description{
  gbFeatureList is an S4 class that extends the
  \code{\link{gbRecord-class}}. This class provides a
  container for feature data retrived from GenBank flat
  files.
}
\details{
  gbFeatureList objects have four slots in addition to the
  slots provided by \code{\link{gbRecord-class}} objects:
  \describe{ \item{.Dir}{The path to the database file
  containing the GenBank record the feature list is part
  of.} \item{.ACCN}{Accession number of the GenBank record
  that the feature list is part of.} \item{.DEF}{The
  definition line (brief description of the sequence) of
  the GenBank record the feature list is part of.}
  \item{.Data}{A list of gbFeature objects} }
}

