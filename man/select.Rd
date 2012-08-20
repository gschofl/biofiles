\name{select}
\alias{select}
\alias{select-methods}
\title{Select elements from a GenBank Record}
\usage{
  select(x, index = NULL, keys = "", cols = "")
}
\arguments{
  \item{x}{A \sQuote{\code{gbRecord}} or
  \sQuote{\code{gbFeatureList}} object}

  \item{index}{(Optional) A numeric vector of feature
  indices.}

  \item{keys}{(Optional) A character string
  \sQuote{\emph{keys=value}}-pairs specifying the elements
  to select. See Details.}

  \item{cols}{(Optional) A character string of
  \sQuote{\emph{keys}} that indicate the data to be
  retrieved from the selected features. Supported
  \sQuote{\emph{keys}} are \dQuote{idx} or \dQuote{index},
  \dQuote{location} or \dQuote{range}, \dQuote{start},
  \dQuote{end}, \dQuote{strand}, \dQuote{key}, or any
  qualifier tag (e.g., \dQuote{locus_tag},
  \dQuote{product}). If no \sQuote{\emph{keys}} is given a
  \code{\link{gbFeature-class}} or
  \code{\link{gbFeatureList-class}} object is returned.}
}
\value{
  Depending on the value of \code{select}.
}
\description{
  This function extracts features or information contained
  in features from \code{gbRecord}s or
  \code{gbFeatureList}s.
}
\details{
  \sQuote{keys} specifications are given in a character
  string separated by semicolons.

  \describe{ \item{by index}{ For example:
  \code{"idx=1,3,4"}, \code{"idx=1:4"},
  \code{"index=1,3,5:8"}, \code{"idx=1;idx=3;idx=4"} }
  \item{by location}{ For example:
  \code{"loc=10000:20000"}, \code{"loc=:20000"},
  \code{"location=:20000,40000:60000"} } \item{by key}{ For
  example: \code{"key=CDS"}, \code{"key=CDS,gene"} }
  \item{by qualifier values}{ For example:
  \code{"product=replication"}, \code{"pseudo"} } \item{by
  any combination of the above}{ For example:
  \code{"location=:10000;key=CDS,gene"} } }
}

