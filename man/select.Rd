\docType{methods}
\name{select}
\alias{select}
\title{Select elements from a GenBank Record}
\arguments{
  \item{x}{A \sQuote{\code{gbRecord}} or
  \sQuote{\code{gbFeatureList}} object}

  \item{subset}{A character string specifying the elements
  to select indicated by index, key, location, or qualifier
  value. See Details.}

  \item{select}{A character string specifying the
  information to be returned from the subsetted elements.
  Supported values are \dQuote{idx} or \dQuote{index},
  \dQuote{location} or \dQuote{range}, \dQuote{start},
  \dQuote{end}, \dQuote{strand}, \dQuote{key}, or any
  qualifier tag (e.g., \dQuote{locus_tag},
  \dQuote{product}). The default will return a
  \code{\link{gbFeature-class}} or
  \code{\link{gbFeatureList-class}} object}
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
  \sQuote{subset} specifications are given in a character
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

