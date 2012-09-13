\name{select}
\alias{select}
\alias{select-methods}
\title{Select features from a GenBank Record}
\usage{
  select(x, ..., keys = "", cols = "")
}
\arguments{
  \item{x}{A \sQuote{\code{gbRecord}} or
  \sQuote{\code{gbFeatureList}} object}

  \item{...}{Named values that specify the features to
  select. These are merged with the values of \code{keys}
  to create the actual query. See Details.}

  \item{keys}{A character string composed of
  \sQuote{\emph{key=value}}-pairs separated by semicolons
  that specify the elements to select. See Details.}

  \item{cols}{A character string of \sQuote{\emph{keys}}
  that indicate the data to be retrieved from the selected
  features. Supported \sQuote{\emph{keys}} are \dQuote{idx}
  or \dQuote{index}, \dQuote{location} or \dQuote{range},
  \dQuote{start}, \dQuote{end}, \dQuote{strand},
  \dQuote{key}, or any qualifier tag (e.g.,
  \dQuote{locus_tag}, \dQuote{product}). If no
  \sQuote{\emph{keys}} is given a
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
  Queries can be are provided as named values using
  predefined keywords and \dQuote{qualifier=value} pairs:

  Possible keywords are:

  \describe{ \item{index} or \item{idx}{ For example:
  \code{idx=c(3,4,5,6)}, \code{idx=100:150},
  \code{index=c(1,12:20)} } \item{range}{ For example:
  \code{range="10000..25000"},
  \code{range="..10000,20000..25000"},
  \code{range="30000.."} } \item{key}{ For example:
  \code{key="CDS"}, \code{key=c("CDS", "gene")},
  \code{key="CDS|gene"} } \item{arbitrary qualifiers}{ For
  example: \code{product="ribosomal"},
  \code{locus_tag=c("CPSIT_0123", "CPSIT_0124",
  "CPSIT_0125")}, \code{pseudo=TRUE} } }

  Alternatively (or additionally) queries can be described
  in a single character string of \dQuote{tag=value} pairs
  passed to the \code{keys} argument.

  Different \dQuote{tag=value} pairs are separated by
  semicolons.

  \describe{ \item{index}{ For example: \code{"idx=1,3,4"},
  \code{"idx=1:4"}, \code{"index=1,3,5:8"},
  \code{"idx=1;idx=3;idx=4"} } \item{by range}{ For
  example: \code{"range=10000..20000"},
  \code{"range=..20000"},
  \code{"range=..20000,40000..60000"} } \item{by key}{ For
  example: \code{"key=CDS"}, \code{"key=CDS,gene"} }
  \item{by qualifier values}{ For example:
  \code{"product=replication"}, \code{"pseudo"} } \item{by
  any combination of the above}{ For example:
  \code{"range=..10000;key=CDS,gene"} } }
}

