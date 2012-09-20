\name{gbLocation-class}
\alias{as.gbLocation,gbLocation-method}
\alias{end,gbLocation-method}
\alias{end<-,gbLocation-method}
\alias{gbLocation-class}
\alias{partial,gbLocation-method}
\alias{range,gbLocation-method}
\alias{show,gbLocation-method}
\alias{start,gbLocation-method}
\alias{start<-,gbLocation-method}
\alias{strand,gbLocation-method}
\alias{strand<-,gbLocation-method}
\alias{width,gbLocation-method}
\title{gbLocation class}
\usage{
  .gbLocation(...)
}
\arguments{
  \item{...}{Slots of gbLocation}
}
\description{
  dQuote{\code{gbLocation}} is a container for GenBank
  Feature Locations. It extends
  \code{\link[intervals]{Intervals_full}} and provides 5
  additional Slots: \describe{ \item{strand}{A single
  integer -1 or 1 (for minus or plus strand, respectively)}
  \item{compound}{A character code specifying how multiple
  segments are joined. One of \sQuote{join} or
  \sQuote{order}.} \item{partial}{A logical matrix
  specifying whether residues are missing from the 5' and
  3' ends respectively.} \item{accession}{} \item{remote}{}
  }
}
\details{
  For more information see the
  \href{ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt}{GenBank
  Release Note}
}

