\name{gbLocation-class}
\alias{gbLocation-class}
\alias{show,gbLocation-method}
\title{gbLocation class}
\usage{
  .gbLocation(...)
}
\description{
  gbLocation is a container for GenBank Feature Locations.
  It extends \code{\link[Intervals]{Intervals_full}} and
  has xxx additional Slots: \describe{ \item{strand}{An
  integer code for minus (-1) or plus (1) strand}
  \item{compound}{A character code specifying how multiple
  segments are joined. One of 'join' or 'order'}
  \item{partial}{A logical matrix specifying whether
  residues are missing from the 5' and 3' ends
  respectively.} \item{accession}{} \item{remote}{} }
}
\details{
  For more information see the
  \href{ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt}{GenBank
  Release Note}
}

