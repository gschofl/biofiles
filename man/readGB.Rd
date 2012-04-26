\name{readGB}
\alias{readGB}
\title{General function for importing GenBank flat files}
\usage{
  readGB(gb, with_sequence=TRUE, force=FALSE)
}
\arguments{
  \item{gb}{Path to a GenBank flat file or an
  \code{\link[rentrez]{efetch-class}} object containing
  GenBank record(s).}

  \item{with_sequence}{If \code{TRUE}, sequence information
  will be included.}

  \item{force}{If \code{TRUE} existing database directories
  are overwritten without prompting.}

  \item{parallel}{Use parallel cores.}
}
\value{
  A (list of) \code{\link{gbRecord-class}} object(s).
}
\description{
  For a description of the GenBank format see
  \url{http://www.ncbi.nlm.nih.gov/collab/FT/}
}

