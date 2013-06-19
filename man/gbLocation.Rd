% Built by roxygen3 0.1
\docType{class}
\name{gbLocation-class}
\alias{gbLocation}
\alias{gbLocation-class}
\title{gbLocation}
\usage{
  gbLocation(gb_base_span, accession = "")
}
\description{
\dQuote{gbLocation} is an S4 class that provides a container for
GenBank feature location descriptors.
}
\details{
For more information see the 
\href{ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt}{GenBank Release Note}
}
\section{Class hierarchy for gbLocation}{
Super classes:
\itemize{
  \item \linkS4class{gbLocationOrNull}
}

}

\section{Generics with methods for gbLocation}{
\itemize{
  \item \code{\link{accession}}
  \item \code{\link{end}}
  \item \code{\link{end<-}}
  \item \code{\link{fuzzy}}
  \item \code{\link{shift}}
  \item \code{\link{show}}
  \item \code{\link{start}}
  \item \code{\link{start<-}}
  \item \code{\link{strand}}
  \item \code{\link{strand<-}}
  \item \code{\link{width}}
}

}

\keyword{classes}
\section{Slots}{
\describe{
  \item{range}{An integer matrix indicating the base numbers delimiting a
sequence span.}

  \item{fuzzy}{A logical matrix indicating fuzzy start and/or end
(e.g. <1..200).}

  \item{strand}{An integer vector containing -1, 1, or NA.}

  \item{compound}{A character code specifying how multiple ranges
are joined. One of \sQuote{join}, \sQuote{order}, or \sQuote{bond}.}

  \item{accession}{A character vector; the accession number of the sequence
of the feature this location lives on.}

  \item{remote}{A logical vector}

  \item{type}{A character vector describing the type of the position. Normally
an "R" for \sQuote{range} (e.g., 1..200 or point position 200) and "B"
for \sQuote{between bases} (e.g., 36^37).}

}
}
