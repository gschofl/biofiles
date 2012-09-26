% Built by roxygen3 0.1
\docType{class}
\name{gbRecord-class}
\alias{gbRecord}
\alias{gbRecord-class}
\title{gbRecord}
\usage{
  gbRecord(gb, with_sequence = TRUE, force = FALSE)
}
\arguments{
  \item{gb}{Path to a valid \code{gbRecord} database, a GenBank flat file or
an \code{\linkS4class{efetch}} object containing GenBank record(s).}

  \item{with_sequence}{Parse with sequence information if avaliable.}

  \item{force}{Overwrite existing database directories without prompting.}
}

\value{
A (list of) \code{\linkS4class{gbRecord}} object(s).
}
\description{
\dQuote{gbRecord} is an S4 class that provides a container for data
parsed from a GenBank record. It is implemented as a
\code{\linkS4class{filehashRDS}} database.

\code{gbRecord} instances can be construced by parsing a GenBank flat
file or an \code{\linkS4class{efetch}} instance containing one or more
GenBank records.
If \code{gb} points to a valid \code{gbRecord} database, a \code{gbRecord}
object is initialised in the global environment.
}
\details{
For a description of the GenBank format see
\url{http://www.ncbi.nlm.nih.gov/collab/FT/}
}
\section{Class hierarchy for gbRecord}{
Super classes:
\itemize{
  \item \linkS4class{filehashRDS}
  \item \linkS4class{filehash}
}

}

\section{Generics with methods for gbRecord}{
\itemize{
  \item \code{\link{[[}}
  \item \code{\link{$}}
  \item \code{\link{accession}}
  \item \code{\link{definition}}
  \item \code{\link{features}}
  \item \code{\link{initialize}}
  \item \code{\link{revcomp}}
  \item \code{\link{select}}
  \item \code{\link{sequence}}
  \item \code{\link{show}}
  \item \code{\link{write}}
}

}

\keyword{classes}
