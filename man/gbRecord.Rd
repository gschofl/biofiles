% Built by roxygen3 0.1
\docType{class}
\name{gbRecord-class}
\alias{gbRecord}
\alias{gbRecord-class}
\title{gbRecord class}
\usage{
  gbRecord(db_dir, header, features, sequence = NULL)
}
\description{
gbRecord is an S4 class that provides a container for data parsed from
a GenBank flat file record. It is implemented as an extension of a
\dQuote{\code{\linkS4class{filehashRDS-class}}} with no additional slots.
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
\keyword{internal}
