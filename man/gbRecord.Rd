% Built by roxygen3 0.1
\docType{class}
\name{gbRecord-class}
\alias{gbRecord}
\alias{gbRecord-class}
\title{gbRecord}
\usage{
  gbRecord(gb, with_sequence = TRUE)
}
\arguments{
  \item{gb}{A vector of paths to GenBank record files,
an \code{\linkS4class{efetch}} object containing GenBank record(s), or
a \code{textConnection} to a character vector that can be parsed as
a Genbank record.}

  \item{with_sequence}{Include sequence information if avaliable.}
}

\value{
A \code{\linkS4class{gbRecord}} or
\code{\linkS4class{gbRecordList}}
}
\description{
\dQuote{gbRecord} is an S4 class that provides a container for data
parsed from a GenBank record.

\code{gbRecord} or \code{\linkS4class{gbRecordList}} objects can be
construced by parsing GenBank flat files or an \code{\linkS4class{efetch}}
object containing one or more GenBank records.
}
\details{
For a description of the GenBank format see
\url{http://www.ncbi.nlm.nih.gov/collab/FT/}
}
\section{Generics with methods for gbRecord}{
\itemize{
  \item \code{\link{[[}}
  \item \code{\link{$}}
  \item \code{\link{accession}}
  \item \code{\link{definition}}
  \item \code{\link{end}}
  \item \code{\link{features}}
  \item \code{\link{geneid}}
  \item \code{\link{listQualif}}
  \item \code{\link{ranges}}
  \item \code{\link{select}}
  \item \code{\link{seqinfo}}
  \item \code{\link{seqlengths}}
  \item \code{\link{sequence}}
  \item \code{\link{show}}
  \item \code{\link{start}}
  \item \code{\link{strand}}
  \item \code{\link{summary}}
  \item \code{\link{width}}
  \item \code{\link{write.FeatureTable}}
  \item \code{\link{write.GenBank}}
}

}

\keyword{classes}
