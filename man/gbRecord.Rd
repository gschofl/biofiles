\name{gbRecord}
\alias{gbRecord}
\title{\code{gbRecord} or \code{\linkS4class{gbRecordList}} objects can be
construced by parsing GenBank flat files or an \code{\linkS4class{efetch}}
object containing one or more GenBank records.}
\usage{
gbRecord(gb, with_sequence = TRUE)
}
\arguments{
  \item{gb}{A vector of paths to GenBank record files, an
  \code{\linkS4class{efetch}} object containing GenBank
  record(s), or a \code{textConnection} to a character
  vector that can be parsed as a Genbank record.}

  \item{with_sequence}{Include sequence information if
  avaliable.}
}
\value{
A \code{\linkS4class{gbRecord}} or
\code{\linkS4class{gbRecordList}}
}
\description{
\code{gbRecord} or \code{\linkS4class{gbRecordList}}
objects can be construced by parsing GenBank flat files or
an \code{\linkS4class{efetch}} object containing one or
more GenBank records.
}
\details{
For a description of the GenBank format see
\url{http://www.ncbi.nlm.nih.gov/collab/FT/}
}

