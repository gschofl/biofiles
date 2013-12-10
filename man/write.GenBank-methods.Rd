\docType{methods}
\name{write.GenBank}
\alias{write.GenBank}
\title{Write GenBank records or features to file in GenBank format}
\usage{
write.GenBank(x, file, append = FALSE, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbRecord}} instance.}

  \item{file}{A connection or a character string naming the
  file to write to.}

  \item{header}{if \code{FALSE} exclude the Genbank
  header.}

  \item{sequence}{if \code{FALSE} exclude the sequence.}

  \item{append}{if \code{TRUE} the data is appended to the
  connection.}
}
\description{
Write GenBank records or features to file in GenBank format
}
\details{
For a description of the GenBank format see
\url{http://www.ncbi.nlm.nih.gov/collab/FT/}
}
\examples{
###
}

