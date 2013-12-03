\name{.gbReference}
\alias{.gbReference}
\title{Generator object for the \code{\linkS4class{gbReference}} class}
\usage{
.gbReference(...)
}
\arguments{
  \item{...}{List of arguments.}
}
\description{
The generator object for the
\code{\linkS4class{gbReference}} reference class.
}
\section{Methods}{
  \describe{ \item{\code{#new(refline, authors, consrtm,
  title, journal, pubmed, remark)}:}{ Create a new
  \code{\linkS4class{gbReference}} object}
  \item{\code{#to_string(write_to_file = FALSE)}:}{
  Generate a character string representation of a GenBank
  reference. } }
}
\seealso{
\code{\linkS4class{gbReference}}
}
\keyword{classes}
\keyword{internal}

