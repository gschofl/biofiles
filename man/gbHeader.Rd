\name{.gbHeader}
\alias{.gbHeader}
\title{Generator object for the \code{\linkS4class{gbHeader}} class}
\usage{
.gbHeader(...)
}
\arguments{
  \item{...}{List of arguments; must be named arguments
  corresponding to the fields of a
  \code{\linkS4class{gbHeader}} object}
}
\description{
The generator object for the \code{\linkS4class{gbHeader}}
reference class.
}
\section{Methods}{
  \describe{ \item{\code{#new(...)}:}{ Create a new
  \code{\linkS4class{gbHeader}} object. }
  \item{\code{#to_string(write_to_file = FALSE)}:}{
  Generate a character string representation of a GenBank
  file header. } \item{\code{#write(file = "", append =
  FALSE, sep = "\n"}:}{ Write a GenBank header to file. } }
}
\seealso{
\code{"\linkS4class{gbHeader}"}
}
\keyword{classes}
\keyword{internal}

