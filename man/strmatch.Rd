\name{strmatch}
\alias{strmatch}
\title{Extract matched group(s) from a string.}
\usage{
  strmatch(pattern, str, capture = TRUE, perl = TRUE,
    global = TRUE, ignore.case = FALSE)
}
\arguments{
  \item{pattern}{character string containing a regular
  expression}

  \item{str}{character vector where matches are sought}

  \item{capture}{if \code{TRUE} capture groups are returned
  in addition to the complete match}

  \item{perl}{if \code{TRUE} perl-compatible regexps are
  used.}

  \item{global}{if \code{TRUE} \code{gregexpr} is used for
  matching otherwise \code{regexpr}.}

  \item{ignore.case}{case sensitive matching}
}
\value{
  a list containing a \code{match} and a \code{capture}
  component
}
\description{
  Extract matched group(s) from a string.
}
\examples{
##
}
\keyword{character}
\keyword{internal}

