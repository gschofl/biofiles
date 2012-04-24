\name{joinLines}
\alias{joinLines}
\title{Join lines of text recursively}
\usage{
  joinLines(lines, extract_pat = ".*", break_pat = NULL,
    sep = TRUE)
}
\arguments{
  \item{lines}{A character vector}

  \item{extract_pat}{A regular expression extracting parts
  of the character strings to join (defaults to \code{.*})}

  \item{break_pat}{A regular expression that sets a break
  condition. The line where this pattern matches is the
  last to be joined. If not set, all lines will be joined}

  \item{sep}{if \code{FALSE} the joined lines are
  concatenated without intervening spaces}
}
\value{
  A list with two elements. The first contains the joined
  character vector, the second the number of lines joined
}
\description{
  Join lines of text recursively
}
\keyword{internal}

