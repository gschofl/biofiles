\name{flatten}
\alias{flatten}
\title{Flatten (Nested) Lists.}
\usage{
  flatten(x, start_after = NULL, stop_at = NULL,
    delim_path = ".", do_warn = TRUE, ...)
}
\arguments{
  \item{x}{An arbitrarily deeply nested \code{list}}

  \item{start_after}{An \code{integer} specifying the layer
  after which to start the flattening. \code{NULL} means to
  start at the very top.}

  \item{stop_at}{An \code{integer} specifying the layer at
  which to stop the flattening. \code{NULL} means there is
  not stop criterion.}

  \item{delim_path}{A \code{character} specifying how the
  names of the resulting flattened list should be pasted.}

  \item{...}{Further args.}
}
\value{
  A named \code{list} that features the desired degree of
  flattening.
}
\description{
  Flatten \code{lists} according to specifications made via
  \code{start_after} and/or \code{stop_at}. When keeping
  the defaults, the function will traverse \code{src} to
  retrieve the values at the respective bottom
  layers/bottom elements. These values are arranged in a
  named \code{list} where the respective names can be
  interpreted as the the paths to the retrieved values.
}
\examples{
##
}
\author{
  Janko Thyson \email{janko.thyson.rstuff@googlemail.com}
}
\keyword{internal}

