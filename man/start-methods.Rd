\docType{methods}
\name{start}
\alias{start}
\alias{start<-}
\title{Get or set the start of genomic features}
\usage{
start(x, join = FALSE, ...)

start(x, check = TRUE) <- value
}
\arguments{
  \item{x}{A \code{gbFeature} or \code{gbFeatureList}
  object.}

  \item{join}{Join compound genomic locations into a single
  range.}

  \item{...}{Further arguments passed to methods.}
}
\value{
An integer vector
}
\description{
Get or set the start of genomic features
}

