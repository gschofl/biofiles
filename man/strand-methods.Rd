\docType{methods}
\name{strand}
\alias{strand}
\alias{strand<-}
\title{Get or set the strand of genomic features}
\usage{
strand(x, join = FALSE, ...)

strand(x, ...) <- value
}
\arguments{
  \item{x}{A \code{gbFeature} or
  \code{gbFeatureList}object.}

  \item{join}{Join compound genomic locations into a single
  range.}

  \item{...}{Further arguments passed to methods.}
}
\value{
An integer vector of 1 (plus strand), -1 (minus strand), or
\code{NA}
}
\description{
Get or set the strand of genomic features
}

