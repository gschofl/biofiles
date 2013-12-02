\name{ranges}
\alias{ranges}
\alias{ranges<-}
\title{Get or set the range of genomic features}
\usage{
ranges(x, join = FALSE, include = "none", exclude = "")

ranges(x, ...) <- value
}
\arguments{
  \item{x}{A \code{gbFeature} or \code{gbFeatureList}
  object.}

  \item{join}{Join compound genomic locations into a single
  range.}

  \item{key}{Include feature keys with ranges.}

  \item{include}{Include qualifiers as metadata columns.
  Can be "none", "all", or a character vector of qualifier
  tags.}

  \item{exclude}{Exclude specific qualifiers.}

  \item{...}{Further arguments passed to methods.}
}
\value{
A \code{\linkS4class{GRanges}} object.
}
\description{
Get or set the range of genomic features
}

