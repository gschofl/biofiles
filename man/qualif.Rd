\docType{methods}
\name{qualif}
\alias{locusTag}
\alias{note}
\alias{product}
\alias{proteinID}
\alias{qualif}
\alias{qualif<-}
\alias{translation}
\title{Get/set qualifiers of GenBank features}
\usage{
qualif(x, which, ...)

qualif(x, which, ...) <- value

locusTag(x)

product(x)

note(x)

proteinID(x)

translation(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}} or
  \code{\linkS4class{gbFeatureList}} instance.}

  \item{which}{(Optional) A character vector giving the
  name(s) of the qualifiers to retrieve.}

  \item{...}{Additional arguments passed to methods.}
}
\description{
Get/set qualifiers of GenBank features
}

