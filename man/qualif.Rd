% Built by roxygen3 0.1
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

  qualif(x, which, value, ...) <- value

  locusTag(x)

  product(x)

  note(x)

  proteinID(x)

  translation(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}} or
\code{\linkS4class{gbFeatureList}} instance.}

  \item{which}{(Optional) A character vector giving the name(s) of the
qualifiers to retrieve.}

  \item{...}{Additional arguments passed to methods.}
}

\description{
Get/set qualifiers of GenBank features
}
\section{Method signatures for generic qualif}{
\itemize{
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
}

}

\section{Method signatures for generic qualif<-}{
\itemize{
  \item \linkS4class{gbFeature}
}

}

\keyword{methods}
