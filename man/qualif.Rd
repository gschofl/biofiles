% Built by roxygen3 0.1
\name{qualif}
\alias{qualif}
\alias{qualif<-}
\title{Get/set feature qualifiers from a GenBank record}
\usage{
  qualif(x, which, attributes = FALSE, ...)

  qualif(x, which, value, ...) <- value
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}} or
\code{\linkS4class{gbFeatureList}} instance.}

  \item{which}{(Optional) A character vector giving the name(s) of the
qualifiers to retrieve.}

  \item{attributes}{Set the \code{accession}, \code{definition},
\code{database} attributes.}

  \item{...}{Additional arguments passed to methods.}
}

\description{
Get/set feature qualifiers from a GenBank record
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
