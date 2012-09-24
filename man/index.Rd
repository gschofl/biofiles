% Built by roxygen3 0.1
\name{index}
\alias{index}
\title{Return feature indices from a GenBank record}
\usage{
  index(x, attributes = FALSE, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}} or
\code{\linkS4class{gbFeatureList}} instance.}

  \item{attributes}{Set the \code{accession}, \code{definition},
\code{database} attributes.}

  \item{...}{Additional arguments passed to methods.}
}

\value{
A numeric vector of feature indeces.
}
\description{
Return feature indices from a GenBank record
}
\section{Method signatures for generic index}{
\itemize{
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
}

}

\keyword{methods}
