% Built by roxygen3 0.1
\name{location}
\alias{location}
\title{Return feature locations from a GenBank record.}
\usage{
  location(x, attributes = FALSE, join = FALSE, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}} or
\code{\linkS4class{gbFeatureList}} instance.}

  \item{attributes}{Set the \code{accession}, \code{definition},
\code{database} attributes.}

  \item{join}{Combine compound locations}

  \item{...}{Additional arguments passed to methods.}
}

\value{
A \code{\linkS4class{gbRange}} object.
}
\description{
Return feature locations from a GenBank record.
}
\section{Method signatures for generic location}{
\itemize{
  \item \linkS4class{gbFeature}
  \item \linkS4class{gbFeatureList}
}

}

\keyword{methods}
