% Built by roxygen3 0.1
\name{revcomp}
\alias{revcomp}
\title{reverse complement features in a GenBank record}
\usage{
  revcomp(x, order = FALSE, updateDb = FALSE, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbRecord}} or
\code{\linkS4class{gbFeatureList}} or instance.
(\code{gbFeatureList} instances must include a \sQuote{source} field).}

  \item{order}{Should the resulting gbFeatureList be reordered.}

  \item{update_db}{Should the sequence and the new feature locations be
updated in the underlying filehash database?}

  \item{...}{Additional arguments passed to methods.}
}

\value{
A \code{\linkS4class{gbFeatureList}} instance
}
\description{
reverse complement features in a GenBank record
}
\section{Method signatures for generic revcomp}{
\itemize{
  \item \linkS4class{gbFeatureList}
  \item \linkS4class{gbRecord}
}

}

\keyword{methods}
