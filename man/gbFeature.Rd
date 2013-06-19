% Built by roxygen3 0.1
\docType{class}
\name{gbFeature-class}
\alias{gbFeature}
\alias{gbFeature-class}
\title{gbFeature}
\usage{
  gbFeature(feature, seqenv, id = 0L)
}
\description{
\dQuote{gbFeature} is an S4 class that provides a container
for GenBank feature tables.
}
\section{Generics with methods for gbFeature}{
\itemize{
  \item \code{\link{[[}}
  \item \code{\link{$}}
  \item \code{\link{accession}}
  \item \code{\link{dbxref}}
  \item \code{\link{definition}}
  \item \code{\link{end}}
  \item \code{\link{end<-}}
  \item \code{\link{fuzzy}}
  \item \code{\link{hasKey}}
  \item \code{\link{hasQualif}}
  \item \code{\link{index}}
  \item \code{\link{key}}
  \item \code{\link{key<-}}
  \item \code{\link{listQualif}}
  \item \code{\link{qualif}}
  \item \code{\link{qualif<-}}
  \item \code{\link{ranges}}
  \item \code{\link{seqinfo}}
  \item \code{\link{seqlengths}}
  \item \code{\link{sequence}}
  \item \code{\link{shift}}
  \item \code{\link{show}}
  \item \code{\link{start}}
  \item \code{\link{start<-}}
  \item \code{\link{strand}}
  \item \code{\link{strand<-}}
  \item \code{\link{summary}}
  \item \code{\link{width}}
}

}

\keyword{classes}
\section{Slots}{
\describe{
  \item{.seqinfo}{An \code{environment} containing the genome sequence as
an \code{\linkS4class{XStringSet}} object and sequence metadata
as a \code{\linkS4class{Seqinfo}} object.}

  \item{.id}{Identifier (index) of the feature in the
GenBank record the feature is part of.}

  \item{key}{The feature key.}

  \item{location}{A \code{\linkS4class{gbLocation}} object.}

  \item{qualifiers}{A named character vector. Name attributes
correspond to GenBank qualifier tags.}

}
}
