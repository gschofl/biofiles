% Built by roxygen3 0.1
\docType{class}
\name{gbFeatureList-class}
\alias{gbFeatureList-class}
\title{gbFeatureList}
\description{
\dQuote{gbFeatureList} is an S4 class that provides a container for
\dQuote{\linkS4class{gbFeature}}s retrived from GenBank flat files.
}
\section{Class hierarchy for gbFeatureList}{
Super classes:
\itemize{
  \item \linkS4class{list}
  \item \linkS4class{vector}
}

}

\section{Generics with methods for gbFeatureList}{
\itemize{
  \item \code{\link{accession}}
  \item \code{\link{dbxref}}
  \item \code{\link{definition}}
  \item \code{\link{end}}
  \item \code{\link{end<-}}
  \item \code{\link{hasKey}}
  \item \code{\link{hasQualif}}
  \item \code{\link{index}}
  \item \code{\link{key}}
  \item \code{\link{listQualif}}
  \item \code{\link{qualif}}
  \item \code{\link{ranges}}
  \item \code{\link{select}}
  \item \code{\link{seqinfo}}
  \item \code{\link{seqlengths}}
  \item \code{\link{sequence}}
  \item \code{\link{show}}
  \item \code{\link{start}}
  \item \code{\link{start<-}}
  \item \code{\link{strand}}
  \item \code{\link{strand<-}}
  \item \code{\link{summary}}
  \item \code{\link{view}}
  \item \code{\link{width}}
}

}

\keyword{classes}
\section{Slots}{
\describe{
  \item{.seqinfo}{An \code{environment} containing the genome sequence as
an \code{\linkS4class{XStringSet}} object and sequence metadata
as a \code{\linkS4class{Seqinfo}} object.}

  \item{.Data}{A list of \code{\linkS4class{gbFeature}} objects.}

}
}
