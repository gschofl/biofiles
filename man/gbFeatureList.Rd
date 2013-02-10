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
  \item \linkS4class{vectorORfactor}
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
  \item \code{\link{qualif}}
  \item \code{\link{range}}
  \item \code{\link{revcomp}}
  \item \code{\link{select}}
  \item \code{\link{sequence}}
  \item \code{\link{show}}
  \item \code{\link{start}}
  \item \code{\link{start<-}}
  \item \code{\link{strand}}
  \item \code{\link{strand<-}}
  \item \code{\link{summary}}
  \item \code{\link{view}}
}

}

\keyword{classes}
\section{Slots}{
\describe{
  \item{.Dir}{The path to the database file containing the GenBank
record the feature list is part of.}

  \item{.ACCN}{Accession number of the GenBank record that the
feature list is part of.}

  \item{.DEF}{The definition line (brief description of the sequence)
of the GenBank record the feature list is part of.}

  \item{.Data}{A list of \dQuote{\code{gbFeature}} instances.}

}
}
