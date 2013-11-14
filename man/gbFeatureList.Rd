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
  \item \code{\link{.dbSource}}
  \item \code{\link{.defline}}
  \item \code{\link{.sequence}}
  \item \code{\link{dbxref}}
  \item \code{\link{end}}
  \item \code{\link{end<-}}
  \item \code{\link{getAccession}}
  \item \code{\link{getComment}}
  \item \code{\link{getDate}}
  \item \code{\link{getDBLink}}
  \item \code{\link{getDBSource}}
  \item \code{\link{getDefinition}}
  \item \code{\link{getDivision}}
  \item \code{\link{getGeneID}}
  \item \code{\link{getKeywords}}
  \item \code{\link{getLength}}
  \item \code{\link{getLocus}}
  \item \code{\link{getMoltype}}
  \item \code{\link{getOrganism}}
  \item \code{\link{getReference}}
  \item \code{\link{getSequence}}
  \item \code{\link{getSource}}
  \item \code{\link{getTaxonomy}}
  \item \code{\link{getTopology}}
  \item \code{\link{getVersion}}
  \item \code{\link{hasKey}}
  \item \code{\link{hasQualif}}
  \item \code{\link{index}}
  \item \code{\link{key}}
  \item \code{\link{listQualif}}
  \item \code{\link{qualif}}
  \item \code{\link{ranges}}
  \item \code{\link{revcomp}}
  \item \code{\link{select}}
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
  \item{.seqinfo}{An \code{seqinfo} containing the genome sequence as
an \code{\linkS4class{XStringSet}} object and sequence metadata
as a \code{\linkS4class{gbHeader}} object.}

  \item{.Data}{A list of \code{\linkS4class{gbFeature}} objects.}

}
}
