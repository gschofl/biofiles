% Built by roxygen3 0.1
\docType{class}
\name{gbFeature-class}
\alias{gbFeature}
\alias{gbFeature-class}
\title{gbFeature}
\usage{
  gbFeature(feature, seqinfo, accession = "", id = 0L)
}
\description{
\dQuote{gbFeature} is an S4 class that provides a container
for GenBank feature tables.
}
\section{Generics with methods for gbFeature}{
\itemize{
  \item \code{\link{.dbSource}}
  \item \code{\link{.defline}}
  \item \code{\link{.sequence}}
  \item \code{\link{[[}}
  \item \code{\link{$}}
  \item \code{\link{dbxref}}
  \item \code{\link{end}}
  \item \code{\link{end<-}}
  \item \code{\link{fuzzy}}
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
  \item \code{\link{key<-}}
  \item \code{\link{listQualif}}
  \item \code{\link{qualif}}
  \item \code{\link{qualif<-}}
  \item \code{\link{ranges}}
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
as a \code{\linkS4class{gbHeader}} object.}

  \item{.id}{Identifier (index) of the feature in the
GenBank record the feature is part of.}

  \item{key}{The feature key.}

  \item{location}{A \code{\linkS4class{gbLocation}} object.}

  \item{qualifiers}{A named character vector. Name attributes
correspond to GenBank qualifier tags.}

}
}
