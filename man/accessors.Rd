\docType{methods}
\name{getLocus}
\alias{getAccession}
\alias{getComment}
\alias{getDBLink}
\alias{getDBSource}
\alias{getDate}
\alias{getDefinition}
\alias{getDivision}
\alias{getFeatures}
\alias{getGeneID}
\alias{getKeywords}
\alias{getLength}
\alias{getLocus}
\alias{getMoltype}
\alias{getOrganism}
\alias{getReference}
\alias{getSequence}
\alias{getSource}
\alias{getTaxonomy}
\alias{getTopology}
\alias{getVersion}
\title{Access slots of a GenBank record.}
\usage{
getLocus(x)

getLength(x)

getMoltype(x)

getTopology(x)

getDivision(x)

getDate(x)

getDefinition(x)

getAccession(x)

getVersion(x)

getGeneID(x)

getDBLink(x)

getDBSource(x)

getSource(x)

getOrganism(x)

getTaxonomy(x)

getKeywords(x)

getReference(x)

getComment(x)

getFeatures(x, ...)

getSequence(x, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbRecord}} instance.}

  \item{...}{Additional arguments passed to methods.}

  \item{x}{A \code{\linkS4class{gbRecord}},
  \code{\linkS4class{gbFeature}}, or
  \code{\linkS4class{gbFeatureList}} instance.}

  \item{...}{Additional arguments passed to methods.}
}
\value{
The \code{\linkS4class{gbFeatureList}} object contained in
a gbRecord database.

An \code{\linkS4class{XStringSet}} object.
}
\description{
Access slots of a GenBank record.

Retrieve feature tables from a GenBank record.

Retrieve sequence from a GenBank record or a Feature.
}

