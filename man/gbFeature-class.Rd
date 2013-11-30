\docType{class}
\name{gbFeature-class}
\alias{gbFeature-class}
\title{gbFeature-class}
\description{
\dQuote{gbFeature} is an S4 class that provides a container
for GenBank feature tables.
}
\section{Slots}{

\describe{
\item{\code{.seqinfo}}{An \code{environment} containing the genome sequence as
an \code{\linkS4class{XStringSet}} object and sequence metadata
as a \code{\linkS4class{gbHeader}} object.}

\item{\code{.id}}{Identifier (index) of the feature in the
GenBank record the feature is part of.}

\item{\code{key}}{The feature key.}

\item{\code{location}}{A \code{\linkS4class{gbLocation}} object.}

\item{\code{qualifiers}}{A named character vector. Name attributes
correspond to GenBank qualifier tags.}
}}
\seealso{
\code{\linkS4class{gbFeatureList}}
}

