\docType{class}
\name{gbFeatureList-class}
\alias{gbFeatureList-class}
\title{Class \code{"gbFeatureList"}}
\description{
\dQuote{gbFeatureList} is an S4 class that provides a
container for \dQuote{\linkS4class{gbFeature}}s retrived
from GenBank flat files.
}
\section{Slots}{

\describe{
\item{\code{.seqinfo}}{A \code{\linkS4class{seqinfo}} object containing the
genome sequence as an \code{\linkS4class{XStringSet}} object and
sequence metadata as a \code{\linkS4class{gbHeader}} object.}

\item{\code{.Data}}{A list of \code{\linkS4class{gbFeature}} objects.}
}}

