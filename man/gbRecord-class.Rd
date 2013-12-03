\name{gbRecord-class}
\alias{gbRecord-class}
\title{Class \code{"gbRecord"}}
\usage{
new_gbRecord(...)
}
\description{
\dQuote{gbRecord} is an S4 class that provides a container
for data parsed from a GenBank or GenPept records. For
instantiation a gbRecord object use the import function
\code{\link{gbRecord}}.
}
\section{Slots}{

\describe{
\item{\code{seqinfo}}{A \code{\linkS4class{seqinfo}} instance; This is a
reference class holding the sequence as an \code{\linkS4class{XStringSet}}
instance and header of the file containing metadata as a
\code{\linkS4class{gbHeader}} object.}

\item{\code{features}}{A \code{\linkS4class{gbFeatureList}} instance.}

\item{\code{contig}}{If present, a CONTIG record.}
}}

