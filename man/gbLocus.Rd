\name{.gbLocus}
\alias{.gbLocus}
\title{Generator object for the \code{\linkS4class{gbLocus}} class}
\usage{
.gbLocus(...)
}
\arguments{
  \item{...}{List of arguments (see NOTE)}
}
\description{
The generator object for the \code{\linkS4class{gbLocus}}
reference class.
}
\note{
Arguments to the \code{#new} method must be named
arguments: \itemize{ \item{lnm}{ Locus name; stored in the
\code{lnm} field. } \item{len}{ Sequence lenght; stored in
the \code{len} field. } \item{mtp}{ Molecule type; stored
in the \code{mtp} field). } \item{div}{ Genbank division;
stored in the \code{div} field. } \item{top}{ Topology;
stored in the \code{top} field. } \item{mdt}{ Modification
date; stored in the \code{mdt} field. } \item{cdt}{ Create
date; stored in the \code{cdt} field. } }
}
\section{Methods}{
  \describe{ \item{\code{#new(lnm, len, mtp, div, top, mdt,
  cdt)}:}{ Create a new \code{\linkS4class{gbLocus}}
  object} }
}
\seealso{
\code{\linkS4class{gbLocus}}
}
\keyword{classes}
\keyword{internal}

