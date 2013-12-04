\name{gbReference-class}
\alias{gbReference-class}
\title{Class \code{"gbReference"}}
\description{
A class representing a GenBank REFERENCE field.
}
\section{Fields}{
  \describe{ \item{\code{refline}:}{ Top line of a
  reference entry.} \item{\code{authors}:}{ Authors of the
  citation. Mandatory or \code{consrtm}. }
  \item{\code{consrtm}:}{ The collective names of
  consortiums. Optional.} \item{\code{title}:}{ Full title
  of citation. Optional. } \item{\code{journal}:}{ The
  journal name, volume, year, and page numbers of the
  citation. Mandatory. } \item{\code{pubmed}:}{ The PubMed
  unique identifier for a citation. Optional. }
  \item{\code{remark}:}{ The relevance of a citation.
  Optional. } }
}

\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.
}
\examples{
showClass("gbReference")
}
\seealso{
\code{\link{.gbReference}}
}
\keyword{classes}
\keyword{internal}

