\name{seqinfo-class}
\alias{seqinfo-class}
\title{Class \code{"seqinfo"}}
\description{
A container for shared date: Header and Sequence.
}
\section{Fields}{
  \describe{ \item{\code{header}:}{ A
  \code{\linkS4class{gbHeader}} object or \code{NULL}. }
  \item{\code{sequence}:}{ A
  \code{\linkS4class[Biostrings]{XStringSet}} object or
  \code{NULL}. } }
}

\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.
}
\examples{
showClass("seqinfo")
}
\seealso{
\code{\link{seqinfo}}
}
\keyword{classes}

