\name{gbLocus-class}
\alias{gbLocus-class}
\title{Class \code{"gbLocus"}}
\description{
A class representing a GenBank LOCUS line.
}
\section{Fields}{
  \describe{ \item{\code{lnm}:}{ Locus name. Usually the
  accession number. } \item{\code{len}:}{ Sequence length;
  In bp or aa, depending on \code{mtp}. }
  \item{\code{mtp}:}{ Molecule type; \emph{NA}, \emph{DNA},
  \emph{RNA}, \emph{tRNA} (transfer RNA), \emph{rRNA}
  (ribosomal RNA), \emph{mRNA} (messenger RNA), \emph{uRNA}
  (small nuclear RNA), or \emph{AA} (protein sequence).
  RNAs can be prefixes ss- (single-stranded), ds-
  (double-stranded), or ms- (mixed-stranded)}
  \item{\code{div}:}{ Genbank division. }
  \item{\code{top}:}{ Topology; linear, circular, or
  missing (\code{NA}). } \item{\code{mdt}:}{ Modification
  date. } \item{\code{cdt}:}{ Create date. } }
}

\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.
}
\examples{
showClass("gbLocus")
}
\seealso{
\code{\link{.gbLocus}}
}
\keyword{classes}
\keyword{internal}

