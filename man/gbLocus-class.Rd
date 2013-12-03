\name{gbLocus-class}
\alias{gbLocus-class}
\title{Class \code{"gbLocus"}}
\description{
A container for GenBank LOCUS records
}
\section{Fields}{
  \describe{ \item{\code{lnm}:}{ Locus name. Usually the
  accession number. } \item{\code{len}:}{ Sequence length;
  In bp or aa, depending on \code{mtp}. }
  \item{\code{mtp}:}{ Molecule type; NA, DNA, RNA, tRNA
  (transfer RNA), rRNA (ribosomal RNA), mRNA (messenger
  RNA), uRNA (small nuclear RNA), or AA (protein sequence).
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

