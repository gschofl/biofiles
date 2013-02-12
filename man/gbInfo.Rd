% Built by roxygen3 0.1
\docType{class}
\name{gbInfo-class}
\alias{gbInfo}
\alias{gbInfo-class}
\title{gbInfo}
\usage{
  gbInfo(object)
}
\description{
\dQuote{gbInfo} is an S4 class that provides basic infomation
about a genomic sequence. It extends the \linkS4class{Seqinfo}.
}
\section{Class hierarchy for gbInfo}{
Super classes:
\itemize{
  \item \linkS4class{Seqinfo}
}

}

\keyword{classes}
\section{Slots}{
\describe{
  \item{db}{A \code{\linkS4class{filehashRDS}} instance.}

  \item{seqnames}{Identifier of the sequence. Usually the Accession number.}

  \item{seqlengths}{Length of the sequence.}

  \item{is_circular}{Is the sequence circular.}

  \item{genome}{The definition of the sequence data.}

}
}
