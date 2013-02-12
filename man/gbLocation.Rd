% Built by roxygen3 0.1
\docType{class}
\name{gbLocation-class}
\alias{gbLocation-class}
\title{gbLocation}
\description{
\dQuote{gbLocation} is an S4 class that provides a container for
GenBank feature locations.
It extends \code{\linkS4class{Intervals_full}} and provides 5
additional Slots.
}
\details{
For more information see the 
\href{ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt}{GenBank Release Note}
}
\section{Class hierarchy for gbLocation}{
Super classes:
\itemize{
  \item \linkS4class{Intervals_full}
  \item \linkS4class{Intervals_virtual}
  \item \linkS4class{matrix}
  \item \linkS4class{Intervals_virtual_or_numeric}
  \item \linkS4class{array}
  \item \linkS4class{structure}
  \item \linkS4class{vector}
}

}

\section{Generics with methods for gbLocation}{
\itemize{
  \item \code{\link{accession}}
  \item \code{\link{end}}
  \item \code{\link{end<-}}
  \item \code{\link{initialize}}
  \item \code{\link{partial}}
  \item \code{\link{show}}
  \item \code{\link{start}}
  \item \code{\link{start<-}}
  \item \code{\link{strand}}
  \item \code{\link{strand<-}}
}

}

\keyword{classes}
\section{Slots}{
\describe{
  \item{strand}{A single integer -1 or 1 (for minus or plus strand,
  respectively)}

  \item{compound}{A character code specifying how multiple segments
  are joined. One of \sQuote{join} or \sQuote{order}.}

  \item{partial}{A logical matrix specifying whether residues are
  missing from the 5' and 3' ends respectively.}

  \item{accession}{An object of class \linkS4class{character}}

  \item{remote}{An object of class \linkS4class{logical}}

}
}
