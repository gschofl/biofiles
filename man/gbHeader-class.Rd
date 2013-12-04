\name{gbHeader-class}
\alias{gbHeader-class}
\title{Class \code{"gbHeader"}}
\description{
A class representing a GenBank/GenPept-format file header.
}
\section{Fields}{
  \describe{ \item{\code{locus}:}{ A
  \code{"\linkS4class{gbLocus}"} object. }
  \item{\code{definition}:}{ \code{character}; Description
  of the sequence. } \item{\code{accession}:}{
  \code{character}; The primary accession number.  A
  unique, unchanging identifier assigned to each GenBank
  sequence record. } \item{\code{version}:}{
  \code{character}; The primary accession number and a
  numeric version number } \item{\code{seqid}:}{
  \code{character}; Gene Identifier ("GI"); An integer key
  assigned to the sequence by NCBI. }
  \item{\code{dblink}:}{ \code{character}; Cross-references
  to resources that support the existence a sequence
  record, such as the Project Database and the NCBI Trace
  Assembly Archive. } \item{\code{dbsource}:}{
  \code{character}; GenPept files only }
  \item{\code{keywords}:}{ code{character}; Short
  description of gene products and other information about
  an entry. } \item{\code{source}:}{ \code{character};
  Common name of the organism. } \item{\code{organism}:}{
  \code{character}; Formal scientific name of the organism.
  } \item{\code{taxonomy}:}{ \code{character}; Taxonomic
  classification levels. } \item{\code{references}:}{
  Citations for all articles containing data reported in
  the entry.A \code{\linkS4class{gbReferenceList}} object.
  } \item{\code{comment}:}{ \code{character}; Remarks. } }
}

\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.
}
\examples{
showClass("gbHeader")
}
\seealso{
\code{\link{.gbHeader}}
}
\keyword{classes}
\keyword{internal}

