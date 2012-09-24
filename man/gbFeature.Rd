% Built by roxygen3 0.1
\docType{class}
\name{gbFeature-class}
\alias{gbFeature-class}
\title{gbFeature class}
\description{
\sQuote{\code{gbFeature}} is an S4 class that extends the class
\code{\linkS4class{gbFeatureList}}. This class provides a container
for feature data retrived from GenBank flat files.
}
\section{Generics with methods for gbFeature}{
\itemize{
  \item \code{\link{[[}}
  \item \code{\link{$}}
  \item \code{\link{accession}}
  \item \code{\link{dbxref}}
  \item \code{\link{definition}}
  \item \code{\link{end}}
  \item \code{\link{end<-}}
  \item \code{\link{hasKey}}
  \item \code{\link{hasQualif}}
  \item \code{\link{index}}
  \item \code{\link{key}}
  \item \code{\link{key<-}}
  \item \code{\link{partial}}
  \item \code{\link{qualif}}
  \item \code{\link{qualif<-}}
  \item \code{\link{range}}
  \item \code{\link{sequence}}
  \item \code{\link{show}}
  \item \code{\link{start}}
  \item \code{\link{start<-}}
  \item \code{\link{summary}}
}

}

\keyword{classes}
\section{Slots}{
\describe{
  \item{.Dir}{The path to the database file containing the GenBank
record the feature is part of.}

  \item{.ACCN}{Accession number of the GenBank record that the
feature is part of.}

  \item{.DEF}{The definition line (brief description of the sequence)
of the GenBank record the feature is part of.}

  \item{.ID}{Identifier (sequential index) of the feature in the
GenBank record the feature is part of.}

  \item{key}{An object of class \linkS4class{character}}

  \item{location}{An object of class \linkS4class{gbLocation}}

  \item{qualifiers}{Named character vector. Name attributes
   correspond to GenBank qualifier tags.}

}
}
