% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbHeader-class.R
\docType{class}
\name{gbReferenceList-class}
\alias{gbReferenceList-class}
\alias{.gbReferenceList}
\title{Generator object for the \code{"gbReferenceList"} reference class}
\arguments{
\item{...}{List of arguments}
}
\description{
Generator object for the \code{"gbReferenceList"} reference class

Generates an object representing a set of GenBank REFERENCE fields.
}
\section{Methods}{

\describe{
\item{\code{to_string(write_to_file = FALSE)}}{Generate a character string representation of a GenBank reference list}
}}

\section{Methods}{

\describe{
\item{\code{#new(ref)}:}{
   Create a new \code{\linkS4class{gbReferenceList}} object }
\item{\code{#to_string(write_to_file = FALSE)}:}{
   Create a string representation of a GenBank REFERENCE list. }
}
}

\section{Fields}{

\describe{
\item{\code{ref}:}{ A list of \code{"\linkS4class{gbReference}"} objects. }
}
}

\section{Extends}{
 All reference classes extend and inherit methods from
   \code{"\linkS4class{envRefClass}"}.
}

\examples{
showClass("gbReferenceList")
}
\keyword{classes}
\keyword{internal}
