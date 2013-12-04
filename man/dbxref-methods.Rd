\docType{methods}
\name{dbxref}
\alias{dbxref}
\alias{dbxref,gbFeature-method}
\alias{dbxref,gbFeatureList-method}
\alias{dbxref,gbRecord-method}
\title{Access the \code{db_xref}s of GenBank features}
\usage{
dbxref(x, db = NULL, ...)

\S4method{dbxref}{gbFeature}(x, db = NULL, ...)

\S4method{dbxref}{gbFeatureList}(x, db = NULL, na.rm = TRUE, ...)

\S4method{dbxref}{gbRecord}(x, db = NULL, na.rm = TRUE, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}},
  \code{\linkS4class{gbFeatureList}}, or
  \code{\linkS4class{gbRecord}} object.}

  \item{db}{(Optional) A character vector giving the
  database names of the desired \code{db_xref}s.}

  \item{...}{Additional arguments passed to methods.}
}
\value{
A named character vector (or list of named character
vectors) of db_xrefs.
}
\description{
Access the \code{db_xref}s of GenBank features
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)

## all db_xrefs associated with CDSs
dbxref(x["CDS"])

## retrieve the TaxId from the "source" field.
dbxref(x[[1]], "taxon")
}

