\docType{methods}
\name{qualif}
\alias{geneID}
\alias{locusTag}
\alias{note}
\alias{product}
\alias{proteinID}
\alias{qualif}
\alias{qualif,gbFeature-method}
\alias{qualif,gbFeatureList-method}
\alias{qualif,gbRecord-method}
\alias{qualif<-}
\alias{translation}
\title{Get/set qualifiers of GenBank features}
\usage{
qualif(x, which, ...)

qualif(x, which, ...) <- value

\S4method{qualif}{gbFeature}(x, which, fixed = FALSE, use.names = TRUE)

\S4method{qualif}{gbFeatureList}(x, which = "", fixed = FALSE,
  use.names = TRUE)

\S4method{qualif}{gbRecord}(x, which = "", fixed = FALSE,
  use.names = TRUE)

locusTag(x)

product(x)

note(x)

proteinID(x)

geneID(x)

translation(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}},
  \code{\linkS4class{gbFeatureList}}, or
  \code{\linkS4class{gbRecord}} object.}

  \item{which}{(Optional) A character vector giving the
  name(s) of the qualifiers to retrieve.}

  \item{fixed}{If \code{TRUE}, \code{which} is matched
  against qualifiers as is, if \code{FALSE} it is treated
  as a regular expression.}

  \item{use.names}{If \code{TRUE}, return a
  \code{data.frame} using \code{which} as column names, if
  \code{FALSE} return, if possible, a character vector or a
  list.}

  \item{...}{Additional arguments passed to methods.}
}
\description{
Get/set qualifiers of GenBank features
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
qualif(x[[1]], 'db_xref')

## use shortcuts to common qualifiers
proteinID(x["CDS"])
}

