\docType{methods}
\name{qualif}
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

translation(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}},
  \code{\linkS4class{gbFeatureList}}, or
  \code{\linkS4class{gbRecord}} object.}

  \item{which}{(Optional) A character vector giving the
  name(s) of the qualifiers to retrieve.}

  \item{...}{Additional arguments passed to methods.}
}
\description{
Get/set qualifiers of GenBank features
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
qualif(x[[1]], 'db_xref')
}

