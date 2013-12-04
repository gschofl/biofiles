\docType{methods}
\name{key}
\alias{key}
\alias{key,gbFeature-method}
\alias{key,gbFeatureList-method}
\alias{key,gbRecord-method}
\alias{key,gbRecordList-method}
\alias{key<-}
\title{Get/set keys of GenBank features}
\usage{
key(x, ...)

key(x, ...) <- value

\S4method{key}{gbFeature}(x)

\S4method{key}{gbFeatureList}(x)

\S4method{key}{gbRecord}(x)

\S4method{key}{gbRecordList}(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}},
  \code{\linkS4class{gbFeatureList}},
  \code{\linkS4class{gbRecord}}, or
  \code{\linkS4class{gbRecordList}} object.}

  \item{...}{Additional arguments passed to methods.}
}
\description{
Get/set keys of GenBank features
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
key(x)
}

