\docType{methods}
\name{hasKey}
\alias{hasKey}
\alias{hasKey,gbFeature-method}
\alias{hasKey,gbFeatureList-method}
\alias{hasKey,gbRecord-method}
\alias{hasKey,gbRecordList-method}
\title{Has a feature a specific key?}
\usage{
hasKey(x, key, ...)

\S4method{hasKey}{gbFeature}(x, key)

\S4method{hasKey}{gbFeatureList}(x, key)

\S4method{hasKey}{gbRecord}(x, key)

\S4method{hasKey}{gbRecordList}(x, key)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}},
  \code{\linkS4class{gbFeatureList}},
  \code{\linkS4class{gbRecord}}, or
  \code{\linkS4class{gbRecordList}} object.}

  \item{...}{Additional arguments to be passed to or from
  methods.}
}
\value{
A logical vector or a list of logical vectors.
}
\description{
Has a feature a specific key?
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
hasKey(x, 'CDS')
}
\seealso{
\code{\link{key}}
}

