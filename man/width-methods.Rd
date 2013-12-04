\docType{methods}
\name{width}
\alias{width}
\alias{width,gbFeature-method}
\alias{width,gbFeatureList-method}
\alias{width,gbRecord-method}
\alias{width,gbRecordList-method}
\title{Get the width of genomic features}
\usage{
width(x, ...)

\S4method{width}{gbFeature}(x, join = FALSE)

\S4method{width}{gbFeatureList}(x, join = FALSE)

\S4method{width}{gbRecord}(x, join = FALSE)

\S4method{width}{gbRecordList}(x, join = FALSE)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}},
  \code{\linkS4class{gbFeatureList}},
  \code{\linkS4class{gbRecord}}, or
  \code{\linkS4class{gbRecordList}} object.}

  \item{...}{Further arguments passed to methods.}
}
\value{
An integer vector or a list of integer vectors.
}
\description{
Get the width of genomic features
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
width(x)
}
\seealso{
\code{\link{start}}, \code{\link{end}},
\code{\link{strand}}, \code{\link{ranges}}
}

