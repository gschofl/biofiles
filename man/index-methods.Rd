\docType{methods}
\name{index}
\alias{index}
\alias{index,gbFeature-method}
\alias{index,gbFeatureList-method}
\alias{index,gbRecord-method}
\alias{index,gbRecordList-method}
\title{Access the indices of GenBank features}
\usage{
index(x, ...)

\S4method{index}{gbFeature}(x)

\S4method{index}{gbFeatureList}(x)

\S4method{index}{gbRecord}(x)

\S4method{index}{gbRecordList}(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}},
  \code{\linkS4class{gbFeatureList}},
  \code{\linkS4class{gbRecord}}, or
  \code{\linkS4class{gbRecordList}} object.}

  \item{...}{Additional arguments passed to methods.}
}
\value{
A numeric vector of feature indeces.
}
\description{
Access the indices of GenBank features
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
index(x)
}

