\docType{methods}
\name{width}
\alias{joint_width}
\alias{joint_width,gbFeature-method}
\alias{joint_width,gbFeatureList-method}
\alias{joint_width,gbRecord-method}
\alias{width}
\alias{width,gbFeature-method}
\alias{width,gbFeatureList-method}
\alias{width,gbRecord-method}
\alias{width,gbRecordList-method}
\title{Get the width of genomic features}
\usage{
width(x)

joint_width(x)

\S4method{width}{gbFeature}(x)

\S4method{joint_width}{gbFeature}(x)

\S4method{width}{gbFeatureList}(x)

\S4method{joint_width}{gbFeatureList}(x)

\S4method{width}{gbRecord}(x)

\S4method{joint_width}{gbRecord}(x)

\S4method{width}{gbRecordList}(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}},
  \code{\linkS4class{gbFeatureList}},
  \code{\linkS4class{gbRecord}}, or
  \code{\linkS4class{gbRecordList}} object.}
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

