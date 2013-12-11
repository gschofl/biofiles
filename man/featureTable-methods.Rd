\docType{methods}
\name{featureTable}
\alias{featureTable}
\alias{featureTable,gbFeatureList-method}
\alias{featureTable,gbRecord-method}
\alias{featureTable,gbRecordList-method}
\title{Tabulate Genbank features}
\usage{
featureTable(x, ...)

\S4method{featureTable}{gbFeatureList}(x)

\S4method{featureTable}{gbRecord}(x)

\S4method{featureTable}{gbRecordList}(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeatureList}},
  \code{\linkS4class{gbRecord}}, or
  \code{\linkS4class{gbRecordList}} object.}

  \item{...}{Additional arguments to be passed to or from
  methods.}
}
\value{
A \code{\link{table}} (or list of \code{table}s) of feature
keys.
}
\description{
Extract a frequency table (or list of tables in the case of
\code{gbRecordList}s) of feature keys.
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
featureTable(x)
}
\seealso{
\code{\link{qualifTable}}
}

