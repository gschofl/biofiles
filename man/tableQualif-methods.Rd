\docType{methods}
\name{tableQualif}
\alias{tableQualif}
\alias{tableQualif,gbFeature-method}
\alias{tableQualif,gbFeatureList-method}
\alias{tableQualif,gbRecord-method}
\alias{tableQualif,gbRecordList-method}
\title{Tabulate the names of Genbank qualifiers}
\usage{
tableQualif(x, ...)

\S4method{tableQualif}{gbFeature}(x)

\S4method{tableQualif}{gbFeatureList}(x)

\S4method{tableQualif}{gbRecord}(x)

\S4method{tableQualif}{gbRecordList}(x)
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
A \code{\link{table}} (or list of \code{table}s) of
qualifiers names.
}
\description{
Extract a frequency table (or list of tables in the case of
\code{gbRecordList}s) of qualifier names.
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
tableQualif(x)
}
\seealso{
\code{\link{listUniqueQualifs}}, \code{\link{hasQualif}}
}

