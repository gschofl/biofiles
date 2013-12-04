\docType{methods}
\name{hasQualif}
\alias{hasQualif}
\alias{hasQualif,gbFeature-method}
\alias{hasQualif,gbFeatureList-method}
\alias{hasQualif,gbRecord-method}
\alias{hasQualif,gbRecordList-method}
\title{Has a feature a specific qualifier?}
\usage{
hasQualif(x, qualifier, ...)

\S4method{hasQualif}{gbFeature}(x, qualifier)

\S4method{hasQualif}{gbFeatureList}(x, qualifier)

\S4method{hasQualif}{gbRecord}(x, qualifier)

\S4method{hasQualif}{gbRecordList}(x, qualifier)
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
Has a feature a specific qualifier?
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
hasQualif(x, 'CDS')
}
\seealso{
\code{\link{listQualif}}, to extract a list of available
qualifiers for each feature;
\code{\link{listUniqueQualif}}, for a vector of all unique
qualifiers present in an object.
}

