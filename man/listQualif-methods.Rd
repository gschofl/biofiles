\docType{methods}
\name{listQualif}
\alias{listQualif}
\alias{listQualif,gbFeature-method}
\alias{listQualif,gbFeatureList-method}
\alias{listQualif,gbRecord-method}
\title{Access the names of Genbank qualifiers}
\usage{
listQualif(x, ...)

\S4method{listQualif}{gbFeature}(x)

\S4method{listQualif}{gbFeatureList}(x)

\S4method{listQualif}{gbRecord}(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}},
  \code{\linkS4class{gbFeatureList}}, or
  \code{\linkS4class{gbRecord}} object.}

  \item{...}{Additional arguments to be passed to or from
  methods.}
}
\value{
A character vector (or list of character vectors) of
qualifiers.
}
\description{
Access the names of Genbank qualifiers
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
listQualif(x["source"])
}
\seealso{
\code{\link{listUniqueQualifs}}, \code{\link{hasQualif}}
}

