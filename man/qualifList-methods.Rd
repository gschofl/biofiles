\docType{methods}
\name{qualifList}
\alias{qualifList}
\alias{qualifList,gbFeature-method}
\alias{qualifList,gbFeatureList-method}
\alias{qualifList,gbRecord-method}
\title{List the names of Genbank qualifiers.}
\usage{
qualifList(x, ...)

\S4method{qualifList}{gbFeature}(x)

\S4method{qualifList}{gbFeatureList}(x)

\S4method{qualifList}{gbRecord}(x)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}},
  \code{\linkS4class{gbFeatureList}},
  \code{\linkS4class{gbRecord}}, or
  \code{\linkS4class{gbRecord}} object.}

  \item{...}{Additional arguments to be passed to or from
  methods.}
}
\value{
A character vector (or list of character vectors) of
qualifier names.
}
\description{
List the names of Genbank qualifiers.
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
qualifList(x["source"])
}
\seealso{
\code{\link{uniqueQualifs}}, \code{\link{hasQualif}}
}

