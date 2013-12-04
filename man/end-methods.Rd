\docType{methods}
\name{end}
\alias{end}
\alias{end,gbFeature-method}
\alias{end,gbFeatureList-method}
\alias{end,gbRecord-method}
\alias{end,gbRecordList-method}
\alias{end<-}
\title{Get or set the end of genomic features}
\usage{
end(x, join = FALSE, ...)

end(x, check = TRUE) <- value

\S4method{end}{gbFeature}(x, join = FALSE, drop = TRUE)

\S4method{end}{gbFeatureList}(x, join = FALSE, drop = TRUE)

\S4method{end}{gbRecord}(x, join = FALSE, drop = TRUE)

\S4method{end}{gbRecordList}(x, join = FALSE, drop = TRUE)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbFeature}},
  \code{\linkS4class{gbFeatureList}},
  \code{\linkS4class{gbRecord}}, or
  \code{\linkS4class{gbRecordList}} object.}

  \item{join}{Join compound genomic locations into a single
  range.}

  \item{...}{Further arguments passed to methods.}
}
\value{
An integer vector or a list of integer vectors.
}
\description{
Get or set the end of genomic features
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
end(x)

cds <- x["CDS"]
end(cds)
}
\seealso{
\code{\link{start}}, \code{\link{strand}},
\code{\link{width}}, \code{\link{ranges}}
}

