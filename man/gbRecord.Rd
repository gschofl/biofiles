\name{gbRecord}
\alias{gbRecord}
\title{Read a GenBank/GenPept-format file.}
\usage{
gbRecord(gbk, with_sequence = TRUE)
}
\arguments{
  \item{gbk}{A vector of paths to GenBank format files, an
  \code{\link[reutils]{efetch}} object containing GenBank
  record(s), or a \code{textConnection} to a character
  vector that can be parsed as a Genbank record.}

  \item{with_sequence}{Include sequence information if
  avaliable.}
}
\value{
An instance of the \code{\linkS4class{gbRecord}} or
\code{\linkS4class{gbRecordList}} classes.
}
\description{
Import data from GenBank/GenPept files into R, represented
as an instance of the \code{\linkS4class{gbRecord}} or
\code{\linkS4class{gbRecordList}} classes.
}
\details{
For a sample GenBank record see
\url{http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html},
for a detailed description of the GenBank feature table
format see \url{http://www.ncbi.nlm.nih.gov/collab/FT/}
}
\examples{
### import from file
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)
getHeader(x)
getFeatures(x)

### quickly extract features as GRanges
ranges(x["CDS"], include=c("product", "note", "protein_id"))

### import directly from NCBI
\dontrun{
require(reutils)
x <- gbRecord(efetch("139189709", "protein", rettype="gp", retmode="text"))
x
}

## Directly subset features
x[[1]]
}

