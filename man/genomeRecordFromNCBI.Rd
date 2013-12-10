\name{genomeRecordFromNCBI}
\alias{genomeRecordFromNCBI}
\title{Retrieve genomes in GenBank format directly form NCBI}
\usage{
genomeRecordFromNCBI(which, ignore.case = TRUE, .parse = TRUE, ...)
}
\arguments{
  \item{which}{Path to ftp directory on NCBI's
  \href{ftp://ftp.ncbi.nih.gov/genomes/}{Genomes FTP site}.
  Examples would be \code{Bacteria/Acetinobacter*} or
  \code{Viruses/Bat_coronavirus*}. If there are multiple
  matching directories the user will be prompted to choose
  one. If there are multiple matching \code{gbk} files the
  user will also be prompted to choose one or more.}

  \item{ignore.case}{Ignore case when matching.}

  \item{...}{Arguments passed on to
  \code{\link[RCurl]{curlOptions}}.}
}
\value{
A \code{\linkS4class{gbRecord}} or
\code{\linkS4class{gbRecordList}} object.
}
\description{
## EXPERIMENTAL ##
}
\examples{
\dontrun{
gbk <- genomeRecordFromNCBI("Bacteria/Chlamydia_muridarum", verbose = TRUE)

}
}

