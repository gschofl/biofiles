% Built by roxygen3 0.1
\name{write.FeatureTable}
\alias{write.FeatureTable}
\title{Write GenBank records or features to file in Feature Table format}
\usage{
  write.FeatureTable(x, file, tablename = "", dbname = "",
    sequence = FALSE, append = FALSE, ...)
}
\arguments{
  \item{x}{A \code{\linkS4class{gbRecord}} instance.}

  \item{file}{A connection or a character string naming the file to write to.}

  \item{tablename}{Optional table name to appear in the first line
of the feature table.}

  \item{dbname}{Data base name associated with the CDS qualifier protein_id.}

  \item{sequence}{if \code{TRUE}, additionally autput fasta file}

  \item{append}{if \code{TRUE} the data is appended to the connection.}
}

\description{
Feature Tables are simple five-column tab-delimited tables specifying the
location and type of each feature. They can be used as input for tbl2asn
or Sequin to generate annotation.
}
\section{Method signatures for generic write.FeatureTable}{
\itemize{
  \item \linkS4class{gbRecord}
}

}

\keyword{methods}
