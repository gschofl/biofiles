% Built by roxygen3 0.1
\name{write_feature_table}
\alias{write_feature_table}
\title{General function for writing out NCBI feature tables}
\usage{
  write_feature_table(db, tablename = "", dbname = "",
    with_sequence = TRUE, outfile = "out.tbl")
}
\arguments{
  \item{db}{A \code{\linkS4class{gbRecord}} instance.}

  \item{tablename}{(Optional) Optional table name to appear in the first line
of the feature table.}

  \item{dbname}{Data base name associated with the CDS qualifier protein_id.}

  \item{with_sequence}{Additionally autput fasta file}

  \item{outfile}{Output file.}
}

\value{
NULL
}
\description{
Feature tables are simple five-column tab-delimited tables specifying the
location and type of each feature. They can be used as input for tbl2asn
or Sequin to generate annotation.
}
