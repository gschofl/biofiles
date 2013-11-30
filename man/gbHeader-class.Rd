\docType{class}
\name{gbHeader-class}
\alias{gbHeader-class}
\title{gbHeader-class}
\description{
Metadata about a sequence record.
}
\section{Slots}{

\describe{
\item{\code{locus}}{\code{gbLocus}; locus name, sequence length, molecule type,
     topology (optional), division code, modification date.}

\item{\code{definition}}{\code{character}; Description of the sequence.}

\item{\code{accession}}{\code{character}; The primary accession number. A unique,
     unchanging identifier assigned to each GenBank sequence record.}

\item{\code{version}}{\code{character}; The primary accession number and a numeric
     version number}

\item{\code{seqid}}{\code{character}; Gene Identifier ("GI"); An integer key assigned
     to the sequence by NCBI.}

\item{\code{dblink}}{\code{character}; Cross-references to resources that
     support the existence a sequence record, such as the Project Database
     and the NCBI Trace Assembly Archive.}

\item{\code{dbsource}}{}

\item{\code{keywords}}{\code{character}; Short description of gene products and other
     information about an entry.}

\item{\code{source}}{\code{character}; Common name of the organism.}

\item{\code{organism}}{\code{character}; Formal scientific name of the organism.}

\item{\code{taxonomy}}{\code{character}; Taxonomic classification levels.}

\item{\code{references}}{\code{character}; Citations for all articles containing data
     reported in the entry.}

\item{\code{comment}}{\code{character}; Remarks.}
}}

