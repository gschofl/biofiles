% Built by roxygen3 0.1
\name{upstream}
\alias{downstream}
\alias{flanking}
\alias{upstream}
\title{Find flanking features}
\usage{
  upstream(query, subject, n=5, include_key="all",
    exclude_key="none")

  downstream(query, subject, n=5, include_key="all",
    exclude_key="none)

  flanking(query, subject, n=5, include_key="all",
    exclude_key="none)
}
\arguments{
  \item{query}{The query \code{\linkS4class{gbRange}},
\code{\linkS4class{gbFeatureList}}, or
\code{\linkS4class{gbFeature}}  instance.}

  \item{subject}{The subject \code{\linkS4class{gbRecord}} or
\code{\linkS4class{gbFeatureList}} instance within which the
n nearest upstream features are searched.}

  \item{n}{The number of upstream features to be included}

  \item{include_key}{Which features should be returned. Defaults to
\dQuote{all}.}

  \item{exclude_key}{Which feature(s) should be excluded from the search.
Defaults to \dQuote{none}.}
}

\description{
Find flanking features
}
