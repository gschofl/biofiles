\docType{methods}
\name{end}
\alias{end}
\alias{end<-}
\title{Get or set the end of genomic features}
\usage{
end(x, join = FALSE, ...)

end(x, check = TRUE) <- value
}
\arguments{
  \item{x}{A \code{gbFeature} or \code{gbFeatureList}
  object.}

  \item{join}{Join compound genomic locations into a single
  range.}

  \item{...}{Further arguments passed to methods.}
}
\value{
An integer vector
}
\description{
Get or set the end of genomic features
}

