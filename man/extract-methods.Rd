\docType{methods}
\name{[,gbRecord,ANY-method}
\alias{[,gbRecord,ANY,ANY,ANY-method}
\alias{[,gbRecord,ANY-method}
\alias{[,gbRecord,missing,ANY,ANY-method}
\alias{[,gbRecord,missing-method}
\alias{[[,gbRecord,ANY,ANY-method}
\alias{[[,gbRecord-method}
\title{Method extensions to extraction operator for gbRecord objects.}
\usage{
\S4method{[}{gbRecord,ANY}(x, i, j, ..., drop = TRUE)

\S4method{[}{gbRecord,missing}(x, i, j, ..., drop = TRUE)

\S4method{[[}{gbRecord,ANY,ANY}(x, i, j, ...)
}
\arguments{
  \item{x}{ object from which to extract element(s) or in
  which to replace element(s).  }

  \item{i}{ indices specifying elements to extract or
  replace.  Indices are \code{numeric} or \code{character}
  vectors or empty (missing) or \code{NULL}.  Numeric
  values are coerced to integer as by
  \code{\link{as.integer}} (and hence truncated towards
  zero).  Character vectors will be matched to the
  \code{\link{names}} of the object (or for
  matrices/arrays, the \code{\link{dimnames}}): see
  \sQuote{Character indices} below for further details.

  For \code{[}-indexing only: \code{i}, \code{j},
  \code{\dots} can be logical vectors, indicating
  elements/slices to select.  Such vectors are recycled if
  necessary to match the corresponding extent.  \code{i},
  \code{j}, \code{\dots} can also be negative integers,
  indicating elements/slices to leave out of the selection.

  When indexing arrays by \code{[} a single argument
  \code{i} can be a matrix with as many columns as there
  are dimensions of \code{x}; the result is then a vector
  with elements corresponding to the sets of indices in
  each row of \code{i}.

  An index value of \code{NULL} is treated as if it were
  \code{integer(0)}.  }

  \item{drop}{For matrices and arrays.  If \code{TRUE} the
  result is coerced to the lowest possible dimension (see
  the examples).  This only works for extracting elements,
  not for the replacement.  See \code{\link{drop}} for
  further details.  }
}
\description{
See the documentation for the \code{\link[base]{Extract}}
generic, defined in the R \code{\link[base]{base-package}}
for the expected behavior.
}
\examples{
gbk_file <- system.file("extdata", "marine_metagenome.gbk", package="biofiles")
x <- gbRecord(gbk_file)

## Extract a gbFeatureList from a gbRecord:
x[1:4]

## Extract a gbFeature
x[[1]]

## Extract ggFeatures by Feature Key
x["CDS"]
}
\seealso{
\code{\link[base]{Extract}}
}

