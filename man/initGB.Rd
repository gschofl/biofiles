\name{initGB}
\alias{initGB}
\alias{initGB-methods}
\title{Create and/or initialise a gbRecord filehash database}
\usage{
  initGB(db_name, create=FALSE, ...)
}
\arguments{
  \item{db_name}{name of the database or database object.}

  \item{create}{if \code{TRUE} a new database is created
  and initialised, else a database object is initialised
  from an existing database.}

  \item{...}{other arguments passed to methods}
}
\value{
  Returns a \sQuote{\code{gbRecord-class}} object
  inheriting from \code{\link[filehash]{filehashRDS-class}}
}
\description{
  Create and/or initialise a gbRecord filehash database
}

