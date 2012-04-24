\docType{methods}
\name{initGB}
\alias{initGB}
\title{Create and/or initialise a gbRecord filehash database}
\arguments{
  \item{db_name}{name of the database or database object.}

  \item{create}{if \code{TRUE} a new database is created
  and initialised, else a database object is initialised
  from an existing database.}

  \item{...}{other arguments passed to methods}
}
\value{
  Returns a \sQuote{\code{gbRecord}} object inheriting from
  \code{\link[filehash]{filehashRDS-class}}
}
\description{
  Create and/or initialise a gbRecord filehash database
}

