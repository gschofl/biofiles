% Built by roxygen3 0.1
\name{initGB}
\alias{initGB}
\title{Create and/or initialise a gbRecord filehash database}
\usage{
  initGB(db_dir, create = FALSE, ...)
}
\arguments{
  \item{db_dir}{Path to database.}

  \item{create}{If \code{TRUE} a new \code{\link{filehashRDS-class}}
instance is created and initialised, else a \code{gbRecord} instance is
initialised from an existing database.}

  \item{...}{Additional arguments passed to methods.}
}

\value{
Returns a \code{\linkS4class{gbRecord}} instance
}
\description{
Create and/or initialise a gbRecord filehash database
}
\section{Method signatures for generic initGB}{
\itemize{
  \item \linkS4class{character}
}

}

\keyword{methods}
