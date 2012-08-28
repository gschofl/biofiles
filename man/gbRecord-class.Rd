\name{gbRecord-class}
\alias{[[,gbRecord-method}
\alias{$,gbRecord-method}
\alias{gbRecord-class}
\alias{getFeatures,gbRecord-method}
\alias{initGB,initGB-method,gbRecord-method}
\alias{select,select-method,gbRecord-method}
\alias{shiftFeatures,gbFeature-method}
\alias{show,gbRecord-method}
\alias{write,gbRecord-method}
\title{gbRecord class}
\usage{
  .gbRecord(...)
}
\arguments{
  \item{...}{Slots of gbRecord}
}
\description{
  gbRecord is an S4 class that provides a container for
  data parsed from a GenBank flat file record. It is
  implemented as an extension of a
  '\link[filehash]{filehashRDS-class}' with no additional
  slots.
}

