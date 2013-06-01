\name{read.mufreq}
\alias{read.mufreq}
\title{Reading a mufreq file in an efficient way}

\description{
  Reading a mufreq file in an efficient way.
}

\usage{
  read.mufreq(file, nrows = -1, gz = TRUE,...)
}

\arguments{
  \item{file}{A mufreq file.}
}

\details{
  \code{read.mufreq} uses read.delim to read in an efficient way the gzipped or uncompressed mufreq file.
}

\seealso{
  \code{read.delim},
}

\examples{
  ## read a mufreq file with the header
  ## m1 <- read.mufreq("sector1.mufreq.gz",header=T)
}
