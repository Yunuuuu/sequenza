\name{get.ci}
\alias{get.ci}
\title{Compute the confidence intervals for CP matrix}

\description{
  Compute the confidence intervals for CP matrix.
}

\usage{
   get.ci(cp.table, interval = 0.95)
}
\arguments{
  \item{cp.table}{A list objects obtained with \code{baf.model.fit} or \code{mufreq/model.fit} functions.}
  \item{interval}{The minimum confidence of interval in decimal. Default 0.95}
}

\value{
  \item{values.x}{A matrix containing the dna.index values with respective likelihood.}
  \item{confint.x}{The two edges of the confidence interval for the correct dna.index estimation.}
  \item{max.x}{The point estimate for the dna.index, with the maximum likelihood.}
  \item{values.y}{A matrix containing the cellularity values with respective likelihood.}
  \item{confint.y}{The two edges of the confidence interval for the correct cellularity estimation.}
  \item{max.y}{The point estimate for the cellularity, with the maximum likelihood.}
}

\seealso{
  \code{baf.model.fit}, \code{mufreq.model.fit}.
}

\examples{
  ## read a mufreq file with the header
  ## m1 <- read.mufreq("sector1.mufreq.gz",header=T)
}
