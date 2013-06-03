\name{get.ci}
\alias{get.ci}
\title{Compute the confidence intervals for CP matrix}

\description{
  Compute the confidence intervals for CP matrix.
}

\usage{
   get.ci(mat, interval = 0.95)
}
\arguments{
  \item{mat}{a matrix contain a column with parameter values (the first OR the second column from a CP matrix) and a column with the log-likelihood (the third column of a CP matrix).}
  \item{interval}{The minimum confidence of interval in decimal. Default 0.95}
}

\value{
  \item{values}{a matrix containing the parameter values with respective likelihood (ideal for plotting the likelihood function).}
  \item{confint}{the confidence intervals}
  \item{max.l}{The value wit maximum likelihood}
  
}

\details{
  
}

\seealso{
  \code{baf.model.fit}, \code{mufreq.model.fit}.
}

\examples{
  ## read a mufreq file with the header
  ## m1 <- read.mufreq("sector1.mufreq.gz",header=T)
}
