\name{baf.dbinom}
\alias{baf.dbinom}
\alias{baf.dbpois}
\alias{mufreq.dbinom}
\alias{mufreq.dbpois}
\alias{depth.ratio.dbinom}
\alias{depth.ratio.dbpois}

\title{Binomial and Poisson distributions}

\description{
  Return the Density of probability using a Binomial or Poisson distribution for the given value and a vector of model point.
}

\usage{
  baf.dbinom(baf, baf.model, depth.t)
  baf.dpois(baf, baf.model, depth.t)
  depth.ratio.dbinom(depth.t, depth.n, depth.ratio, depth.ratio.model)
  depth.ratio.dpois(depth.t, depth.n, depth.ratio, depth.ratio.model)
  mufreq.dbinom(mufreq, mufreq.model, depth.t, seq.errors = 0.01)
  mufreq.dpois(mufreq, mufreq.model, depth.t, seq.errors = 0.01)

}

\arguments{
  \item{baf}{The B-allele frequency in the tested position}
  \item{mufreq}{The mutation frequency in the tested position}
  \item{depth.n}{The number of reads covering the tested position in the normal sample}
  \item{depth.t}{The number of reads covering the tested position in the tumor sample}
  \item{depth.ratio}{The adjusted sequencing depth ratio in the tested position}
  \item{baf.model}{A vector of the possibles values of B-allele frequency of the model}
  \item{mufreq.model}{A vector of the possibles values of mutation frequency of the model}
  \item{depth.ratio.model}{A vector of the possible values of depth ratio of the model}
}

\details{
 Those functions are used internally by the various Bayesian model fitting function. They shouldn't be used for anything other than that.
  
  \code{baf.dbinom} and \code{baf.dpois} return the density of the probability for each values of BAF in the model, relative to a given BAF value in the tested sample.
  The function use respectively a Binomial and a Poisson distribution to compute the density values.
  \code{depth.ratio.dbinom} and \code{depth.ratio.dpois} return the density of the probability for each values of depth ratio in the model, relative to a given depth ratio value in the tested sample.
  The function use respectively a Binomial and a Poisson distribution to compute the density values.
  \code{mufreq.dbinom} and \code{mufreq.dpois} return the density of the probability for each values of mutation frequency in the model, relative to a given mutation frequency value in the tested sample.
  The function use respectively a Binomial and a Poisson distribution to compute the density values.  
}

\seealso{
  \code{theoric.baf}, \code{theoric.mufreq}, \code{theoric.depth.ratio}
}

\examples{
  ## depth.ratio.dpois(depth.t = mat[1,]$depth.sample, depth.n = mat[1,]$depth.normal, depth.ratio = mat[1,]$adjusted.ratio, test.ratio)
}
