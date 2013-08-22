\name{baf.model.fit}
\alias{baf.model.fit}
\alias{mufreq.model.fit}
\title{Model fitting using Bayesian inference}

\description{
  Compute the log-likelihood distribution for the specified range of cellularity and DNA-content parameters
}

\usage{
  mufreq.model.fit(cellularity = seq(0.3, 1, by = 0.01), dna.content = seq(0.7, 4, by = 0.01), mc.cores = 2, ...)
  baf.model.fit(cellularity = seq(0.3, 1, by = 0.01), dna.content = seq(0.7, 4, by = 0.01), mc.cores = 2, ...)                    
}

\arguments{
  \item{cellularity}{Vector of values of cellularity parameters to be tested. Default = seq(0.3, 1, by = 0.01)}
  \item{dna.content}{Vector of values of DNA-content parameters to be tested. Default = seq(0.7, 4, by = 0.01)}
  \item{mc.cores}{number of core to use, defined as in the \code{parallel} packages. Default 2.}
  \item{...}{Any argument accepted by \code{mufreq.bayes} or \code{baf.bayes}}
}

\value{
  \item{dna.content}{Values of dna.content tested}
  \item{cellularity}{Values of cellularity tested}
  \item{L}{Log-likelihood of the pair cellularity/dna-content tested}
}

\details{
  \code{baf.model.fit} uses the function \code{baf.bayes} to infer the log-likelihood of model-fitting using the possible combination of cellularity and DNA-content values provided in the arguments.
  Similarly \code{mufreq.model.fit} fit the mutation/depth ratio model using the function code{mufreq.bayes}.
  \code{baf.model.fit} Is generally used to infer cellularity and dna-contents on segmented chromosomes, while \code{mufreq.model.fit} can be also used to estimate the same parameters, in addition can be useful to roughly estimate sub-clonal fractions. Be aware  that the mutation frequency is more sensible to sampling bias, and is generally more noisy compared to the segmented B-allele frequency. 
}

\seealso{
  \code{baf.bayes}, \code{mufreq.bayes}, \code{cp.plot}.
}

\examples{
# example with few segment artificially generated with cellularity 0.75 (75%) and dna.contents 1 (diploid cell)
seg.s1<- data.frame(start.pos = c(1,4e6,10e6,40e6,70e6), end.pos = c(4e6,10e6,40e6,70e6,200e6), 
                    Bf=c(0.5,0.1250,0.5,0.1250,0.5), depth.ratio = c(1,0.625,1,0.625,1), N.BAF = 300)

# Scan over different values of cellularity and DNA-content
CP <- baf.model.fit(Bf = seg.s1$Bf, depth.ratio = seg.s1$depth.ratio, weight.ratio = 300,
                    weight.Bf = 300, avg.depth.ratio = 1, cellularity = seq(0.1, 1, 0.01),
                    dna.content = seq(0.5, 3, 0.05), mc.cores = 2)                   
}
