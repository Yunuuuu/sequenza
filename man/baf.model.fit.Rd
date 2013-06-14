\name{baf.model.fit}
\alias{baf.model.fit}
\alias{mufreq.model.fit}
\title{Model fitting using Bayesian inference}

\description{
  Compute the log-likelihood distribution for the specified range of cellularity and DNA-content parameters
}

\usage{
  mufreq.model.fit(mufreq, depth.ratio, weight.mufreq = 10, weight.ratio = 10, avg.depth.ratio,
                   cellularity.range = c(0.3,0.8), dna.content.range = c(0.5,4), by.c = 0.01, 
                   by.p = 0.01, mc.cores = 2, CNt.max = 7, CNr = 2, priors.labels = 2, priors.values = 3)
  baf.model.fit(Bf, depth.ratio, weight.Bf = 10, weight.ratio = 10, cellularity.range = c(0.3,1),
                dna.content.range = c(0.7,4), by.c = 0.01, by.p = 0.01, avg.depth.ratio,
                mc.cores = 2, CNt.max = 7, CNr = 2, priors.labels = 2, priors.values = 3)
}

\arguments{
  \item{Bf}{A vector containing the B-allele frequency, values range is from 0 to 0.5}
  \item{mufreq}{A vector containing the mutation frequency, values range is from 0 to 1}
  \item{depth.ratio}{A vector containing the depth ratios of tumor versus normal}
  \item{weight.Bf}{A vector of positive integers (bigger then 0) corresponding to the weight of the corresponding B-allele frequency value. Can be a single integer to give same weight to all BAF.}
  \item{weight.mufreq}{A vector of positive integers (bigger then 0) corresponding to the weight of the corresponding mutation frequency value. Can be a single integer to give same weight to all frequency values.}  
  \item{weight.ratio}{A vector of positive integers (bigger then 0) corresponding to the weight of the corresponding ratio value. Can be a single integer to give same weight to all the ratios.}
  \item{cellularity.range}{Range of the values of the cellularity parameters to be tested. Default = c(0.3,0.8)}
  \item{dna.content.range}{Range of the values of the DNA-content parameters to be tested. Default = c(0.5,4)}
  \item{by.c}{Increment value to use to generate the vector of cellularities to test. Default = 0.01}
  \item{by.p}{Increment value to use to generate the vector of dna-contents to test. Default = 0.01}
  \item{avg.depth.ratio}{The median value of depth ratio}
  \item{mc.cores}{Number of core to use to compute the possible combinations. Default = 2}
  \item{CNt.max}{Maximum number of possible copy number. default 7.}
  \item{CNr}{Copy number status of the germline reference. Default 2.}
  \item{priors.labels}{Label for the copy number state priors.}
  \item{priors.values}{Value for the copynumber state priors. By default there is a increased priors on the diploid status.}
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
CP <- baf.model.fit(Bf = seg.s1$Bf, depth.ratio = seg.s1$depth.ratio, 
                    weight.ratio = 300, weight.Bf = 300,
                    by.p = 0.05, by.c = 0.01, avg.depth.ratio = 1,
                    cellularity.range = c(0.1,1), dna.content.range = c(0.5,3), mc.cores = 2)                   
}
