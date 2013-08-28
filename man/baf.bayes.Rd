\name{baf.bayes}
\alias{baf.bayes}
\alias{mufreq.bayes}
\title{Fit the model with specified cellularity and DNA-content parameters}

\description{
  Given two values for cellularity and DNA content parameters, the function returns, for each given values, the copy number and relative B-allele/mutated allele with the log-likelihood of the goodness of fit.
}

\usage{
  baf.bayes(Bf, depth.ratio, cellularity, dna.content, avg.depth.ratio,
            weight.Bf = 100, weight.ratio = 100, CNt.min = 0, CNt.max = 7, CNr = 2,
            priors.labels = CNt.min:CNt.max, priors.values = 1, ratio.priority = FALSE, skew.baf = 0.95)
  mufreq.bayes(mufreq, depth.ratio, cellularity, dna.content, avg.depth.ratio,
            weight.mufreq = 100, weight.ratio = 100, CNt.min = 1, CNt.max = 7, CNr = 2,
            priors.labels = CNt.min:CNt.max, priors.values = 1)
}

\arguments{
  \item{Bf}{A vector containing B-allele frequencies (values can range from 0 to 0.5).}
  \item{mufreq}{A vector containing mutation frequencies (values can range from 0 to 1).}
  \item{depth.ratio}{A vector containing the depth ratios of tumor versus normal.}
  \item{weight.Bf}{A vector of positive integers corresponding to the weight of the corresponding B-allele frequency value. Can be a single integer to give same weight to all BAFs.}
  \item{weight.mufreq}{A vector of positive integers corresponding to the weight of the corresponding mutation frequency value. Can be a single integer to give same weight to all frequency values.}  
  \item{weight.ratio}{A vector of positive integers corresponding to the weight of the corresponding ratio value. Can be a single integer to give same weight to all the ratios.}
  \item{cellularity}{The fraction of the sample comprising tumor cells (from 0 to 1).}
  \item{dna.content}{The amount of DNA per tumor cell, relative to a normal cell. E.g. dna.content = 2 indicates a tetraploid tumor.}
  \item{avg.depth.ratio}{The median value of depth ratio.}
  \item{CNt.min}{Minimum copy number to consider in the model.}
  \item{CNt.max}{Maximum copy number to consider in the model.}
  \item{CNr}{Copy number of the normal genome. }
  \item{priors.labels}{Label for the copy number state priors.}
  \item{priors.values}{Value for the copy number state priors. Default priors is 1 for every copy number state.}
  \item{ratio.priority}{Boolean, if true only the depth ratio would be considered to determine the copy number state, while the Bf value would select the numbers of B-allele.}
  \item{skew.baf}{The observed B-allele frequency often doen not converge to the 0.5 frequency, for normal diploid position. This argument indicate the value of the percentile of all the observed B-allele frequency to skew the theoric value below 0.5. Default 0.95}
}

\value{
  \item{CNt}{Copy number of the tumor cell at the tested point}
  \item{A}{Number of A-alleles at the tested point}
  \item{B}{Number of B-alleles at the tested point}
  \item{CNr}{Copy number of the normal cell at the tested point (It is the same as the given CNr in the arguments}
  \item{Mt}{Number of mutated alleles at the tested point}
  \item{L}{Log-likelihood of model fitting at the given point}
}
\details{
  \code{baf.bayes} and \code{mufreq.depth.ratio.fit} uses a naive Bayesian approach to calculate the likelihood of fitness of the data point with the model-point resulting from the given values of cellularity and DNA-contents.
}

\seealso{
   \code{baf.model.fit}, \code{mufreq.model.fit}.
}

\examples{
# example with few segment artificially generated with cellularity 0.75 (75%) and dna.contents 1 (diploid cell)

seg.s1<- data.frame(start.pos = c(1,4e6,10e6,40e6,70e6), end.pos = c(4e6,10e6,40e6,70e6,200e6), 
                    Bf=c(0.5,0.1250,0.5,0.1250,0.5), depth.ratio = c(1,0.625,1,0.625,1), N.BAF = 300)


# Trying to fit the model with different values
baf.bayes(Bf = seg.s1$Bf, depth.ratio = seg.s1$depth.ratio,
          weight.ratio = 300, weight.Bf = 300, 
          cellularity = 0.75, dna.content = 1, avg.depth.ratio = 1)

baf.bayes(Bf = seg.s1$Bf, depth.ratio = seg.s1$depth.ratio,
          weight.ratio = 300, weight.Bf = 300,
          cellularity = 0.5, dna.content = 1, avg.depth.ratio = 1)
          
baf.bayes(Bf = seg.s1$Bf, depth.ratio = seg.s1$depth.ratio,
          weight.ratio = 300, weight.Bf = 300,
          cellularity = 0.75, dna.content = 1.5, avg.depth.ratio = 1)
}