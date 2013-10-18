\name{baf.model.fit}
\alias{baf.model.fit}
\alias{mufreq.model.fit}
\title{Model fitting using Bayesian inference}

\description{
  Compute the log-likelihood distribution for the specified range of cellularity and DNA-index parameters
}

\usage{
  mufreq.model.fit(cellularity = seq(0.3, 1, by = 0.01),
                   dna.index = seq(0.7, 4, by = 0.01), mc.cores = getOption("mc.cores", 2L), ...)
  baf.model.fit(cellularity = seq(0.3, 1, by = 0.01),
                dna.index = seq(0.7, 4, by = 0.01), mc.cores = getOption("mc.cores", 2L), ...)                    
}

\arguments{
  \item{cellularity}{Vector of values of cellularity parameters to be tested.}
  \item{dna.index}{Vector of values of DNA-index parameters to be tested.}
  \item{mc.cores}{number of core to use, defined as in the \code{parallel} packages.}
  \item{...}{Any argument accepted by \code{mufreq.bayes} or \code{baf.bayes}}
}

\value{
  \item{dna.index}{Values of dna.index tested}
  \item{cellularity}{Values of cellularity tested}
  \item{L}{Log-likelihood of the pair cellularity/dna.index tested}
}

\details{
  \code{baf.model.fit} uses the function \code{baf.bayes} to infer the log-likelihood of model-fitting using the possible combinations of cellularity and DNA-index values provided in the arguments.
  Similarly \code{mufreq.model.fit} fit the mutation/depth ratio model using the function \code{mufreq.bayes}.
  \code{baf.model.fit} Is generally used to infer cellularity and dna.index on segmented chromosomes, while \code{mufreq.model.fit} can be also used to estimate the same parameters from mutations, and in addition can be useful to roughly estimate sub-clonal fractions. Be aware  that the mutation frequency is more sensible to sampling bias, and is generally more noisy compared to the segmented B-allele frequency. 
}

\seealso{
  \code{baf.bayes}, \code{mufreq.bayes}, \code{cp.plot}, \code{get.ci}.
}

\examples{
  \dontrun{

data.file <-  system.file("data", "abf.data.abfreq.txt.gz", package = "sequenza")
# read all the chromosomes:
abf.data  <- read.abfreq(data.file)
# Gather genome wide GC-stats from raw file:
gc.stats <- gc.sample.stats(data.file)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
# Read only one chromosome:
abf.data  <- read.abfreq(data.file, chr.name = 1)
abf.data$good.s.reads <- abf.data$depth.sample * 
                         abf.data$sample.reads.above.quality
# Correct the coverage of the loaded chromosome:
abf.data$adjusted.ratio <- abf.data$depth.ratio / 
                           gc.vect[as.character(abf.data$GC.percent)]
# Select the heterozygous positions
abf.hom  <- abf.data$ref.zygosity == 'hom'
abf.het  <- abf.data[!abf.hom, ]
# Detect breakpoints
breaks <- find.breaks(abf.het, gamma = 80, kmin = 10, baf.thres = c(0, 0.5))
# use heterozygous and homozygous position to measure segment values
seg.s1 <- segment.breaks(abf.data, breaks = breaks)

# filter out small ambiguous segments, and conveniently weight the segments by size:
seg.filtered <- seg.s1[(seg.s1$end.pos - seg.s1$start.pos) > 10e6, ]
weights.seg  <- 150 + round((seg.filtered$end.pos - 
                             seg.filtered$start.pos) / 1e6, 0)
# get the genome wide mean of the normalized depth ratio:
avg.depth.ratio <- mean(gc.stats$adj[,2])
# run the BAF model fit

CP <- baf.model.fit(Bf = seg.filtered$Bf, depth.ratio = seg.filtered$depth.ratio, 
                    weight.ratio = weights.seg,
                    weight.Bf = weights.seg,
                    avg.depth.ratio = avg.depth.ratio,
                    cellularity = seq(0.1,1,0.01), 
                    dna.index = seq(0.5,3,0.05), mc.cores = 4)  
                    
confint <- get.ci(CP)
dna.index   <- confint$max.x
cellularity <- confint$max.y

   }                 
}
