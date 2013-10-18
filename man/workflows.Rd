\name{sequenza}
\alias{sequenza.extract}
\alias{sequenza.fit}
\title{Use sequenza to estimate tumor purity and ploidy.}

\description{
  Main interface of the program, aim to standardize and simplify \code{sequenza} use..
}

\usage{
  sequenza.extract(file, gz = TRUE, window = 1e6, overlap = 1, gamma = 80, kmin = 10,
                    mufreq.treshold = 0.10, min.reads = 40, max.mut.types = 1,
                    min.type.freq = 0.9)
  sequenza.fit(sequenza.extract, female = TRUE, segment.filter = 1e7, XY = c(X = "X", Y = "Y"),
           cellularity = seq(0.1,1,0.01), dna.index = seq(0.5,3,0.05),
           priors.labels = "2", priors.values = 2, mc.cores = getOption("mc.cores", 2L))                    
}

\arguments{
  \item{file}{An ABfreq file.}
  \item{gz}{Boolean indicating if the file is in the gzip format.}
  \item{window}{Size of the window to be gathered. Smaller window would take more time to be compute.}
  \item{overlap}{Integer defining the number of overlapping windows.}
  \item{gamma, kmin}{Arguments passed to \code{aspcf} from the \code{copynumber} package.}
  \item{mufreq.treshold}{Mutation frequency threshold.}
  \item{min.reads}{Minimal number of reads above the quality threshold to accept the mutation call.}
  \item{max.mut.types}{Maximum number of mutation type per position. Integer from 1 to 3 (since the base are only 4). Default 3, accept "noisy" mutation calls.}
  \item{min.type.freq}{Minimal frequency of aberrant types.}
  \item{sequenza.extract}{A list obtained from the \code{sequenza.extract} function.}
  \item{female}{Boolean indicating if the sample is male or female, to properly handle the X and Y chromosomes.}
  \item{segment.filter}{Threshold in base pairs to filter out the noise caused by short segment while testing for cellularity and ploidy. Will not effect the segmentation results.}
  \item{XY}{Vector containing the character used for the X and Y chromosome. Default c(X = "X", Y = "Y").}
  \item{cellularity}{Vector of values of cellularity parameters to be tested.}
  \item{dna.index}{Vector of values of DNA-index parameters to be tested.}
  \item{priors.labels}{Label for the copy number state priors.}
  \item{priors.values}{Value for the copy number state priors. Default priors is 1 for every copy number state.}
  \item{mc.cores}{number of core to use, defined as in the \code{parallel} packages.}
}

\details{
   The function \code{sequenza.extract} take care of efficiently read the raw data, normalize the depth.ratio considering the CG-content bias, the allele specific segmentation, filtering noisy mutations and binning the raw data. And it returns all the computed object in a single list object.
   Every object is itself a list, containing all the relative information, divided chromosome by chromosome.
}

\seealso{
  \code{genome.view}, \code{baf.bayes}, \code{cp.plot}, \code{get.ci}.
}

\examples{
  \dontrun{

data.file <-  system.file("data", "abf.data.abfreq.txt.gz",
              package = "sequenza")
test <- sequenza.extract(data.file)
CP   <- sequenza.fit(test, mc.cores = 4)
   }                 
}
