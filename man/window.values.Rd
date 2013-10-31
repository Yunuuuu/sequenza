\name{windowValues}
\alias{windowValues}
\title{Binning the data on within windows of a specified size}

\description{
  Given value with corresponding genomic position with chromosomes and positions information, the function will return a list containing as much matrices as the chromosomes in the input vector. Each matrix would contains base-pair windows describing the chromosome, for each window the weighted mean, 25 and 75\% quartiles of the input values and the number of observation are reported.
}

\usage{
windowValues(x, positions, chromosomes, window = 1e6, overlap = 0, verbose = TRUE,
             weight = rep.int( x = 1, times = length(x)), start.coord = 1)
}
\arguments{
  \item{x}{The value to be windowed.}
  \item{positions}{base-pairs position.}
  \item{chromosomes}{Names or numbers of the chromosome.}
  \item{window}{size of windows used when plotting mean and quartile ranges of depth ratios and B-allele frequencies. Smaller windows will take more time to compute.}
  \item{overlap}{Integer defining the number of overlapping windows. Default is 0, no overlap.}
  \item{verbose}{Display information about the chromosome being processed?}
  \item{weight}{The weight to be assigned to each \code{x} value, usually related to the depth of reads.}
  \item{start.coord}{The coordinate to start computing the windows. If NULL, the first position available would be considered.}
}

\details{
  Sequencing data produces a big amount of data difficult to be handled by standard graphic representation. In addition, sample analyzed with older machines and with low or middle depth of coverage (20x to 50x), measures like read depth are subject of big variations, unrelated to biological reasons. Using \code{windowValues} it is possible to reduce either the amount of data as well the measure noise.
}

\seealso{
  \code{plot.windows}, \code{mutation.table}
}

\examples{
   \dontrun{

data.file <-  system.file("data", "abf.data.abfreq.txt.gz", package = "sequenza")
abf.data <- read.abfreq(data.file)
# 1Mb windows, each window is overlapping with 1 other adjacent window: depth ratio
abf.ratio <- windowValues(x = abf.data$depth.ratio,
                           positions = abf.data$n.base,
                           chromosomes = abf.data$chromosome, window = 1e6,
                           weight = abf.data$depth.normal, start.coord = 1, overlap = 1)
abf.hom  <- abf.data$ref.zygosity == 'hom'
abf.het  <- abf.data[!abf.hom, ]
# 1Mb windows, each window is overlapping with 1 other adjacent window: B-allele frequency
abf.bafs  <- windowValues(x = abf.het$Bf,
                           positions = abf.het$n.base,
                           chromosomes = abf.het$chromosome, window = 1e6,
                           weight = abf.het$depth.sample,
                           start.coord = 1, overlap = 1)
   }
}
