\name{window.values}
\alias{window.values}
\title{Transform the data on base-pair windows of a specified size}

\description{
  Given value with corresponding genomic position with chromosomes and positions information, the function will return a list containing as much matrices as the chromosomes in the input vector. Each matrix would contains base-pair windows describing the chromosome, for each window the weighted mean, 25 and 75\% quartiles of the input values and the number of observation are reported.
}

\usage{
window.values(x, positions, chromosomes, window = 1e6, overlap = 0,
              weight = rep.int( x = 1, times = length(x)), start.coord = NULL)}

\arguments{
  \item{x}{The value to be windowed.}
  \item{positions}{base-pairs position.}
  \item{chromosomes}{Names or numbers of the chromosome.}
  \item{window}{Size of the window to be gathered. Default 1e6. Smaller window would take more time to be compute.}
  \item{overlap}{Base-pairs overlapping between adjacent windows. Default 0, no overlap.}
  \item{weight}{The weight to be assigned to each \code{x} value. generally is related to the depth of reads.}
  \item{start.coord}{If assigned the windowing would start from that position instead that from the first available coordinate. Default NULL.}
}

\details{
  Sequencing data, even exome seq. produce a big amount of data difficult to be handled by graphic representation. Also, especially sample analyzed with older machines and with low or middle depth of coverage (20x to 70x), value like depth of reads er position have a lot of variation. Using \code{window.values} it is possible to reduce the both the amount of data and the noise, in order to easily store and graphically represents chromosome by chromosome, or chunks of chromosomes.
}

\seealso{
  \code{plot.windows}, \code{mutation.table}
}

\examples{

# 1Mb windows with 500kb of overlap of depth ratio
#abf.ratio <- window.values(x = round(abf.data$adjusted.ratio, 3),
#                           positions = abf.data$n.base,
#                           chromosomes = abf.data$chromosome, window = 1e6,
#                           weight = abf.data$depth.normal, start.coord = 1, overlap = 5e5)

# 1Mb windows with 500kb of overlap of B-allele frequency
#abf.bafs  <- window.values(x = abf.het$Bf,
#                           positions = abf.het$n.base,
#                           chromosomes = abf.het$chromosome, window = 1e6,
#                           weight = round(x = abf.het$good.s.reads, digits = 0),
#                           start.coord = 1, overlap = 5e5)
}
