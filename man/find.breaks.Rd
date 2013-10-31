\name{find.breaks}
\alias{find.breaks}
\alias{segment.breaks}
\title{Segment sequencing data using allele-specific copy number algorithm}

\description{
  This function uses \code{\link{aspcf}} from the package \code{copynumber} to segment depth ratio and B-allele frequency obtained from sequencing data.
}

\usage{
find.breaks(abf.baf, gamma = 80, kmin = 10, baf.thres = c(0, 0.5), verbose = FALSE, ...)
segment.breaks(abf.tab, breaks) 
}

\arguments{
  \item{abf.baf}{An abfreq file containing only the heterozygous positions.}
  \item{abf.tab}{A complete abfreq file.}  
  \item{gamma, kmin, baf.thres, verbose}{Arguments passed to \code{aspcf}.}
  \item{breaks}{breaks obtained with the function \code{find.breaks}.}
  \item{...}{Additional arguments passed to \code{aspcf}.}
}

\details{
  \code{copynumber} is a package to perform an efficient segmentation of SNP-array data. The function \code{find.breaks} uses an allele specific algorithm to find break points. The default parameters are adjusted to behave better on sequencing data, but a careful choice of an optimal \code{gamma} value is advised.
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
abf.data  <- read.abfreq(data.file, chr.name = 12)

# Correct the coverage of the loaded chromosome:
abf.data$adjusted.ratio <- abf.data$depth.ratio / 
                           gc.vect[as.character(abf.data$GC.percent)]
# Select the heterozygous positions
abf.hom  <- abf.data$ref.zygosity == 'hom'
abf.het  <- abf.data[!abf.hom, ]
# Detect breakpoints
breaks <- find.breaks(abf.het, gamma = 80, kmin = 10, baf.thres = c(0, 0.5))
# use heterozygous and homozygous position to measure segment values
segment.breaks(abf.data, breaks = breaks)
   }
}
