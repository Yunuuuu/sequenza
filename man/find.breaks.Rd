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
  \code{copynumber} is a package to make efficient copy number estimation in SNP-array data. The function \code{find.breaks} uses its algorithm to find break points. The default parameters are adjusted to behave better on sequencing data, but a careful choice of better \code{gamma} value is advised.
}

\examples{
# Find breaks from heterozygous positions:
# breaks <- find.breaks(abf.het, gamma = 80, kmin = 10, baf.thres = c(0, 0.5))

# find the mean B-allele frequency and depth tratio values for the segments using all the available data. 
# seg.s1 <- segment.breaks(abf.data, breaks = breaks)
}
