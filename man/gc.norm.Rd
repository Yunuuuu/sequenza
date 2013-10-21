\name{gc.sample.stats}
\alias{gc.sample.stats}
\alias{gc.norm}

\title{Correct the depth of coverage ratio value considering the biases derived from GC-contents.}

\description{
  Detect GC-contents driven bias, collect chromosomes coordinate in the file to improve I/O performance. 
}

\usage{
gc.sample.stats(file, gz = TRUE)
gc.norm(x, gc)
}

\arguments{
  \item{file}{A file in the ABfreq format.}
  \item{x}{A vector of values in relation with the GC-content.}
  \item{gc}{A vector corresponding to \code{x} with the relative GC-content information.}
  \item{gz}{Boolean, indicating if the file is compressed as gzip file or not. Default TRUE.}
}

\details{
  The function \code{gc.norm} given a vector of values \code{x} and a corresponding vector of GC-percentage, \code{gc.value}, performs some statistic on the relationship between the value imputed and the GC-content. The function returns a list object with two matrices: \code{raw} and \code{adj}, containing respectively the quartiles information for each GC percentage provided of the raw and the median-normalized data; and three vectors: \code{gc.values}, \code{raw.mean}, \code{raw.median}, containing respectively the values of the CG\% and corresponding mean or median value of \code{x}.
  The function \code{gc.sample.stats} reads the content of the whole file and compute the  \code{gc.norm} on the depth.ratio. In addition it compute statistics on the file contents to effectively access single chromosome reading.
  Such information is included in the list object returned by the function, in the matrix \code{file.metrics}.
}

\examples{

  \dontrun{

data.file <-  system.file("data", "abf.data.abfreq.txt.gz", package = "sequenza")
# read all the chromosomes:
abf.data  <- read.abfreq(data.file)
# Normalize coverage by GC-content
gc.stats <- gc.norm(x = abf.data$depth.ratio,
                    gc = abf.data$GC.percent)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
abf.data$adjusted.ratio <- abf.data$depth.ratio / 
                           gc.vect[as.character(abf.data$GC.percent)]

# Alternatively gather genome wide GC-stats from raw file:
gc.stats <- gc.sample.stats(data.file)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
# Read only one chromosome:
abf.data  <- read.abfreq(data.file, chr.name = 12)
# Correct the coverage of the loaded chromosome:
abf.data$adjusted.ratio <- abf.data$depth.ratio / 
                           gc.vect[as.character(abf.data$GC.percent)]

   }
}
