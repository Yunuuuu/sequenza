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
  \item{file}{A file in the abfreq format.}
  \item{x}{A vector of values in relation with the GC-content.}
  \item{gc}{A vector corresponding to \code{x} with the relative GC-content information.}
  \item{gz}{Boolean, indicating if the file is compressed as gzip file or not. Default TRUE.}
}

\details{
  The function \code{gc.sample.stats} reads the content of the whole file and compute the normalization.
  It returns a list containing three objects: \code{raw} and \code{adj}, containing respectively the quartiles information for each GC percentage provided of the raw and the normalized data, and \code{gc.value} containing the GC\% values. In order to normalize the data is possible to use the median of the raw data and the GC-content information for each point.
  
  The function \code{gc.norm} uses vectors already loaded in the environment. For small systems is not recommended to load an entire abfreq file in memory. it returns a list containing two objects: \code{ratio} with the normalized depth.ratio values and \code{gc.value} containing the GC\% values
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
