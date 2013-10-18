\name{mutation.table}
\alias{mutation.table}
\alias{mut.fractions}
\title{Extract mutations on homozygous position from an abfreq file.}

\description{
  In an abfreq files are included the bases that diverge from the reference sample genome. The function \code{mutation.table} allow to easily extract such divergent positions.
}

\usage{
mutation.table(abf.tab, mufreq.treshold = 0.15, min.reads = 40,
               max.mut.types = 3, min.type.freq = 0.9, segments = NULL)
mut.fractions(AB.sample, Af)                
}

\arguments{
  \item{abf.tab}{An abfreq file.}
  \item{mufreq.treshold}{Mutation frequency threshold.}
  \item{min.reads}{Minimal number of reads above the quality threshold to accept the mutation call.}
  \item{max.mut.types}{Maximum number of mutation type per position. Integer from 1 to 3 (since the base are only 4). Default 3, accept "noisy" mutation calls.}
  \item{min.type.freq}{Minimal frequency of aberrant types.}
  \item{segments}{If specified, the values of depth ratio would be taken from the segments instead that from the raw data.}
  \item{AB.sample}{Base composition of the aberrant sample.}
  \item{Af}{A-allele frequency.}
}

\details{
  Call for mutation within impure tumor samples is a difficult tasks. The mutation frequency would be different comparing samples with different contamination rate. When the normal tissue is the major component of the sample (highly impure sample) it will even be difficult to distinguish a mutation from the sequencing error.
  
  The function \code{mutation.table} allow to set a threshold to extract the highest probable mutations, trying to separate the real mutation, from the sub-clonal and the sequencing error. However it would still be possible to catch sub-clonal mutation if needed (eg: infer for sub-clonal fraction in high purity samples.)
  The \code{data.frame} resulting from the function \code{mutation.table} contains some of the column of its original abfreq file, and two additional columns: \code{F} with the mutation frequency, and \code{mutation} which include 3 characters: germline nucleotide, the symbol '>' and the sample mutated nucleotide. For example, a sample that from germline \code{adednine} is mutated to \code{guanine} would be annotate as A>G.
  
  The function \code{mutation.fraction} is designed to split the different mutation types for each position (in the case of position with heterogeneous base contents), and return information on how many nucleotides were found (other then the germline nuceotide) in the column \code{base.count} and information on the major bases among the aberrant nucleotides, including \code{max.base.freq} column, indicating the frequency of the major aberrant base, considering only the aberrant reads; \code{base} indication the base of the major aberrant nucleotide; and \code{freq} which is the mutation frequency, considering the total reads in a given position.
}

\examples{

   \dontrun{

data.file <-  system.file("data", "abf.data.abfreq.txt.gz", package = "sequenza")
abf.data  <- read.abfreq(data.file)
# Detect how many reads passed the quality treshold
abf.data$good.s.reads <- abf.data$depth.sample * 
                         abf.data$sample.reads.above.quality
# Normalize coverage by GC-content
gc.stats <- gc.norm(x = abf.data$depth.ratio,
                    gc = abf.data$GC.percent)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
abf.data$adjusted.ratio <- abf.data$depth.ratio / 
                           gc.vect[as.character(abf.data$GC.percent)]     
# Subset mutations, apply mutation frequency treshold.
mut.tab   <- mutation.table(abf.data, mufreq.treshold = 0.15, 
                            min.reads = 40, max.mut.types = 1,
                            min.type.freq = 0.9)
mut.tab <- na.exclude(mut.tab)
   }
}
