\name{mutation.table}
\alias{mutation.table}
\alias{mut.fractions}
\title{Extract mutations on homozygous position from an ABfreq file.}

\description{
 \code{mutation.table} extracts positions from an ABfreq file that differ from the reference genome.
}

\usage{
mutation.table(abf.tab, mufreq.treshold = 0.15, min.reads = 40,
               max.mut.types = 3, min.type.freq = 0.9, segments = NULL)
mut.fractions(AB.sample, Af)                
}

\arguments{
  \item{abf.tab}{an ABfreq table, as output from \code{read.abfreq}.}
  \item{mufreq.treshold}{mutation frequency threshold.}
  \item{min.reads}{minimal number of reads above the quality threshold to accept the mutation call.}
  \item{max.mut.types}{maximum number of different base substitutions per position. Integer from 1 to 3 (since there are only 4 different bases). Default is 3, to accept "noisy" mutation calls.}
  \item{min.type.freq}{minimal frequency of aberrant types.}
  \item{segments}{if specified, the values of depth ratio would be taken from the segments rather than from the raw data.}
  \item{AB.sample}{base composition of the aberrant sample.}
  \item{Af}{A-allele frequency.}
}

\details{
  Calling mutations in impure tumor samples is a difficult task, because the degree of contamination by normal cells affects the measured mutation frequency. In highly impure samples, where the normal cells comprise the major component of the sample, mutations can be so diluted that it will even be difficult to distinguish them from sequencing errors.
  
  The function \code{mutation.table} tries to separate true mutations from sequencing errors, based on the given threshold. In samples with low contamination, it should even be possible to catch sub-clonal mutations using this function.
  
  The function \code{mut.fractions} is designed to split the different mutation types for each position (in the case of position with heterogeneous base contents), and return information on how many nucleotides were found (other then the germline nuceotide) in the column \code{base.count} and information on the major bases among the aberrant nucleotides, including \code{max.base.freq} column, indicating the frequency of the major aberrant base, considering only the aberrant reads; \code{base} indication the base of the major aberrant nucleotide; and \code{freq} which is the mutation frequency, considering the total reads in a given position.
}

\value{
A data frame, which in addition to some of the columns of the ABfreq table, contains the following two columns:
\item{F}{the mutation frequency}
\item{mutation}{a character representation of the mutation. For example, a mutation from A in the germline to G in the tumor is annotated as "A>G".}

}
\examples{

   \dontrun{

data.file <-  system.file("data", "abf.data.abfreq.txt.gz", package = "sequenza")
abf.data  <- read.abfreq(data.file)
# Detect how many reads passed the quality treshold

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
