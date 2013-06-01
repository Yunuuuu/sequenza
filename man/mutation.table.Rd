\name{mutation.table}
\alias{mut.fractions}
\title{Extract mutations on homozygous position from an abfreq file.}

\description{
  In an abfreq files are recorded the base that diverge from the reference sample genome. The function \code{mutation.table} allow to easily extract such divergen positions.
}

\usage{
mutation.table(abf.tab, mufreq.treshold = 0.15, min.reads = 40,
               max.mut.types = 3, min.type.freq = 0.9, segments = NULL)
mut.fractions(AB.sample, Af)                
}

\arguments{
  \item{abf.tab}{A abfreq file.}
  \item{mufreq.treshold}{Mutation frequency treshold.}
  \item{min.reads}{Minimal number of reads above the quality treshold to accept the mutation call.}
  \item{max.mut.types}{Maximum number of mutation type per position. Integer from 1 to 3 (since the base are only 4). Default 3, accept "noisy" mutation calls.}
  \item{min.type.freq}{Minimal frequency of aberrant types.}
  \item{segments}{If specified, the values would be taken from the segments instead that from the raw data.}
  \item{AB.sample}{Base composition of the aberrant sample.}
  \item{Af}{A-allele frequency.}
}

\details{
  Call for mutation in impure tumor sample is a delicate tasks. The mutation frequency can be really low for more contaminated samples have mutation frequency Sequncing data have olways a noise background.\code{mutation.table}
}

\examples{
  ## read a mufreq file with the header
  ## m1 <- read.mufreq("sector1.mufreq.gz",header=T)
}
