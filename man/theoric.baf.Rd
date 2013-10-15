\name{theoretical.baf}
\alias{theoretical.baf}
\alias{theoretical.mufreq}
\alias{theoretical.depth.ratio}

\title{Calculate the numeric value of the model points}

\description{
  Calculate the numeric value of the model point given the parameters of cellularity and DNA-content. 
}

\usage{
theoretical.baf(CNr, CNt, cellularity)
theoretical.depth.ratio(cellularity = 0.5, dna.index = 1, copy.number.ratio = 1.5, avg.depth.ratio = 1)
theoretical.mufreq(cellularity, CNr = 2, CNt = 2, Mt = 1)

}

\arguments{
  \item{CNr}{Copy number of the reference (usually a normal non-tumor sample).}
  \item{CNt}{Copy number of aberrant sample (usually cancer sample).}
  \item{Mt}{Number of alleles carrying the mutation.}
  \item{cellularity}{Tumor fraction of the sample.}
  \item{dna.index}{Ratio of the total tumor DNA versus the total normal DNA.}
  \item{copy.number.ratio}{Ratio of the copy number in the aberrant sample versus the normal reference copy number.}
  \item{avg.depth.ratio}{The average value of the normalized depth ratio.}
}

\details{
  The B-allele frequency value is depending of the normal contamination (or tumor cellularity), and on the ratio of B-alleles and A-alleles. Given a value of cellularity and the total copy number of the tumor, the function \code{theoretical.baf} returns the possible configuration of A and B allele frequency and the corresponding theoretic value. The \code{data.frame} produced by the function contains four columns: \code{A} number of A-alleles, \code{B} number of B-alleles, \code{BAF} the B-allele frequency value (always between 0 and 0.5), \code{CNt} the total copy-number state (it is always the sum of A+B).
  
  The ratio of read depth in a specific position, tumor versus normal, is proportional to the tumor copy number state. The value is depending on the normal contamination, the ploidy of the tumor gebnome, and the local copy number.
  The function \code{theoretical.depth.ratio} return the theoretic depth ratio, given the values of \code{cellularity}, \code{dna.index}, the ratio between the tumor copy number and the normal copy number (normal copy number is 2 for all the chromosome except X and Y in male) and the average value of depth ratio.
  
  The mutation frequency is depending on the tumor purity and on the number of alleles carrying the mutation.
  By Increasing the copy number states increases as well the number of alleles that might be mutated. So for every N copy number, we can have N+1 possible mutation frequencies (considering also 0 as none of the alleles can be mutated).
  The function \code{theoretical.mufreq} return the value of the corresponding mutation frequency, given the values of \code{cellularity}, the copy number of the normal cell (\code{CNr}), the copy number of the tumor (\code{CNt}) at that position and the numbers of mutated alleles (\code{Mt}). 
}

\seealso{
  \code{model.points}.
}