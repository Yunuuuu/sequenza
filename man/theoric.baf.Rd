\name{theoretical.baf}
\alias{theoretical.baf}
\alias{theoretical.mufreq}
\alias{theoretical.depth.ratio}

\title{Calculate the numeric value of the model points}

\description{
  Calculate the numeric value of the model point given the parameters of cellularity and DNA-index. 
}

\usage{
theoretical.baf(CNr, CNt, cellularity)
theoretical.depth.ratio(cellularity = 0.5, dna.index = 1, copy.number.ratio = 1.5, avg.depth.ratio = 1)
theoretical.mufreq(cellularity, CNr = 2, CNt = 2, Mt = 1)

}

\arguments{
  \item{CNr}{Copy number of the reference (usually a normal non-tumor sample). Default = 2}
  \item{CNt}{Copy number of aberrant sample (usually cancer sample). Default = 2}
  \item{Mt}{Number of alleles carrying the mutation. Default 1}
  \item{cellularity}{Tumor purity of the sample. Default 0.5}
  \item{dna.index}{Ratio of the total tumor DNA versus the total normal DNA. Default 1.}
  \item{copy.number.ratio}{Ratio of the copy number in the aberrant sample versus the normal reference copy number. Default 1.}
  \item{avg.depth.ratio}{median or mean value of the normalized depth ratios. Default 1.}
}

\details{
  The B-allele frequency value is depending of the normal contamination (or the tumor cellularity), and on the number of B-alleles versus the number of A-alleles. Given a value of cellularity and the total copy number of the tumor, the function \code{theoretical.baf} returns the possible configuration of A and B allele frequency and the corresponding theoretic value. The \code{data.frame} produced by the function contains four columns: \code{A} number of A-alleles, \code{B} number of B-allele, \code{BAF} the B-allele frequency value (always between 0 and 0.5), \code{CNt} the total copy-number state, it is always the sum of A+B.
  
  The ratio of the reads depth in a specific position, tumor versus normal, is proportional to the tumor copy number state. The value is depending on the normal contamination, the total DNA-index of the tumor (respect to the DNA in the normal sample. It is, in our opinion, a better definition of the "ploidy" term), and the tumor copy number.
  The function \code{theoretical.depth.ratio} return the theoretic depth ratio, given the values of \code{cellularity}, \code{dna.index}, the ratio between the tumor copy number and the normal copy number (normal copy number is 2 for all the chromosome except X and Y in male) and the median value of depth ratio.
  
  The mutation frequency is depending on the tumor purity and on the number of alleles carrying the mutation.
  By Increasing the copy number states increases as well the number of alleles that might be mutated. So for every N copy number, we can have N+1 possible mutation frequencies (considering also 0 as also none of the alleles can be mutated..
  The function \code{theoretical.mufreq} return the value of the corresponding mutation frequency, given the values of \code{cellularity}, the copy number of the normal cell (\code{CNr}), the copy number of the tumor (\code{CNt}) at that position and the numbers of mutated alleles (\code{Mt}). 
}

\seealso{
  \code{model.points}.
}