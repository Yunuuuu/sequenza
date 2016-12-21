\name{theoretical_baf}
\alias{theoretical_baf}
\alias{theoretical_mufreq}
\alias{theoretical_depth_ratio}

\title{Calculates cellularity-dependent model points}

\description{
  Calculates the theoretically expected values of BAF, mutation frequency or depth ratio for given values of cellularity, DNA-content and copy number.
}

\usage{
theoretical_baf(CNt, CNn = 2, cellularity)
theoretical_depth_ratio(CNt, CNn = 2, cellularity, ploidy,
                        normal_ploidy = 2, avg_depth_ratio = 1)
theoretical_mufreq(Mt, CNt, CNn = 2, cellularity)

}

\arguments{
  \item{CNn}{copy number in the normal sample.}
  \item{CNt}{copy number in the tumor sample.}
  \item{Mt}{number of alleles carrying the mutation.}
  \item{cellularity}{fraction of tumor cells in the sample.}
  \item{ploidy}{2 * ratio between total DNA content in a tumor cell and a normal cell.}
  \item{normal_ploidy}{ploidy value in the normal sample. Default in human 2.}
  \item{avg_depth_ratio}{average normalized depth ratio.}
}

\details{
  The observed B-allele frequency, depth ratio and mutation frequency are affected by the cellularity of the tumor sample, which is the inverse of the degree of contamination by normal cells.
  Three functions are included, which for a given cellularity produce the values of B-allele frequency, mutation frequency or depth ratio that can be theoretically expected.

  \code{theoretical_baf} returns a dataframe with the possible copy numbers of A and B alleles, along with their corresponding B-allele frequency and the total copy number state (always the sum of A+B).

  \code{theoretical_depth_ratio} returns the theoretical depth ratio at a single specific position, given values of cellularity, ploidy, the ratio between the tumor copy number and the normal copy number at that position, and the average depth ratio of the sample.

  \code{theoretical_mufreq} returns the theoretical mutation frequency at a single specific position, given values of cellularity, copy number in the normal and tumor samples at that position, and the number of mutated alleles.
}

\seealso{
  \code{\link{model_points}}
}
