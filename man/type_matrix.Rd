\name{types_matrix}
\alias{types_matrix}
\alias{baf_types_matrix}
\alias{mufreq_types_matrix}
\title{Creates a matrix of type tags}

\description{
  Type tags are a utensil to distinguish genomic positions by their copy number state, number A and B alleles and the number of mutated alleles.
  This function creates a matrix of all possible type tags, given the copy number of the normal sample and the range of possible copy numbers in the tumor sample.
}

\usage{
  baf_types_matrix(CNt_min, CNt_max, CNn = 2)
  mufreq_types_matrix(CNt_min, CNt_max, CNn = 2)
}

\arguments{
  \item{CNt_min}{minimum copy number in the tumor.}
  \item{CNt_max}{maximum copy number in the tumor.}
  \item{CNn}{copy number of the normal sample.}
}

\value{
\code{baf_types_matrix} returns a matrix with the following columns:
\item{CNn}{number of alleles in the normal sample.}
\item{CNt}{numbers of alleles in the aberrant sample}
\item{B}{number of B alleles}

\code{mufreq_types_matrix} returns a matrix with the following columns:
\item{CNn}{number of alleles in the normal sample.}
\item{CNt}{numbers of alleles in the aberrant sample}
\item{Mt}{number of mutated alleles}
}
\details{
    A type tag consists of 3 integers signifying the copy number in the normal and tumor samples and the number of B alleles (\code{baf_types_matrix}) or the mutated alleles (\code{mufreq_types_matrix}).
    The two functions return all the possible type tags that can be created within the given range of tumor copy numbers.
}

\seealso{
  \code{theoretical_mufreq}, \code{theoretical_depth_ratio}, \code{theoretical_baf}, \code{model_points}.
}

\examples{
  ## Generate matrix types from 0 to 4 copy number, being the non-tumor chromosome diploid.
  baf_types_matrix(CNt_min = 0, CNt_max = 4, CNn = 2 )

  ## Generate matrix types from 0 to 4 copy number, being the non-tumor chromosome monoploid.
  mufreq_types_matrix(CNt_min = 0, CNt_max = 4, CNn = 1 )

}
