\name{types.matrix}
\alias{types.matrix}
\title{Crate a data.frame of types tags}

\description{
  Types tags are quick way to identify a position by its copy number status and number of mutated alleles. this function helps on create a matrix of possible types.
}

\usage{
  types.matrix(CNt.min = 1, CNt.max = 7, CNr = 2)
}

\arguments{
  \item{CNt.min}{Minimum copy number in the tumor.}
  \item{CNt.max}{Maximum copy number in the tumor.}
  \item{CNr}{Copy number of the reference sample.}
}

\value{
\item{CNr}{Number of alleles in the normal sample.}
\item{CNt}{Numbers of allele in the aberrant sample}
\item{Mt}{Number of mutated alleles}
}
\details{
  \code{types.matrix} would return all the combinations of mutated allele from the minimum copy number state \code{CNt.min} to the maximum copy number \code{CNt.max}
}

\seealso{
  \code{theoretical.mufreq}, \code{theoretical.depth.ratio}, \code{theoretical.baf}, \code{model.points}.
}

\examples{
  ## Generate matrix types from 0 to 4 copy number, being the non-tumor chromosome diploid.
  types.matrix(CNt.min = 0, CNt.max = 4, CNr = 2 )
  
  ## Generate matrix types from 0 to 4 copy number, being the non-tumor chromosome monoploid.
  types.matrix(CNt.min = 0, CNt.max = 4, CNr = 1 )

}
