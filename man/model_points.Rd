\name{model_points}
\alias{model_points}
\alias{baf_model_points}
\alias{mufreq_model_points}
\title{Generate B-allele frequency, mutation frequency and depth ratios at given model points, cellularity and ploidy values}

\description{
  The \code{baf_model_points} and \code{mufreq_model_points} combine \code{theoretical_baf}, \code{theoretical_mufreq} and \code{theoretical_depth_ratio} to model the theoretical respective values based on given values of cellularity and ploidy.
}

\usage{
  baf_model_points(cellularity, ploidy, baf_types, avg_depth_ratio)
  mufreq_model_points(cellularity, ploidy, mufreq_types, avg_depth_ratio)
}

\arguments{
  \item{cellularity}{fraction of tumor cells in the sample.}
  \item{ploidy}{2 * ratio between total DNA content in a tumor cell and a normal cell.}
  \item{baf_types}{matrix with the sets of copy numbers and number of mutated alleles over which to model mutation frequency and depth ratio. The matrix can be generated with \code{\link{baf_types_matrix}}.}
  \item{mufreq_types}{matrix with the sets of copy numbers and number of mutated alleles over which to model mutation frequency and depth ratio. The matrix can be generated with \code{\link{mufreq_types_matrix}}.}
  \item{avg_depth_ratio}{average normalized depth ratio.}
}
\value{
For \alias{baf_model_points}
  \item{BAF}{modelled values of B-allele frequency for the given types.}
  \item{depth_ratio}{modelled values of depth ratio for the given types.}
For \alias{mufreq_model_points} a matrix of two columns:
  \item{mufreqs}{modelled values of mutation frequency for the given types.}
  \item{depth_ratio}{modelled values of depth ratio for the given types.}
}

\details{
  The \code{baf_model_points} and \code{mufreq_model_points} functions generate the theoretical
  values of B-allele frequency, mutation frequency and depth ratio for the given type tags. To
  learn more about type tags see \code{\link{types_matrix}}.
}

\seealso{
  \code{\link{types_matrix}}, \code{\link{theoretical_depth_ratio}}, \code{\link{theoretical_mufreq}}.
}

\examples{

# Simulate a cellularity of 0.5, ploidy of 2 and types from min CNt 0
# and max = 4 on an originally diploid genome:
types <- types_matrix(CNt_min = 0, CNt_max = 4, CNn = 2)
cbind(types, model_points(cellularity = 0.5, ploidy = 2,
                          types = types, avg_depth_ratio = 1))
# Simulate a cellularity of 0.5, ploidy of 2 and types from min CNt 0
# and max = 4 on an originally monoallelic genome:
types <- types_matrix(CNt_min = 0, CNt_max = 4, CNn = 1)
cbind(types, model_points(cellularity = 0.5, ploidy = 2,
                          types = types, avg_depth_ratio = 1))

}
