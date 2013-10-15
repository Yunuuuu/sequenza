\name{model.points}
\alias{model.points}
\title{Generate mutation frequency and depth ratio model point at given cellularity and DNA-index value.}

\description{
  \code{model.points} function is a wrapper of various theoretic model generation functions to easily create points.
}

\usage{
model.points(cellularity = 0.5, dna.index = 1, types = cbind(CNr = 2, CNt = 2, Mt = 1),
             avg.depth.ratio = 1)}

\arguments{
  \item{cellularity}{A specified value (from 0 to 1) indicating the purity for the given tumor sample.}
  \item{dna.index}{A specified value (starting from 0.5) indicating the DNA index in the tumor sample. eg.: dna.index = 1 means diploid tumor.}
  \item{types}{desired types matrix obtained with the function \code{\link{types.matrix}}.}
  \item{avg.depth.ratio}{The median value of depth ratio.}
}
\value{
  \item{mufreq}{Theoretical value of mutation frequency at the corresponding "type tag".}
  \item{depth.ratio}{Theoretical value of depth ratio at the corresponding "type tag".}
}

\details{
  \code{model.points} generate the values of the given types tag. To learn more about the types see the manual for the function \code{\link{types.matrix}}.
}

\seealso{
  \code{\link{types.matrix}}, \code{\link{theoretical.depth.ratio}}, \code{\link{theoretical.mufreq}}.
}

\examples{

# Simulate a cellularity of 0.5, DNA index of 1 and types from min CNt 0 and max = 4 on an originally diploid chromosome:
types <- types.matrix(CNt.min = 0, CNt.max = 4, CNr = 2)
cbind(types, model.points(cellularity = 0.5, dna.index = 1,
                          types = types, avg.depth.ratio = 1))
# Simulate a cellularity of 0.5, DNA index of 1 and types from min CNt 0 and max = 4 on an originally monoallelic chromosome:
types <- types.matrix(CNt.min = 0, CNt.max = 4, CNr = 1)
cbind(types, model.points(cellularity = 0.5, dna.index = 1,
                          types = types, avg.depth.ratio = 1))

}
