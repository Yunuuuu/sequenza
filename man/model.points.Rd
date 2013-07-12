\name{model.points}
\alias{model.points}
\title{Function to generate mutation frequency and depth ratio model point at given cellularity and DNA-content value.}

\description{
  \code{model.points} function is a wrapper of various theoretic model generation functions to easily create points.
}

\usage{
model.points(cellularity = 0.5, dna.content = 1, types = cbind(CNr = 2, CNt = 2, Mt = 1),
             avg.depth.ratio = 1)}

\arguments{
  \item{cellularity}{A specified value (from 0 to 1) indicating the purity for the given tumor sample.}
  \item{dna.content}{A specified value (starting from 0.5) indicating the DNA-content in the tumor sample. eg.: dna.content = 1 it means diploid tumor.}
  \item{types}{desired types matrix obtained with the function \code{types.matrix}.}
  \item{avg.depth.ratio}{The median value of depth ratio}
}
\value{
  \item{mufreq}{theoretic value of mutation frequency at the corresponding "type tag".}
  \item{depth.ratio}{theoretic value of depth ratio at the corresponding "type tag"}
}

\details{
  \code{model.points} generate the values of the given types tag given. To learn more about the types look the manual for the function \code{tpes.matrix}.
}

\seealso{
  \code{types.matrix}, \code{theoretical.depth.ratio}, \code{theoretical.mufreq}.
}

\examples{
model.points()
}
