\name{cp.plot}
\alias{cp.plot}
\alias{cp.plot.contours}
\title{Plot log-likelihood for the tested values of cellularity and ploidy}

\description{
  This function uses the \code{colorgram} function from the package \code{squash} to plot log-likelihood for the tested combinations of cellularity and ploidy
}

\usage{
   cp.plot(cp.table, ...)
   cp.plot.contours(cp.table, likThresh = c(0.95),
                    col = palette(), legend.pos = "bottomright", pch = 18, ...)
}

\arguments{
  \item{cp.table}{matrix, as output from \code{\link{baf.model.fit}} or \code{\link{mufreq.model.fit}}.}
  \item{likThresh}{vector of selected likelihood tresholds.}
  \item{col}{vector of colours.}
  \item{legend.pos}{position for placing the legend.}
  \item{pch}{character used to indicate the point estimate.}
  \item{...}{additional arguments accepted by the function \code{\link{colorgram}} for \code{cp.plot}, or \code{\link{contour}} for \code{cp.plot.contours}.}
}


\seealso{
  \code{baf.model.fit}, \code{mufreq.model.fit}.
}

\examples{

data(CP.example)
cp.plot(CP.example)
cp.plot.contours(CP.example, add = TRUE)

data(CPsnp.example)
cp.plot(CPsnp.example)
cp.plot.contours(CPsnp.example, likThresh = c(0.95, 0.9999), add = TRUE)

}
