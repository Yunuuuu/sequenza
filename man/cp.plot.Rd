\name{cp.plot}
\alias{cp.plot}
\alias{cp.plot.contours}
\title{Plot log-likelihood for the tested values of cellularity and DNA-index}

\description{
  This function uses function \code{colorgram} of the package \code{squash} to plot log-likelihood for the tested combinations of cellularity and DNA-index
}

\usage{
   cp.plot(cp.table, ...)
   cp.plot.contours(cp.table, likThresh = c(0.95),
                    col = palette(), legend.pos = 'bottomright', pch = 18, ...)          
}

\arguments{
  \item{cp.table}{A matrix obtained using \code{baf.model.fit} or \code{mufreq.model.fit}.}
  \item{likThresh}{Vector for the likelihood tresholds to select.}
  \item{col}{color palette, or vector of colours.}
  \item{legend.pos}{Position for the legent in the function \code{cp.plot.contours}.}
  \item{pch}{Charachter to use to draw the extimated point in the \code{cp.plot.contour} function.}
  \item{...}{Additional arguments accepted by the function \code{colorgram} for \code{cp.plot}, or \code{countours} for \code{cp.plot.contours}.}
}


\seealso{
  \code{baf.model.fit}, \code{mufreq.model.fit}.
}

\examples{
# example with few segment artificially generated with cellularity 0.75 (75%) and dna.indexs 1 (diploid cell)
seg.s1<- data.frame(start.pos = c(1,4e6,10e6,40e6,70e6), end.pos = c(4e6,10e6,40e6,70e6,200e6), 
                    Bf=c(0.5,0.1250,0.5,0.1250,0.5), depth.ratio = c(1,0.625,1,0.625,1), N.BAF = 300)

# Scan over different values of cellularity and DNA-index
CP <- baf.model.fit(Bf = seg.s1$Bf, depth.ratio = seg.s1$depth.ratio,
                    weight.ratio = 300, weight.Bf = 300, avg.depth.ratio = 1,
                    cellularity = seq(0.1,1,0.01), dna.index = seq(0.5,3,0.05), mc.cores = 2)
# Plot the results                    
cp.plot(cp.table = CP)
}
