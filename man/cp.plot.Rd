\name{cp.plot}
\alias{cp.plot}
\title{Plot log-likelihood for the tested values of cellularity and DNA-content}

\description{
  This function uses function \code{colorgram} of the package \code{squash} to plot log-likelihood for the tested combinations of cellularity and DNA-content
}

\usage{
   cp.plot(cp.table, map = makecmap(seq(from = median(cp.table$L, na.rm = TRUE),
                                        to = max(cp.table$L,na.rm = TRUE),
                                        by = 0.1), n = 10),
           outlier = "white", ...)                 
}

\arguments{
  \item{cp.table}{A matrix obtained using \code{baf.model.fit} or \code{mufreq.model.fit}.}
  \item{map}{\code{map} argument as used in the \code{colorgram} function. Default to display nicely the likelihood colors.}
  \item{outlier}{\code{outlier} argument as used in the \code{colorgram} function. Default value to do not color low likelihood positions.}
  \item{...}{Additional arguments accepted by the function \code{colorgram}}
}


\seealso{
  \code{baf.model.fit}, \code{mufreq.model.fit}.
}

\examples{
# example with few segment artificially generated with cellularity 0.75 (75%) and dna.contents 1 (diploid cell)
seg.s1<- data.frame(start.pos = c(1,4e6,10e6,40e6,70e6), end.pos = c(4e6,10e6,40e6,70e6,200e6), 
                    Bf=c(0.5,0.1250,0.5,0.1250,0.5), depth.ratio = c(1,0.625,1,0.625,1), N.BAF = 300)

# Scan over different values of cellularity and DNA-content
CP <- baf.model.fit(Bf = seg.s1$Bf, depth.ratio = seg.s1$depth.ratio,
                    weight.ratio = 300, weight.Bf = 300, avg.depth.ratio = 1,
                    cellularity = seq(0.1,1,0.01), dna.content = seq(0.5,3,0.05), mc.cores = 2)
# Plot the results                    
cp.plot(cp.table = CP)
}
