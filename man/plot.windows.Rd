\name{plotWindows}
\alias{plotWindows}
\title{Plot the specified windowed chromosome.}

\description{
  \code{plotWindows} is a convenient plot function that visualizes a specified chromosome attribute using the matrix produced by the function \code{windowValues}.
}

\usage{
plotWindows(abf.window, m.lty = 1, m.lwd = 3,
             m.col = "black", q.bg = "lightblue", log2.plot = FALSE,
             n.min = 1, xlim, ylim, add = FALSE, ...)
}

\arguments{
  \item{abf.window}{A matrix of base-pair windos and corresponding statistics to be plotted. A list of such matrices can be output from \code{windowValues}}
  \item{m.lty}{The line type used for plotting mean values.}
  \item{m.lwd}{The line width used for plotting mean values.}
  \item{m.col}{The line color used for plotting mean values.}
  \item{q.bg}{The background color for the area between the 0.25 and 0.75 quartiles.}
  \item{log2.plot}{Boolean, if TRUE values are log2 scaled.}
  \item{n.min}{The minimum number of observations required for a point to be plotted.}
  \item{xlim}{Range of the x axis.}
  \item{ylim}{Range of the y axis.}
  \item{add}{Boolean, if TRUE the plot will be added to an existing opened device.}
  \item{...}{Any other arguments accepted by \code{plot}.}
}

\seealso{
  \code{chromosome.view},
}

\examples{

   \dontrun{

data.file <-  system.file("data", "abf.data.abfreq.txt.gz", package = "sequenza")
abf.data <- read.abfreq(data.file)
# 1Mb windows, each window is overlapping with 1 other adjacent window: depth ratio
abf.ratio <- windowValues(x = abf.data$depth.ratio,
                           positions = abf.data$n.base,
                           chromosomes = abf.data$chromosome, window = 1e6,
                           weight = abf.data$depth.normal, start.coord = 1, overlap = 1)
                           
plotWindows(abf.r.win[[1]], log2.plot = FALSE, 
            ylab = "Depth ratio", xlab = "Position (bases)", 
            main = names(abf.r.win)[1], las = 1, n.min = 1, 
            ylim = c(0, 2.5))

plotWindows(abf.r.win[[17]], log2.plot = FALSE, 
            ylab = "Depth ratio", xlab = "Position (bases)", 
            main = names(abf.r.win)[1], las = 1, n.min = 1, 
            ylim = c(0, 2.5))
   }
}
