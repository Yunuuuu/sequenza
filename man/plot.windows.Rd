\name{plotWindows}
\alias{plotWindows}
\title{Plot the specified windowed chromosome.}

\description{
  \code{plotWindows} is a convenient plot function that consent to visualize a specified chromosome attribute using the matrix produced by the function \code{windows.value}.
}

\usage{
plotWindows(abf.window, m.lty = 1, m.lwd = 3,
             m.col = "black", q.bg = "lightblue", log2.plot = FALSE,
             n.min = 1, xlim, ylim, add = FALSE, ...)
}

\arguments{
  \item{abf.window}{A matrix in the list object produced by \code{windows.value}}
  \item{m.lty}{The type of lines to use for the segments indicating the mean values.}
  \item{m.lwd}{The width of the segments indicating the mean values.}
  \item{m.col}{The color of the segments indicating the mean values.}
  \item{q.bg}{The color to use to fill the background between the 0.25 and 0.75 quartiles for each values.}
  \item{log2.plot}{Boolean arguments, it would scale the value in log2 if TRUE.}
  \item{n.min}{Integer to filter the points to plot by the minimum number of observation required.}
  \item{xlim}{\code{xlim} of the plot}
  \item{ylim}{\code{ylim} of the plot}
  \item{add}{Boolean, If TRUE the plot will be added to an existing opened device.}
  \item{...}{Any other argument accepted by the function \code{plot}.}
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
