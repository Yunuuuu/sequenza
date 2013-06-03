\name{plot.windows}
\alias{plot.windows}
\title{Plot the specified windowed chromosome.}

\description{
  \code{plot.windows} is a wrapper to easily plot a specified chromosome included in the list object returned by the function \code{windows.value}.
}

\usage{
plot.windows(abf.window, m.lty = 1, m.lwd = 3,
             m.col = "black", q.bg = "lightblue", log2.plot = FALSE,
             n.min = 1, xlim = NULL, ylim = NULL, add = FALSE, ...)
}

\arguments{
  \item{abf.window}{A specified matrix in the list object produced by \code{windows.value}}
  \item{m.lty}{The type of lines to use for the segments indicating the mean values. Default 1.}
  \item{m.lwd}{The width of the segments indicating the mean values. Default 3.}
  \item{m.col}{The color of the segments indicating the mean values. Default "black".}
  \item{q.bg}{The color to use to fill the background between the 0.25 and 0.75 quartiles for each values. Default "lightblue"}
  \item{log2.plot}{Boolean arguments, it would scale the value in log2 if TRUE. Default FALSE.}
  \item{n.min}{Filter the points to plot by the minimum number of observation required. Integer, default 1. }
  \item{xlim}{\code{xlim} of the plot}
  \item{ylim}{\code{ylim} of the plot}
  \item{add}{Boolean, If TRUE the plot will be added to an existing opened device. Default FALSE}
  \item{...}{Any other argument accepted by the function \code{plot}.}
}

\seealso{
  \code{chromosome.view},
}

\examples{

# Use windows of 1Mb on a sample 

# abf.ratio <- window.values(x = round(abf.data$adjusted.ratio, 3),
#                           positions = abf.data$n.base,
#                           chromosomes = abf.data$chromosome, window = 1e6,
#                           overlap = 5e5, weight = abf.data$depth.normal, start.coord = 1)
#

# Plot the first of the resulting chromosome 
#

# plot.windows(abf.ratio$ windows[[1]], ylab = "Depth ratio", xlab = "Bp", 
#              main = "Chromosome 1", las =1, n.min = 1e4, ylim = c(0, 2))

# Plot all the chromosome:

# for (i in 1:lenght(abf.ratio)) {
#    plot.windows(abf.ratio$ windows[[i]], ylab = "Depth ratio", xlab = "Bp", 
#                 main = [i], las =1, n.min = 1e4, ylim = c(0, 2))
# }
}
