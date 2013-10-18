\name{cp.plot}
\alias{cp.plot}
\alias{cp.plot.contours}
\title{Plot log-likelihood for the tested values of cellularity and DNA-index}

\description{
  This function uses the \code{colorgram} function from the package \code{squash} to plot log-likelihood for the tested combinations of cellularity and DNA-index
}

\usage{
   cp.plot(cp.table, ...)
   cp.plot.contours(cp.table, likThresh = c(0.95),
                    col = palette(), legend.pos = 'bottomright', pch = 18, ...)          
}

\arguments{
  \item{cp.table}{A matrix obtained using \code{baf.model.fit} or \code{mufreq.model.fit}.}
  \item{likThresh}{Vector containing selected likelihood tresholds.}
  \item{col}{Vector of colours.}
  \item{legend.pos}{Position for the legend in the function \code{cp.plot.contours}.}
  \item{pch}{Charachter used to draw the extimated point in the \code{cp.plot.contour} function.}
  \item{...}{Additional arguments accepted by the function \code{colorgram} for \code{cp.plot}, or \code{countours} for \code{cp.plot.contours}.}
}


\seealso{
  \code{baf.model.fit}, \code{mufreq.model.fit}.
}

\examples{
  \dontrun{

data.file <-  system.file("data", "abf.data.abfreq.txt.gz", package = "sequenza")
# read all the chromosomes:
abf.data  <- read.abfreq(data.file)
# Gather genome wide GC-stats from raw file:
gc.stats <- gc.sample.stats(data.file)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
# Read only one chromosome:
abf.data  <- read.abfreq(data.file, chr.name = 1)
abf.data$good.s.reads <- abf.data$depth.sample * 
                         abf.data$sample.reads.above.quality
# Correct the coverage of the loaded chromosome:
abf.data$adjusted.ratio <- abf.data$depth.ratio / 
                           gc.vect[as.character(abf.data$GC.percent)]
# Select the heterozygous positions
abf.hom  <- abf.data$ref.zygosity == 'hom'
abf.het  <- abf.data[!abf.hom, ]
# Detect breakpoints
breaks <- find.breaks(abf.het, gamma = 80, kmin = 10, baf.thres = c(0, 0.5))
# use heterozygous and homozygous position to measure segment values
seg.s1 <- segment.breaks(abf.data, breaks = breaks)
# filter out small ambiguous segments, and conveniently weight the segments by size:
seg.filtered <- seg.s1[(seg.s1$end.pos - seg.s1$start.pos) > 10e6, ]
weights.seg  <- 150 + round((seg.filtered$end.pos - 
                             seg.filtered$start.pos) / 1e6, 0)
# get the genome wide mean of the normalized depth ratio:
avg.depth.ratio <- mean(gc.stats$adj[,2])
# run the BAF model fit

CP <- baf.model.fit(Bf = seg.filtered$Bf, depth.ratio = seg.filtered$depth.ratio, 
                    weight.ratio = weights.seg,
                    weight.Bf = weights.seg,
                    avg.depth.ratio = avg.depth.ratio,
                    cellularity = seq(0.1,1,0.01), 
                    dna.index = seq(0.5,3,0.05), mc.cores = 4)  
                    
#plot the results
cp.plot(CP)
cp.plot.contours(CP,add = TRUE)
   }
}
