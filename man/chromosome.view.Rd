\name{chromosome.view}
\alias{chromosome.view}
\alias{genome.view}

\title{Graphically represent multiple chromosomal features in parallel}

\description{
  Graphically represent multiple chromosomal features in parallel.
}

\usage{
chromosome.view(baf.windows, ratio.windows, mut.tab = NULL, segments = NULL,
                min.N.baf = 1, min.N.ratio = 10000, main = "", vlines = FALSE,
                legend.inset = c(-20 * strwidth("a", units = "figure"), 0),
                BAF.style = "none", CNr = 2, cellularity = NULL, dna.index = NULL,
                avg.depth.ratio = NULL, model.lwd = 1, model.lty = "24",
                model.col = 1, x.chr.space = 10)
genome.view(seg.cn, info.type = "AB", ...)
}

\arguments{
  \item{baf.windows}{Matrix containing the windowed B-allele frequency values for one chromosome.}
  \item{ratio.windows}{Matrix containing the windowed depth ratios value for one chromosome.}
  \item{mut.tab}{Mutation table of one single chromosome. If specified the mutation would be drawn in a top panel. \code{mut.tab} must be generated from the \code{mutation.table} command.}
  \item{segments}{Segmentation for one single chromosome. If specified the segment for B-allele frequency and depth ratio would be displayed with a red line.}
  \item{min.N.baf}{Minimum amount of observations in a BAF window to be plotted.}
  \item{min.N.ratio}{Minimum amount of observations in a depth ratio window to be plotted.}
  \item{CNr}{Copy number of the germline genome. }
  \item{vlines}{Boolean if TRUE the plot would include dotted vertical lines corresponding to breaks.}
  \item{cellularity}{Value of the estimated \code{cellularity} parameter.}
  \item{dna.index}{Value of the estimated \code{dna.index} parameter.}
  \item{avg.depth.ratio}{The average value of the normalized depth ratio.}
  \item{main}{Main title of the plot.}
  \item{legend.inset}{The inset argument to pass to the \code{legend} function Define the distance of the mutation legend from the plot border.}
  \item{BAF.style}{The style to visualize the theoretical thresholds for the B allele frequency panel. Possible values are "lines", "blocks", "none". Default is "none".}
  \item{model.lwd}{Width of the theoretical lines, if the segments matrix contains the columns A, B and CNt.}
  \item{model.lty}{Line type of the theoretical lines, if the segments matrix contains the columns A, B and CNt.}
  \item{model.col}{Color of the theoretical lines, if the segments matrix contains the columns A, B and CNt.}
  \item{x.chr.space}{Step in megabase on the positions to visualize on the x-axis.}
  \item{seg.cn}{Genome wide segments, with the columns A, B and CNt.}
  \item{info.type}{Information to plot in \code{genome.view}. Available options are "CNt" to draw the copy numbers and "AB" (default) to draw the two alleles.}
  \item{...}{Optional arguments from \code{plot}.}
}

\details{
  \code{chromosome.view} is a convenient plotting function using \code{plot.windows} and some handy \code{par} parameter, to nicely display the sequencing results, and the model fitting.
}

\seealso{
  \code{windowValues}, \code{plot.windows}, \code{mutation.table}, \code{find.breaks}.
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

# Binning the values of depth ratio and B allele frequency
abf.r.win <- windowValues(x = abf.data$adjusted.ratio,
             positions = abf.data$n.base,
             chromosomes = abf.data$chromosome,
             window = 1e6, overlap = 1,
             weight = abf.data$depth.normal)

abf.b.win <- windowValues(x = abf.het$Bf,
             positions = abf.het$n.base,
             chromosomes = abf.het$chromosome,
             window = 1e6, overlap = 1,
             weight = round(x = abf.het$good.s.reads, digits = 0))
# create mutation table:
mut.tab   <- mutation.table(abf.data, mufreq.treshold = 0.15, 
                            min.reads = 40, max.mut.types = 1,
                            min.type.freq = 0.9, segments = seg.s1)
# chromosome view without parametes:
chromosome.view(mut.tab = mut.tab[mut.tab$chromosome == "1",],
                baf.windows = abf.b.win[[1]], 
                ratio.windows = abf.r.win[[1]], min.N.ratio = 1,
                segments = seg.s1[seg.s1$chromosome == "1",],
                main = "Chromosome 1")

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
                    
confint <- get.ci(CP)
dna.index   <- confint$max.x
cellularity <- confint$max.y
#detect copy number alteration on the segments:

cn.alleles <- baf.bayes(Bf = seg.s1$Bf, depth.ratio = seg.s1$depth.ratio,
                        cellularity = cellularity, dna.index = dna.index,
                        avg.depth.ratio = 1)

seg.s1 <- cbind(seg.s1, cn.alleles)

# Chromosome view with estimated paramenters:
chromosome.view(mut.tab = mut.tab[mut.tab$chromosome == "1",],
                baf.windows = abf.b.win[[1]], 
                ratio.windows = abf.r.win[[1]], min.N.ratio = 1,
                segments = seg.s1[seg.s1$chromosome == "1",],
                main = "Chromosome 1", cellularity = cellularity,
                dna.index = dna.index, avg.depth.ratio = 1,
                BAF.style = "lines")
   }
}
