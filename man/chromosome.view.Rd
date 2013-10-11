\name{chromosome.view}
\alias{chromosome.view}
\alias{genome.view}

\title{Graphically represent multiple chromosomal features in parallel}

\description{
  Reading a mufreq file in an efficient way.
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
  \item{mut.tab}{If specified the mutation would be drawn in a top panel. \code{mut.tab} must be generated from the \code{mutation.table} command. Default NULL.}
  \item{segments}{If specified the segment for B-allele frequency and depth ratio would be displayed with a red line, and vertical dotted lines would be drawn corresponding to the breaks points. Default NULL.}
  \item{min.N.baf}{Minimum amount of observations in a BAF window to be plotted. Default 1.}
  \item{min.N.ratio}{Minimum amount of observations in a depth ratio window to be plotted. Default 1e4.}
  \item{CNr}{Copy number of the germline genome. Default = 2.}
  \item{vlines}{Boolean if TRUE the plot would include dotted vertical lines corresponding to breaks. Default FALSE.}
  \item{cellularity}{Value of the estimated \code{cellularity} parameter.}
  \item{dna.index}{Value of the estimated \code{dna.index} parameter.}
  \item{avg.depth.ratio}{The median value of depth ratio.}
  \item{main}{Main title of the plot. Default "".}
  \item{legend.inset}{The inset argument to pass to the \code{legend} function, to place the mutations lagen in the plot.}
  \item{BAF.style}{The style of drawing the theoretical lines for the B allele frequency panel. Possible values are "lines", "blocks", "none". Default is "none".}
  \item{model.lwd}{Width of the theoretical lines, if the segments matrix contains the columns A, B and CNt.}
  \item{model.lty}{Line type of the theoretical lines, if the segments matrix contains the columns A, B and CNt.}
  \item{model.col}{Color of the theoretical lines, if the segments matrix contains the columns A, B and CNt.}
  \item{x.chr.space}{Step in megabase on the positions to visualize on the x-axis. Default is 10.}
  \item{seg.cn}{Genome wide segments, with the columns A, B and CNt.}
  \item{info.type}{Information to plot in \code{genome.view}. Available oprioan are "CNt", for copy number and "AB" (default) for the two alleles.}
  \item{...}{Oprional arguments from \code{plot}.}
}

\details{
  \code{chromosome.view} is a versatile wrapper using \code{plot.windows} and some handy \code{par} parameter, to nicely display the sequencing results, and the model fitting.
}

\seealso{
  \code{windowValues}, \code{plot.windows}, \code{mutation.table}, \code{find.breaks}.
}

\examples{
## Plot the chromosome view for chroimosome 1

#chromosome.view(mut.tab = mut.tab[mut.tab$chromosome == "1"],
#                baf.windows = abf.bafs$windows[[1]],
#                ratio.windows = abf.ratio$windows[[1]],
#                segments = seg.s1, main = "TCGA-42-2591 chr 1",
#                cellularity = 0.85, dna.index = 1.1, avg.depth.ratio = 1)


}
