\name{baf.bayes}
\alias{baf.bayes}
\alias{mufreq.bayes}
\title{Fit the model with specified cellularity and DNA-content parameters}

\description{
  Given the estimated cellularity and DNA index parameters, the function returns, for every given segment (or position) the most likely copy number and B-allele (or mutated allele) with the corresponding log-likelihood of the fit.
}

\usage{
  baf.bayes(Bf, depth.ratio, cellularity, dna.index, avg.depth.ratio,
            weight.Bf = 100, weight.ratio = 100, CNt.min = 0, CNt.max = 7, CNr = 2,
            priors.table = data.frame(CN = CNt.min:CNt.max, value = 1), ratio.priority = FALSE, skew.baf = 0.95)
  mufreq.bayes(mufreq, depth.ratio, cellularity, dna.index, avg.depth.ratio,
            weight.mufreq = 100, weight.ratio = 100, CNt.min = 1, CNt.max = 7, CNr = 2,
            priors.table = data.frame(CN = CNt.min:CNt.max, value = 1))
}

\arguments{
  \item{Bf}{A vector containing B-allele frequencies (values can range from 0 to 0.5).}
  \item{mufreq}{A vector containing mutation frequencies (values can range from 0 to 1).}
  \item{depth.ratio}{A vector containing the depth ratios of tumor versus normal.}
  \item{weight.Bf}{A vector corresponding to the weight of B-allele frequency value.}
  \item{weight.mufreq}{A vector corresponding to the weight of the mutation frequency value.}  
  \item{weight.ratio}{A vector corresponding to the weight of the ratio value.}
  \item{cellularity}{The fraction of the sample comprising tumor cells (from 0 to 1).}
  \item{dna.index}{The amount of DNA per tumor cell, relative to a normal cell. E.g. dna.index = 2 indicates a tetraploid cell.}
  \item{avg.depth.ratio}{The average value of the normalized depth ratio.}
  \item{CNt.min}{Minimum copy number to consider in the model.}
  \item{CNt.max}{Maximum copy number to consider in the model.}
  \item{CNr}{Copy number of the normal genome. }
  \item{priors.table}{A data frame containing the columns \code{CN} and \code{value}, respectively indicating the copy number and the corresponding weight. To every copy number is assigned the value 1 as default, so every values different then 1 will change the corresponding weight.}
  \item{ratio.priority}{Boolean, if true only the depth ratio would be considered to determine the copy number state, while the Bf value would select the numbers of B-allele.}
  \item{skew.baf}{The observed B-allele frequency often does not converge to the 0.5 frequency, for normal diploid positions. This argument indicate at which percentile of the observed B-allele frequency to skew the theoretic value below 0.5.}
}

\value{
  \item{CNt}{Copy number of the tumor cell at the tested point}
  \item{A}{Number of A-alleles at the tested point}
  \item{B}{Number of B-alleles at the tested point}
  \item{CNr}{Copy number of the normal cell at the tested point (It is the same as the given CNr in the arguments}
  \item{Mt}{Number of mutated alleles at the tested point}
  \item{L}{Log-likelihood of model fitting at the given point}
}
\details{
  \code{baf.bayes} and \code{mufreq.bayes} uses a naive Bayesian approach to calculate the likelihood of fitness of the data point with the model-point resulting from the given values of cellularity and DNA-contents.
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
                    dna.index = seq(0.5,3,0.05))  
                    
confint <- get.ci(CP)
dna.index   <- confint$max.x
cellularity <- confint$max.y

#detect copy number alteration on the segments:

cn.alleles <- baf.bayes(Bf = seg.s1$Bf, depth.ratio = seg.s1$depth.ratio,
                        cellularity = cellularity, dna.index = dna.index,
                        avg.depth.ratio = 1)

head(cbind(seg.s1, cn.alleles))

# create mutation table:
mut.tab   <- mutation.table(abf.data, mufreq.treshold = 0.15, 
                            min.reads = 40, max.mut.types = 1,
                            min.type.freq = 0.9, segments = seg.s1)
                            
mut.tab.clean <- na.exclude(mut.tab)

# Detect mutated alleles:
mut.alleles <- mufreq.bayes(mufreq = mut.tab.clean$F,
                            depth.ratio = mut.tab.clean$adjusted.ratio,
                            cellularity = cellularity, dna.index = dna.index,
                            avg.depth.ratio = avg.depth.ratio)
head(cbind(mut.tab.clean[,c("chromosome","n.base","F","adjusted.ratio", "mutation")], mut.alleles))


   }
}