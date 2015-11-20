\name{sequenza.subclonal}
\alias{sequenza.subclonal}
\title{Sequenza basic functionalities to infer subclonal populations}

\description{
  These three functions are intended to be the main user interface of the package, to run several of the functions of \code{sequenza} in a standardized pipeline.
}

\usage{
  sequenza.results(sequenza.extract, cp.table = NULL, sample.id, out.dir = getwd(),
                   cellularity = NULL, ploidy = NULL, female = TRUE, CNt.max = 20,
                   ratio.priority = FALSE, XY = c(X = "X", Y = "Y"),
                   chromosome.list = 1:24)
}

\arguments{
  \item{sequenza.extract}{a list of objects as output from the \code{sequenza.extract} function.}
  \item{cp.table}{a list of objects as output from the \code{sequenza.fit} function.}
  \item{female}{logical, indicating whether the sample is male or female, to properly handle the X and Y chromosomes. Implementation only works for the human normal karyotype.}
  \item{CNt.max}{maximum copy number to consider in the model.}
  \item{XY}{character vector of length 2 specifying the labels used for the X and Y chromosomes.}
  \item{cellularity}{vector of candidate cellularity parameters.}
  \item{ploidy}{vector candidate ploidy parameters.}
  \item{ratio.priority}{logical, if TRUE only the depth ratio will be used to determine the copy number state, while the Bf value will be used to determine the number of B-alleles.}
  \item{chromosome.list}{vector containing the index or the names of the chromosome to include in the model fitting.}
  \item{sample.id}{identifier of the sample, to be used as a prefix for saved objects.}
  \item{out.dir}{output directory where the files and objects will be saved.}
}

\details{   
   The third function, \code{sequenza.results}, saves a number of objects in a specified directory (default is the working directory). The objects are:
\itemize{
 \item The list of segments with resulting copy numbers and major and minor alleles. 
 \item The candidate mutation list with variant allele frequency, and copy number and number of mutated allele, in relation of the clonal population (for sub-clonal population it needs to be processed with further methods). 
 \item A plot of all the chromosomes in one image, representing the major and minor alleles and the absolute copy number changes (genome_view).
 \item Multiple plots with one chromosome per image, representing copy-number, B-allele frequency and mutation in parallel (chromosome_view).
 \item Results of the model fitting (CP_contours and confints_CP)
 \item A summary of the copy number state of the sample (CN_bars).
}
}

\seealso{
  \code{\link{genome.view}}, \code{\link{baf.bayes}}, \code{\link{cp.plot}}, \code{\link{get.ci}}.
}

\examples{
  \dontrun{

data.file <-  system.file("data", "example.seqz.txt.gz",
              package = "sequenza")
test <- sequenza.extract(data.file)
test.CP   <- sequenza.fit(test)
sequenza.results(test, test.CP, out.dir = "example",
                 sample.id = "example")

   }
}
