\name{sequenza.subclonal}
\alias{sequenza.subclonal}
\title{Sequenza basic functionalities to infer subclonal populations}

\description{
  Function intended to be the main a wrapper for many subclonal population tests.
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
   \code{sequenza.subclonal}, run basic functionalities to retrive information useful for subclonal populations inference. It return the following information:
\itemize{
 \item The list of segments with resulting copy numbers and major and minor alleles, and non integer copy number value. 
 \item The candidate mutation list with variant allele frequency, and copy number and number of mutated allele, in relation of the clonal population and cancer cell fraction (ccf). 
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
sequenza.subclonal(test, test.CP, out.dir = "example",
                 sample.id = "example")

   }
}
