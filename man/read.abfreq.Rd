\name{read.abfreq}
\alias{read.abfreq}
\alias{read.acgt}
\title{Read an abfreq or acgt format file}

\description{
  Efficiently reads an abfreq or acgt format file into R.
}

\usage{
read.abfreq(file, nrows = -1, fast = FALSE, gz = TRUE, header = TRUE,
            colClasses = c('character', 'integer', 'character', 'integer', 
            'integer', 'numeric', 'numeric', 'numeric', 'character', 
            'numeric', 'numeric', "character", "character"), chr.name = NULL, n.lines = NULL, ...)

read.acgt(file, colClasses = c('character', 'integer', 'character', 'integer', 
          'integer', 'integer', 'integer', 'integer'), ...)
}

\arguments{
\item{file}{file name}
\item{nrows}{number of rows to read from the file. Default is -1 (all rows).}
\item{fast}{logical. If TRUE the file will be pre-parsed to count the number of rows; in some cases this can speed up the file reading.}
\item{gz}{logical. If TRUE (the default) the function expects a gzipped file.}
\item{header}{a logical value indicating whether the file contains the names of the variables as its first line.}
\item{colClasses}{character.  A vector of classes to be assumed for the columns. By default the acgt and abfreq format is expected.}           
\item{chr.name}{If specified, only the selected chromosome will be extracted instead of the entire file.}
\item{n.lines}{A vector of length 2 specifying the first and last line to read from the file. If specified, only the selected portion of the file will be used. It depends to the \code{sed} UNIX commands.}
\item{...}{Any arguments accepted by \code{read.delim}. For \code{read.acgt}, also any arguments accepted by \code{read.abfreq}.}
}

\details{
  abfreq is a tab separated text file with column headers. The file has currently 13 columns. The first 3 columns are derived from the original \code{pileup} file and contain:
\tabular{rll}{
1  \tab \code{chromosome} \tab with the chromosome names;\cr
2  \tab \code{n.base} \tab with the base position;\cr
3  \tab\code{base.ref} \tab with the base in the reference genome used (usually hg19). Note the \code{base.ref} is NOT the base of the germline.
}
The remaining 10 columns contain the following information:
\tabular{rll}{
4  \tab \code{depth.normal} \tab read depth observed in the normal sample;\cr
5  \tab \code{depth.sample} \tab read depth observed in the tumor sample;\cr
6  \tab \code{depth.ratio} \tab ratio of \code{depth.sample} and \code{depth.normal};\cr    
7  \tab \code{Af} \tab A-allele frequency observed in the tumor sample;\cr
8  \tab \code{Bf} \tab B-allele frequency observed in the tumor sample;\cr
9  \tab \code{ref.zygosity} \tab zygosity of the reference sample. "hom" corresponds to AA or BB, whereas "het" corresponds to AB or BA;\cr
10  \tab \code{GC.percent} \tab GC-content (percent), calculated from the reference genome in fixed nucleotide windows ;\cr
11  \tab \code{sample.reads.above.quality} \tab fraction of reads that passed the quality threshold (threshold specified in the pre-processing software);\cr
12  \tab \code{AB.germline} \tab base found in the germline sample;\cr
13  \tab\code{AB.sample} \tab base found in the tumor sample;
}

The \code{acgt} file format is similar to the \code{abfreq} format, but contains only 8 columns. The first 3 are the same as in the \code{abfreq} file, derived from the pileup format. The remaining 5 columns contain the following information:
\tabular{rll}{
4  \tab \code{read.depth} \tab read depth. The column is derived from the pileup file;\cr
5  \tab \code{A} \tab number of times A was observed among the reads that were above the quality threshold;\cr
6  \tab \code{C} \tab number of times C was observed among the reads that were above the quality threshold;\cr
7  \tab \code{G} \tab number of times G was observed among the reads that were above the quality threshold;\cr
8  \tab \code{T} \tab number of times T was observed among the reads that were above the quality threshold;\cr
}

}

\seealso{
  \code{read.delim}.
}

\examples{
   \dontrun{

data.file <-  system.file("data", "abf.data.abfreq.txt.gz", package = "sequenza")
## read chromosome 1 from an abfreq file.
abf.data <- read.abfreq(data.file, chr.name = 1)

## fast accessing cromosome 17 using the file metrics
gc.stats <- gc.sample.stats(data.file)
chrX <- gc.stats$file.metrics[gc.stats$file.metrics$chr == "X", ]
abf.data <- read.abfreq(data.file, n.lines = c(chrX$start, chrX$end))

## Comparison the running time of the two different methods.
system.time(read.abfreq(data.file, n.lines = c(chrX$start, chrX$end)))
system.time(abf.data <- read.abfreq(data.file,chr.name="X"))
   }
}
