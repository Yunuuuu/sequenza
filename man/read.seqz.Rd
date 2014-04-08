\name{read.seqz}
\alias{read.seqz}
\alias{read.acgt}
\title{Read an seqz or acgt format file}

\description{
  Efficiently reads an seqz or acgt file into R.
}

\usage{
read.seqz(file, nrows = -1, fast = FALSE, gz = TRUE, header = TRUE,
            colClasses = c("character", "integer", "character", "integer",
            "integer", "numeric", "numeric", "numeric", "character",
            "numeric", "numeric", "character", "character", "character"),
            chr.name = NULL, n.lines = NULL, ...)

read.acgt(file, colClasses = c("character", "integer", "character", "integer",
          "integer", "integer", "integer", "integer", "character"), ...)
}

\arguments{
\item{file}{file name}
\item{nrows}{number of rows to read from the file. Default is -1 (all rows).}
\item{fast}{logical. If TRUE the file will be pre-parsed to count the number of rows; on some systems this can speed up the file reading.}
\item{gz}{logical. If TRUE (the default) the function expects a gzipped file.}
\item{header}{logical, indicating whether the file contains the names of the variables as its first line.}
\item{colClasses}{character.  A vector of classes to be assumed for the columns. By default the acgt and seqz format is expected.}
\item{chr.name}{if specified, only the selected chromosome will be extracted instead of the entire file.}
\item{n.lines}{vector of length 2 specifying the first and last line to read from the file. If specified, only the selected portion of the file will be used. Requires the \code{sed} UNIX utility.}
\item{...}{any arguments accepted by \code{read.delim}. For \code{read.acgt}, also any arguments accepted by \code{read.seqz}.}
}

\details{
\code{read.seqz} is a function that allows to efficiently access a file by chromosome or by number of line. The specific content of a \code{seqz} file or an \code{acgt} is explained in the \code{value} section.
}

\format{
seqz is a tab separated text file with column headers. The file has currently 14 columns. The first 3 columns are derived from the original \code{pileup} file and contain:
\describe{
\item{chromosome}{with the chromosome names}
\item{position}{with the base positions}
\item{base.ref}{with the base in the reference genome used (usually hg19). Note the \code{base.ref} is NOT the base of the germline.}
The remaining 10 columns contain the following information:
\item{depth.normal}{read depth observed in the normal sample}
\item{depth.tumor}{read depth observed in the tumor sample}
\item{depth.ratio}{ratio of \code{depth.tumor} and \code{depth.normal}}
\item{Af}{A-allele frequency observed in the tumor sample}
\item{Bf}{B-allele frequency observed in the tumor sample}
\item{zygosity.normal}{zygosity of the reference sample. "hom" corresponds to AA or BB, whereas "het" corresponds to AB or BA}
\item{GC.percent}{GC-content (percent), calculated from the reference genome in fixed nucleotide windows }
\item{good.reads}{number of reads that passed the quality threshold (threshold specified in the pre-processing software)}
\item{AB.normal}{base found in the germline sample}
\item{AB.tumor}{base found in the tumor sample}
\item{tumor.strand}{frequency of the variant nucleotides detected on the forward orientation.}
}

The \code{acgt} file format is similar to the \code{seqz} format, but contains only 8 columns. The first 3 are the same as in the \code{seqz} file, derived from the pileup format. The remaining 5 columns contain the following information:
\describe{
\item{read.depth}{read depth. The column is derived from the pileup file}
\item{A}{number of times A was observed among the reads that were above the quality threshold}
\item{C}{number of times C was observed among the reads that were above the quality threshold}
\item{G}{number of times G was observed among the reads that were above the quality threshold}
\item{T}{number of times T was observed among the reads that were above the quality threshold}
\item{strand}{string indicating the frequencies of reads in the forward strand for A, C, G and T, respectively, separated by ":".}
}
}

\seealso{
  \code{read.delim}.
}

\examples{
   \dontrun{

data.file <-  system.file("data", "example.seqz.txt.gz", package = "sequenza")

## read chromosome 1 from an seqz file.
seqz.data <- read.seqz(data.file, chr.name = 1)

## Fast access to chromosome X using the file metrics
gc.stats <- gc.sample.stats(data.file)
chrX <- gc.stats$file.metrics[gc.stats$file.metrics$chr == "X", ]
seqz.data <- read.seqz(data.file, n.lines = c(chrX$start, chrX$end))

## Compare the running time of the two different methods.
system.time(read.seqz(data.file, n.lines = c(chrX$start, chrX$end)))
system.time(seqz.data <- read.seqz(data.file,chr.name="X"))

   }
}
