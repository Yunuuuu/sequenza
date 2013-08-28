\name{read.abfreq}
\alias{read.abfreq}
\alias{read.acgt}
\title{Read an abfreq or an acgt format file}

\description{
  Reading efficiently an abfreq or an acgt format file.
}

\usage{
read.abfreq(file, nrows = -1, fast = FALSE, gz = TRUE, header = TRUE,
            colClasses = c('factor', 'integer', 'factor', 'integer', 
            'integer', 'numeric', 'numeric', 'numeric', 'factor', 
            'numeric', 'numeric', "factor", "factor"), chr.name = NULL, n.lines = NULL, ...)

read.acgt(file, colClasses = c('factor', 'integer', 'factor', 'integer', 
          'integer', 'integer', 'integer', 'integer'), ...)
}

\arguments{
\item{file}{file name}
\item{nrows}{number of rows to reads from the files. Default -1, all rows.}
\item{fast}{Boolean, if true the file would be pre-parsed to count the number of rows; in some cases this can speed up the file reading.}
\item{gz}{Boolean if the input file is gzipped or not.}
\item{header}{Boolean to specify if the file have the header or not. Default TRUE.}
\item{colClasses}{Class of the columns expected in the file. Default is to match the acgt and the abfreq format.}           
\item{chr.name}{If chromosome name is specified, only the selected chromosome would be extracted instead of the whole file. Default NULL.}
\item{n.lines}{A vector of two numbes specifine first and last file line to read. If the n.lines is specified, only the selected portion of file would be used. It depends to the \code{sed} UNIX commands. Default NULL.}
\item{...}{Any arguments accepted by \code{read.delim}, or \code{read.abfreq} for the \code{read.acgt} command.}
}

\details{
  \code{abfreq} is a tab separated text file with column names headers. The current column order is still important for few functions, however most of the functions takes as argument separated column from the abfreq file. The file has currently 13 columns. The first 3 columns are derived from the original \code{pileup} file and contains:
\tabular{rll}{
1  \tab \code{chromosome} \tab with the chromosome names;\cr
2  \tab \code{n.base} \tab with the base coordinate in bp;\cr
3  \tab\code{base.ref} \tab with the base in the \code{chromosome} and \code{n.base} position of the reference genome used (usually hg19). Note the \code{base.ref} is NOT the base of the germline.
}
The remaining 10 columns contains the following information:
\tabular{rll}{
4  \tab \code{depth.normal} \tab read depth on the position observed in the normal sample;\cr
5  \tab \code{depth.sample} \tab read depth on the position observed in the aberrant sample;\cr
6  \tab \code{depth.ratio} \tab ratio of \code{depth.sample} and \code{depth.normal};\cr    
7  \tab \code{Af} \tab A-allele frequency observed in the aberrant sample;\cr
8  \tab \code{Bf} \tab B-allele frequency observed in the aberrant sample;\cr
9  \tab \code{ref.zygosity} \tab zigosity of the reference sample at the position. AA it "hom", AB and BA are indicate with "het";\cr
10  \tab \code{GC.percent} \tab genomic GC-content percent at the specific position ;\cr
11  \tab \code{sample.reads.above.quality} \tab percentage (in decimal) of reads that passed the quality threshold (threshold specified in the pre-processing software);\cr
12  \tab \code{AB.germline} \tab base of pairs found in the germline genome;\cr
13  \tab\code{AB.sample} \tab base of pairs found in the aberrant sample genome;
}

The \code{acgt} file format is similar to the \code{abfreq} file, but contains only 8 columns. The first 3 are the same as in the \code{abfreq} file. deriving from the pileup format. The remaining 5 column contain the following information:
\tabular{rll}{
4  \tab \code{read.depth} \tab read depth on the position observed in the sample. The column is deriving from the pileup file;\cr
5  \tab \code{A} \tab Number of A observed among the reads above the quality threshold;\cr
6  \tab \code{C} \tab Number of C observed among the reads above the quality threshold;\cr
7  \tab \code{G} \tab Number of G observed among the reads above the quality threshold;\cr
8  \tab \code{T} \tab Number of T observed among the reads above the quality threshold;\cr
}

}

\seealso{
  \code{read.delim}.
}

\examples{
  ## read chromosome 1 from an abfreq file.
  ## abf.data <- read.abfreq("abfreq/TCGA-42-2591_abfreq.txt.gz", header = TRUE, chr.name = "1")
}
