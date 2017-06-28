\name{read.seqz}
\alias{read.seqz}
\alias{read.acgt}
\title{Read an seqz or acgt format file}

\description{
    Efficiently reads an seqz file into R.
}

\usage{
    read.seqz(file, n_lines = NULL, gzip = TRUE, col_types = "ciciidddcddccc",
        chr_name = NULL, buffer = 33554432, parallel = 2L,
        col_names = c("chromosome", "position", "base.ref", "depth.normal",
            "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal",
            "GC.percent", "good.reads", "AB.normal", "AB.tumor",
            "tumor.strand"),...)
}

\arguments{
    \item{file}{file name}
    \item{gzip}{a logical inidicating if the input file is a gzipped
        file (default), or is a plain text file.}
    \item{col_types}{a string describing the classes of each columns of the
        input file (see \code{\link{read_tsv}}). The default value corresponds
        to the columns of a seqz file.}
    \item{chr_name}{if specified, only the selected chromosome will be
        extracted instead of the entire file. For \code{tabix} indexed files
        this argument can also be used to extract coordinated-selected genomic
        regions, eg \code{chr_name="1:1-1000000"}, will select the first
        megabase of chromosome 1.}
    \item{n_lines}{vector of length 2 specifying the first and last line to
        read from the file. If specified, only the selected portion of the
        file will be used.}
    \item{buffer}{maximal size of each chunk in bytes(see
        \code{\link{chunk.apply}}).}
    \item{parallel}{integer, number of threads used to process a seqz file
        (see \code{\link{chunk.apply}}).}
    \item{col_names}{names of the columns of the seqz file. The default
        corresponds to the columns names of a seqz file.}
    \item{...}{any arguments accepted by \code{read_tsv}.}
}

\details{
    \code{read.seqz} is a function that allows to efficiently access a
    \code{seqz} file by chromosome or by number of line. The function can also
    access coordinate scpecific regions with \code{tabix} indexed \code{seqz}
    files.
    The specific content of a \code{seqz} file is explained in the \code{value}
    section.
}

\format{
    A seqz file is a tab separated text file with column headers. The file has
    14 columns. The first 3 columns are derived from the original \code{pileup}
    file and contain:
    \describe{
        \item{chromosome}{with the chromosome name}
        \item{position}{with the base position}
        \item{base.ref}{with the base in the reference genome used.
            Note the \code{base.ref} is NOT necessarily the same base as in
            the normal specimen, but it is the nucleotide in the used genome
            assembly.}
    The remaining 10 columns contain the following information:
        \item{depth.normal}{read depth observed in the normal sample}
        \item{depth.tumor}{read depth observed in the tumor sample}
        \item{depth.ratio}{ratio of \code{depth.tumor} and \code{depth.normal}}
        \item{Af}{A-allele frequency observed in the tumor sample}
        \item{Bf}{B-allele frequency observed in the tumor sample in
            heterozygous positions}
        \item{zygosity.normal}{zygosity of the reference sample. "hom"
            corresponds to AA or BB, whereas "het" corresponds to AB or BA}
        \item{GC.percent}{GC-content (percent), calculated from the reference
            genome in fixed nucleotide windows }
        \item{good.reads}{number of reads that passed the quality threshold
            (threshold specified in the pre-processing software), in the
            tumor specimen}
        \item{AB.normal}{base(s) found in the germline sample; for heterozygous
            positions AB are sorted using the values of Af and Bf respectively}
        \item{AB.tumor}{base(s) found in the tumor sample not present in the
            normal specimen. The field include all the variants found in the
            tumor alignment, separated by a colon. Each variant contains the
            base and the observed frequency}
        \item{tumor.strand}{frequency of the variant nucleotides detected on
            the forward orientation. The field have a consistent structure with
            \code{AB.tumor}, indicating the fraction, relative to the total
            number of reads presenting the specific variant, orientated
            in the forward direction}
    }
}

\seealso{
  \code{read_delim}.
}

\examples{
   \dontrun{

    data_file <-  system.file("data", "example.seqz.txt.gz", package = "sequenza")

    ## read chromosome 1 from an seqz file.
    seqz_data <- read.seqz(data_file, chr_name = 1)

    ## Fast access to chromosome X using the file metrics
    gc.stats <- gc.sample.stats(data.file)
    chrX <- gc.stats$file.metrics[gc.stats$file.metrics$chr == "X", ]
    seqz.data <- read.seqz(data.file, n.lines = c(chrX$start, chrX$end))

    ## Compare the running time of the two different methods.
    system.time(read.seqz(data.file, n.lines = c(chrX$start, chrX$end)))
    system.time(seqz.data <- read.seqz(data.file,chr.name="X"))

   }
}
