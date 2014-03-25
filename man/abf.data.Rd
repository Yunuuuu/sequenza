\name{example.seqz}
\docType{data}
\alias{example.seqz}

\title{Example \dQuote{seqz} data}

\description{
 The \dQuote{seqz} file is produced by \command{sequenza-utils.py} and typically has the file extension \file{.seqz}. The data here is representative of an exome-sequenced tumor sample, such as could be obtained from TCGA. 
}

\usage{
  data(example.seqz)
}

\details{
\code{example.seqz} can be loaded in the standard R way via \code{data(example.seqz)}, or it can be read from a text file using \code{\link{read.seqz}}. The former is useful for examples and testing, whereas the latter is representative of the standard workflow.
}

\format{
A data frame with a header row and 13 columns:
\describe{
  \item{chromosome}{chromosome names}
  \item{n.base}{base positions}
  \item{base.ref}{the base in the reference genome used (usually hg19). Note that \code{base.ref} is NOT necessarily the base of the germline.}
  \item{depth.normal}{read depth observed in the normal sample}
  \item{depth.sample}{read depth observed in the tumor sample}
  \item{depth.ratio}{ratio of \code{depth.sample} and \code{depth.normal}}
  \item{Af}{A-allele frequency observed in the tumor sample}
  \item{Bf}{B-allele frequency observed in the tumor sample}
  \item{ref.zygosity}{zygosity of the reference sample. "hom" corresponds to AA or BB, whereas "het" corresponds to AB or BA}
  \item{GC.percent}{GC-content (percent), calculated from the reference genome in fixed nucleotide windows }
  \item{good.s.reads}{number of reads that passed the quality threshold (threshold specified in the pre-processing software)}
  \item{AB.germline}{base found in the germline sample}
  \item{AB.sample}{base found in the tumor sample}
  \item{sample.strand}{frequency of the variant nucleotides detected on the forward orientation.}
}
}

\keyword{datasets}
