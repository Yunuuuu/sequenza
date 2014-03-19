\name{VarScan2seqz}
\alias{VarScan2seqz}
\title{Import varscan output files, to be able to use it with sequenza algorithm.}

\description{
   Create an seqz file from VarScan2 output files.
}

\usage{
  VarScan2seqz(varscan.somatic, varscan.copynumber = NULL)
}

\arguments{
  \item{varscan.somatic}{data.frame from the output file of varscan somatic.}
  \item{varscan.copynumber}{data.frame from the output file of varscan copynumber. it is not required. default NULL}
}

\details{
  It is possible to run sequenza with the output of VarScan2, as version 2.3.6, but also previous release might be working as well.
  The function will create a data.frame with the seqz format, as explained in \code{\link{read.seqz}} page.
}
