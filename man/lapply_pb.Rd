\name{mclapplyPb}
\alias{mclapplyPb}

\title{mclapply enhanced with a progress bar}

\description{
  This function combines \code{mclapply} (from the \code{multicore} package) with \code{\link{txtProgressBar}} (from \code{utils}).
}

\usage{
  mclapplyPb(X, FUN, mc.cores = 2, ...)
}

\arguments{
  \item{X, FUN, mc.cores, ...}{Arguments passed to \code{\link{mclapply}}.}

}

\details{
 Misc function borrowed from around the internet.
}