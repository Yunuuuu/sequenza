\name{lapply_pb}

\title{A wrapper to use progress bar on multi core lapply function.}

\description{
  A wrapper to use progress bar on multi core lapply function.
}

\usage{
  lapply_pb(X, FUN, mc.cores = 2, ...
}

\arguments{
  \item{X}{X vector as used by \code{mclapply}.}
  \item{FUN}{FUN function as used by \code{mclapply}.}
  \item{mc.cores}{mc.cores as used by \code{mclapply}.}
  \item{...}{Other arguments accepted by \code{mclapply}.}
}

\details{
 Misc function borrowed from around the internet.
}