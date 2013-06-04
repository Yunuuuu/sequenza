\name{mutation.colors}
\alias{mutation.colors}
\title{Assign a color to mutation type}

\description{
  Function to assign a color to a mutation type, or draw a legned with mutation type and corresponding color
}

\usage{
  mutation.colors(mut.type = NULL, l.pos = "topright") 
}

\arguments{
  \item{mut.type}{Mutation type in the formata as 'A>G'. If NULL a legend would be drawn.}
  \item{l.pos}{In case of NULL \code{mut.type}, the parameter \code{l.pos} would specify the position of the legend}
}

\details{
  Graphic parameter used on \code{chromosome.view} plot function.
}