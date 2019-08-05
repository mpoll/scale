\name{un_scale_transform}
\alias{un_scale_transform}
\title{
  Reverse Simple Preconditioning Transformation
}
\description{
Applies the reverse of the `simple' preconditioning transformation described within the `Implementation Details' section of Pollock et al. (2016) to points in the transformed space (\code{eta.pos}), which are centred such that the origin maps to the point at which the control variate is computed (\code{beta.star}). Outputs a vector giving the position within the original space (\code{beta.pos}).
}
\usage{
un_scale_transform(eta.pos)
}
\arguments{
  \item{eta.pos}{Point within the transformed space to be reverse transformed to the original space}
  \item{beta.star}{Point within the original space that the control variate is computed at}
  \item{n.sigma}{Simple preconditioning matrix}
}
\value{
\item{beta.pos }{Point within the original space}
}
\references{
  Pollock, Murray, Paul Fearnhead, Adam M. Johansen, and Gareth O. Roberts. "The scalable Langevin exact algorithm: Bayesian inference for big data." arXiv preprint arXiv:1609.03436 (2016).
}
\author{
\packageAuthor{scale}
}
\note{
As discussed within Pollock et al. (2016) the simple preconditioning procedure uses a diagonal matrix. Modification to this function would be required (along with associated package functions) if the full preconditioning matrix were to be used (\code{precon} / \code{precon.inv}).
}

\seealso{
scale_transform
}
\keyword{ Transformation }% 
