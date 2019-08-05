\name{scale_transform}
\alias{scale_transform}
\title{
  Simple Preconditioning Transformation
}
\description{
Applies the `simple' preconditioning transformation described within the `Implementation Details' section of Pollock et al. (2016) to points in the original space (\code{beta.pos}), centred at the point at which the control variate is computed (\code{beta.star}). Outputs a vector giving the position within the transformed space (\code{eta.pos}).
}
\usage{
scale_transform(beta.pos)
}
\arguments{
  \item{beta.pos}{Point within the original space to be transformed}
  \item{beta.star}{Point within the original space that the control variate is computed at}
  \item{n.sigma.inv}{Simple inverse preconditioning matrix}
}
\value{
\item{eta.pos }{Point within the transformed space}
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
un_scale_transform
}

\keyword{ Transformation }