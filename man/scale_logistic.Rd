\name{scale_logistic}
\alias{scale_logistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scale for a logistic data set
}
\description{
Function for running ScaLE on a logisitic data set.
}
\usage{
scale_logistic(fnm = "logistic_default.RData", p.num = 2^10, t.inc = 0.01, T.extend = 0.01, run.length = 1e+05, ss.size = 10, ss.on = TRUE, seed.default = 1, data.precompute = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnm}{File name output to be saved to.}
  \item{p.num}{Specify the number of particles.}
  \item{t.inc}{Specify the auxiliary mesh increment size.}
  \item{T.extend}{Specify the frequency (in diffusion time) the code should be saved.}
  \item{run.length}{Specify the total diffusion time the code has to be executed for. By default set to be 10.}
  \item{ss.size}{Specify the number of data points required within the sub-sampling construction. By default set to be 2.}
  \item{ss.on}{Specify whether to use sub-sampling for evaluation of the killing intensity. By default is set to TRUE (i.e. ON), otherwise set to FALSE (i.e. OFF).}
  \item{seed.default}{For reproducibility enter the seed for the random number generator. By default is set to 1.}
  \item{data.precompute}{Whether there exists a preloaded data file with example specification. By default is TRUE. Otherwise, specify the example you wish to run ("big.logistic.example", "rare.logistic.example" or  "airline.logistic.example").}
}
\value{
  \item{curr.seed}{Current state of the RNG.}
  \item{time.elapse}{Total elapsed system time.}
  \item{simn}{List containing full ScaLE output (see below for list contents).}
  \item{simn$anc.times}{Vector of auxiliary mesh times}
  \item{simn$dimen}{Dimension}
  \item{simn$ess.thresh}{Percentage effective sample size threshold for conduting resampling.}
  \item{simn$log.p.wei}{Current log particle weights.}
  \item{simn$neg.wei.mech}{Mechanism by which to deal with negative particle weights (note: testing only).}
  \item{simn$p.anc.ess}{Vector of ESS calculation at auxiliary mesh times.}
  \item{simn$p.anc.neg}{Record of particles which had negative weights (note: testing only).}
  \item{simn$p.cyc.arr}{Current particle set layers with dimensions ordered by threshold passage times. Also see simn$p.pass.arr.}
  \item{simn$p.dapts}{Current data points associated with each particle.}
  \item{simn$p.idx}{Current index of each particle.}
  \item{simn$p.layer}{Current particle set ordered by next event times.}
  \item{simn$p.mat}{Current locations of particle set.}
  \item{simn$p.mu}{Initialisation point(s) of particle set.}
  \item{simn$p.num}{Number of particles}
  \item{simn$p.pass.arr}{Current particle set layers with dimensions ordered by dimension. Also see simn$p.cyc.arr.}
  \item{simn$p.path.renew}{Function used to update the particle set from one time to the next.}
  \item{simn$scale.transform}{Function used to transform points in original parametrisation space to preconditioned.}
  \item{simn$ss.phi}{Function used to evaluate true killing intensity (phi).}
  \item{simn$ss.phiC}{Function used to evaluate for a given hypercube the dominating Poisson killing intensity.}
  \item{simn$ss.size}{Subsample size.}
  \item{simn$T.fin}{Maximum diffusion time attained.}
  \item{simn$t.inc}{Auxiliary mesh increment size.}
  \item{simn$T.start}{Minimum diffusion time attained.}
  \item{simn$theta}{Vector of theta, specifying the size of the hypercube.}
  \item{simn$un.scale.transform}{Function used to transform points in preconditioned space to original parametrisation.}

}
\references{
Pollock et al. (2016)
}
\author{
M Pollock
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{scale_logistic()}
}

