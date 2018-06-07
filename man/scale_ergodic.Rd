\name{scale_ergodic}
\alias{scale_ergodic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ScaLE Ergodic Extraction
}
\description{
Functional that takes the simulation output from running ScaLE (scale_exact or scale_extend or scale_logistic) and extracts weighted draws for use within the QSMC estimator.

User is free to specify what proportion of the output is to be retained as a fraction of the total output, or as a function of diffusion time. By default the entire simulation is extracted for use within the QSMC estimator.
}
\usage{
scale_ergodic(simn, retain.frac = NULL, retain.frac.range = NULL, trunc.time = NULL, trunc.time.range = NULL, even.weights = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simn}{Output from running scale_exact or scale_extend or scale_logistic}
  \item{retain.frac}{Specify the fraction of the simulation to be retained, otherwise do not specify}
  \item{retain.frac.range}{Specify a vector of length 2 containing the lower and upper diffusion time of the simulation to be extracted, otherwise do not specify}
  \item{trunc.time}{Specify the diffusion time from which the simulation is to be retained, otherwise do not specify}
  \item{trunc.time.range}{Specify a vector of length 2 containing the lower and upper diffusion time of the simulation to be extracted, otherwise do not specify}
  \item{even.weights}{TRUE returns a (re-sampled) particle set with equal weights. FALSE returns weighted particle set. FALSE by default}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{idxs}{Mesh indicies extracted}
  \item{times}{Vector of times extracted}
  \item{weights}{Vector of particle weights}
  \item{draws}{Matrix of particle draws (draws x dimension)}
}\references{
See Section 3.3 of Pollock et al. (2016)
}
\author{
M Pollock
}

\seealso{
scale_exact
scale_extend
scale_logistic
}
