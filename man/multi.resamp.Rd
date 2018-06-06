\name{multi.resamp}
\alias{multi.resamp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multinomial Resampling%%  ~~function to do ... ~~
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  A simple multinomial resampling algorithm; produces a list of indices
  corresponding to n independent samples with replacement from the
  categorical distribution with the provided weight vector.
}
\usage{
multi.resamp(p.wei, n = length(p.wei))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.wei}{
Normalized particle weights
}
  \item{n}{
Number of particles (defaults to the length of the weight vector)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List of particle indices.
}

\seealso{
  resid.resamp
  strat.resamp
  system.resamp
}
\examples{
weight <- runif(10)
weight <- weight / sum(weight)
resample.idx <- multi.resamp( weight)
}

