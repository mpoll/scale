\name{resid.resamp}
\alias{resid.resamp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Residual Resampling%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
resid.resamp(p.wei, n = length(p.wei), nest.resamp = strat.resamp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.wei}{
%%     ~~Describe \code{p.wei} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{nest.resamp}{
%%     ~~Describe \code{nest.resamp} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (p.wei, n = length(p.wei), nest.resamp = strat.resamp) 
{
    if ((sum(p.wei) > 0) & (n > 0)) {
        det.resamp <- floor(n * p.wei/sum(p.wei))
        p.idx <- c(rep(1:length(p.wei), det.resamp), nest.resamp(p.wei - 
            det.resamp * (sum(p.wei)/n), n = n - sum(det.resamp))$p.idx)
    }
    else {
        p.idx <- numeric(0)
    }
    list(p.idx = p.idx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
