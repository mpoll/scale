\name{strat.resamp}
\alias{strat.resamp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
strat.resamp(p.wei, n = length(p.wei))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.wei}{
%%     ~~Describe \code{p.wei} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (p.wei, n = length(p.wei)) 
{
    vec.length <- length(p.wei)
    cum.wei <- cumsum(p.wei)
    if (cum.wei[vec.length] > 0) {
        cum.wei <- cum.wei/cum.wei[vec.length]
        samp.vec <- seq(0, (n - 1)/n, 1/n) + runif(n, 0, 1/n)
        p.idx <- findInterval(samp.vec, cum.wei) + 1
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
