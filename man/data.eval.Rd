\name{data.eval}
\alias{data.eval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
data.eval(beta.pos, data.idxs, factor = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta.pos}{
%%     ~~Describe \code{beta.pos} here~~
}
  \item{data.idxs}{
%%     ~~Describe \code{data.idxs} here~~
}
  \item{factor}{
%%     ~~Describe \code{factor} here~~
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
function (beta.pos, data.idxs, factor = 1) 
{
    grad.log.pi <- lap.log.pi <- numeric(dimen)
    for (i in 1:length(data.idxs)) {
        datum.call <- datum.eval(beta.pos, data.idxs[i])
        grad.log.pi <- grad.log.pi + datum.call$grad.log.pi
        lap.log.pi <- lap.log.pi + datum.call$lap.log.pi
    }
    list(grad.log.pi = factor * grad.log.pi, lap.log.pi = factor * 
        lap.log.pi)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
