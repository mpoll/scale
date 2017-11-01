\name{datum.eval}
\alias{datum.eval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
datum.eval(beta.pos, data.idx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta.pos}{
%%     ~~Describe \code{beta.pos} here~~
}
  \item{data.idx}{
%%     ~~Describe \code{data.idx} here~~
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
function (beta.pos, data.idx) 
{
    datum.call <- datum(data.idx)
    datum.z <- exp(sum(datum.call$x * beta.pos))
    log.pi <- datum.call$x * beta.pos * datum.call$y - log(1 + 
        datum.z)
    grad.log.pi <- diag(n.sigma) * datum.call$x * (datum.call$y - 
        datum.z/(1 + datum.z))
    lap.log.pi <- -(diag(n.sigma) * datum.call$x)^2 * datum.z/(1 + 
        datum.z)^2
    list(log.pi = log.pi, grad.log.pi = grad.log.pi, lap.log.pi = lap.log.pi)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
