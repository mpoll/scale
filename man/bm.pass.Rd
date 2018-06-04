\name{bm.pass}
\alias{bm.pass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Brownian motion passage time%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bm.pass(s = 0, x = 0, theta = 1, Jst.t = 0.64, Jst.rat = 0.5776972)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{Jst.t}{
%%     ~~Describe \code{Jst.t} here~~
}
  \item{Jst.rat}{
%%     ~~Describe \code{Jst.rat} here~~
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
function (s = 0, x = 0, theta = 1, Jst.t = 0.64, Jst.rat = 0.5776972) 
{
    repeat {
        sim <- dev.rej(dev.pr()$X)
        if (sim$Acc == 1) {
            break
        }
    }
    tau <- s + theta^2 * sim$X
    minI <- 2 * rbinom(1, 1, 0.5) - 1
    y <- x - theta * minI
    list(tau = tau, y = y, minI = minI)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
