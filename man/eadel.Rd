\name{eadel}
\alias{eadel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eadel(n, s, t, x, y, m, u)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{u}{
%%     ~~Describe \code{u} here~~
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
function (n, s, t, x, y, m, u) 
{
    if (x == m) {
        xI <- 1
    }
    else {
        xI <- 0
    }
    if (y == m) {
        yI <- 1
    }
    else {
        yI <- 0
    }
    if (max(xI, yI) == 1) {
        delT <- 2
    }
    else {
        delT <- 1
    }
    if (m > max(x, y)) {
        x <- -x
        y <- -y
        m <- -m
        u <- -u
    }
    if (max(x - u, y - u, m - x, m - y) >= 0) {
        out <- 0
    }
    if (delT == 1) {
        out <- eagamma_cpp(n, s, t, x, y, m, u)/(1 - exp(-2 * 
            (x - m) * (y - m)/(t - s)))
    }
    if (delT == 2) {
        if (xI * yI == 0) {
            xoy <- max(x, y)
            out <- eadel2_cpp(n, s, t, m, xoy, u)
        }
        else {
            out <- 0
        }
    }
    if (out < 0) {
        out <- 0
    }
    if (out > 1) {
        out <- 1
    }
    if ((t - s) == 0) {
        out <- 1
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
