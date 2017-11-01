\name{eabesmid}
\alias{eabesmid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eabesmid(q, s, tau, t, x, m, y, minI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{minI}{
%%     ~~Describe \code{minI} here~~
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
function (q, s, tau, t, x, m, y, minI) 
{
    if (minI == -1) {
        x <- -x
        y <- -y
        m <- -m
    }
    bI <- 0
    if (q == s) {
        bI <- 1
        w <- x
    }
    if (q == tau) {
        bI <- 1
        w <- m
    }
    if (q == t) {
        bI <- 1
        w <- y
    }
    t <- t - s
    tau <- tau - s
    q <- q - s
    if (bI == 0) {
        if (q < tau) {
            Ra1 <- sqrt(tau)
            Ra2 <- (x - m) * (tau - q)/((tau)^(3/2))
            Ra3 <- (tau - q)/(tau)
        }
        else {
            Ra1 <- sqrt(t - tau)
            Ra2 <- (y - m) * (q - tau)/((t - tau)^(3/2))
            Ra3 <- (q - tau)/(t - tau)
        }
        BB3 <- rnorm(3, 0, sqrt(Ra3 * (1 - Ra3)))
        w <- m + Ra1 * sqrt((Ra2 + BB3[1])^2 + (BB3[2])^2 + (BB3[3])^2)
    }
    list(w = minI * w)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
