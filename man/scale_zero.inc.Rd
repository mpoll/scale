\name{scale_zero.inc}
\alias{scale_zero.inc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scale_zero.inc(p.event.ch, p.curr.idx, p.curr, p.anc.neg, log.p.wei, ss.phiC, ss.phiC.up)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.event.ch}{
%%     ~~Describe \code{p.event.ch} here~~
}
  \item{p.curr.idx}{
%%     ~~Describe \code{p.curr.idx} here~~
}
  \item{p.curr}{
%%     ~~Describe \code{p.curr} here~~
}
  \item{p.anc.neg}{
%%     ~~Describe \code{p.anc.neg} here~~
}
  \item{log.p.wei}{
%%     ~~Describe \code{log.p.wei} here~~
}
  \item{ss.phiC}{
%%     ~~Describe \code{ss.phiC} here~~
}
  \item{ss.phiC.up}{
%%     ~~Describe \code{ss.phiC.up} here~~
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
function (p.event.ch, p.curr.idx, p.curr, p.anc.neg, log.p.wei, 
    ss.phiC, ss.phiC.up) 
{
    p.anc.neg.append <- c(p.event.ch, exp(log.p.wei[p.curr.idx]), 
        p.curr[c("t", "PhiL", "PhiU")])
    p.phi <- p.curr["PhiU"] - p.event.ch * p.curr["Delta"]
    print(sprintf("Old PhiL:\%f, Old PhiU:\%f", p.curr["PhiL"], 
        p.curr["PhiU"]))
    ss.phiC <- ss.phiC.up(p.phi, p.curr["PhiL"], p.curr["PhiU"])
    print(sprintf("New PhiL:\%f, New PhiU:\%f", ss.phiC(1)$phiL, 
        ss.phiC(1)$phiU))
    list(p.anc.neg.append = p.anc.neg.append, ss.phiC = ss.phiC)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
