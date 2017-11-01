\name{ss.phiC}
\alias{ss.phiC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ss.phiC(p.mat.curr, l.bound = NULL, u.bound = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.mat.curr}{
%%     ~~Describe \code{p.mat.curr} here~~
}
  \item{l.bound}{
%%     ~~Describe \code{l.bound} here~~
}
  \item{u.bound}{
%%     ~~Describe \code{u.bound} here~~
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
function (p.mat.curr, l.bound = NULL, u.bound = NULL) 
{
    l.bds <- pmin(0, l.bound)
    u.bds <- pmax(0, u.bound)
    distance <- pmax(abs(l.bds), abs(u.bds))
    alpha.bds <- alpha.prime.bds <- numeric(dimen)
    for (j in 1:dimen) {
        alpha.bds[j] <- dsz * sum(distance * dim.grad.max[j, 
            , 1])
        alpha.prime.bds[j] <- dsz * sum(distance * dim.grad.max[j, 
            , 2])
    }
    A.bds <- alpha.bds^2
    phi.bds <- (sum(2 * abs(alpha.cent) * alpha.bds) + sum(A.bds) + 
        sum(alpha.prime.bds))/2
    phiL <- phi.cent - phi.bds
    phiU <- phi.cent + phi.bds
    intensity <- phiU - phiL
    list(distance = distance, alpha.bds = alpha.bds, alpha.prime.bds = alpha.prime.bds, 
        A.bds = A.bds, phiL = phiL, phiU = phiU, intensity = intensity)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
