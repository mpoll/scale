\name{scale_logistic}
\alias{scale_logistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scale_logistic(fnm = "logistic_default.RData", p.num = 2^10, t.inc = 0.01, T.extend = 0.01, run.length = 1e+05, ss.size = 10, ss.on = TRUE, seed.default = 1, data.precompute = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnm}{
%%     ~~Describe \code{fnm} here~~
}
  \item{p.num}{
%%     ~~Describe \code{p.num} here~~
}
  \item{t.inc}{
%%     ~~Describe \code{t.inc} here~~
}
  \item{T.extend}{
%%     ~~Describe \code{T.extend} here~~
}
  \item{run.length}{
%%     ~~Describe \code{run.length} here~~
}
  \item{ss.size}{
%%     ~~Describe \code{ss.size} here~~
}
  \item{ss.on}{
%%     ~~Describe \code{ss.on} here~~
}
  \item{seed.default}{
%%     ~~Describe \code{seed.default} here~~
}
  \item{data.precompute}{
%%     ~~Describe \code{data.precompute} here~~
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
function (fnm = "logistic_default.RData", p.num = 2^10, t.inc = 0.01, 
    T.extend = 0.01, run.length = 1e+05, ss.size = 10, ss.on = TRUE, 
    seed.default = 1, data.precompute = TRUE) 
{
    set.seed(seed.default)
    if (data.precompute == "big.logistic.example") {
        big.logistic.example()
    }
    if (data.precompute == "rare.logistic.example") {
        rare.logistic.example()
    }
    if (data.precompute == "airline.logistic.example") {
        airline.logistic.example()
    }
    curr.seed <- .Random.seed
    save(file = fnm, dim.grad.max = dim.grad.max, alpha.cent = alpha.cent, 
        alpha.cent.sq = alpha.cent.sq, alpha.p.cent = alpha.p.cent, 
        phi.cent = phi.cent, curr.seed = curr.seed)
    time.elapse <- system.time(simn <<- scale_exact(p.num = p.num, 
        t.inc = t.inc, T.fin = T.extend, ss.phi = ss.phi, ss.phiC = ss.phiC, 
        dimen, scale.transform, un.scale.transform, T.start = 0, 
        x.init = NULL, ss.size = ss.size, ess.thresh = 0.5, resamp.method = resid.resamp, 
        neg.wei.mech = scale_zero.wei, prev.simn = NULL, progress.check = FALSE, 
        phi.record = FALSE, resamp.freq = p.num - 1, theta = NULL, 
        p.path.renew = p.path.renew))[3]
    curr.seed <- .Random.seed
    save(file = fnm, simn = simn, dim.grad.max = dim.grad.max, 
        alpha.cent = alpha.cent, alpha.cent.sq = alpha.cent.sq, 
        alpha.p.cent = alpha.p.cent, phi.cent = phi.cent, curr.seed = curr.seed, 
        time.elapse = time.elapse)
    print(time.elapse)
    while (simn$T.fin < run.length) {
        time.inc <- system.time(simn <<- scale_extend(simn, t.inc, 
            T.extend))[3]
        time.elapse <- time.elapse + time.inc
        curr.seed <- .Random.seed
        save(file = fnm, simn = simn, dim.grad.max = dim.grad.max, 
            alpha.cent = alpha.cent, alpha.cent.sq = alpha.cent.sq, 
            alpha.p.cent = alpha.p.cent, phi.cent = phi.cent, 
            curr.seed = curr.seed, time.elapse = time.elapse)
        print(time.elapse)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
