\name{scale_resamp}
\alias{scale_resamp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scale_resamp(p.num, p.idx, log.p.wei, p.mat, p.layer, p.layer.sor.I, p.dapts, ss.size, resamp.method, ess.thresh, p.pass.arr = NULL, p.cyc.arr = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.num}{
%%     ~~Describe \code{p.num} here~~
}
  \item{p.idx}{
%%     ~~Describe \code{p.idx} here~~
}
  \item{log.p.wei}{
%%     ~~Describe \code{log.p.wei} here~~
}
  \item{p.mat}{
%%     ~~Describe \code{p.mat} here~~
}
  \item{p.layer}{
%%     ~~Describe \code{p.layer} here~~
}
  \item{p.layer.sor.I}{
%%     ~~Describe \code{p.layer.sor.I} here~~
}
  \item{p.dapts}{
%%     ~~Describe \code{p.dapts} here~~
}
  \item{ss.size}{
%%     ~~Describe \code{ss.size} here~~
}
  \item{resamp.method}{
%%     ~~Describe \code{resamp.method} here~~
}
  \item{ess.thresh}{
%%     ~~Describe \code{ess.thresh} here~~
}
  \item{p.pass.arr}{
%%     ~~Describe \code{p.pass.arr} here~~
}
  \item{p.cyc.arr}{
%%     ~~Describe \code{p.cyc.arr} here~~
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
function (p.num, p.idx, log.p.wei, p.mat, p.layer, p.layer.sor.I, 
    p.dapts, ss.size, resamp.method, ess.thresh, p.pass.arr = NULL, 
    p.cyc.arr = NULL) 
{
    wei.control <- part.ess(log.p.wei, p.num)
    if (wei.control$ess <= ess.thresh) {
        r.idx <- resamp.method(wei.control$p.wei)$p.idx
        p.idx <- p.idx[r.idx]
        p.mat <- p.mat[, r.idx, drop = FALSE]
        p.dapts[] <- sample(dsz, 3 * p.num * ss.size, replace = TRUE)
        if (is.null(p.pass.arr) == FALSE) {
            p.pass.arr <- p.pass.arr[, , r.idx]
            dimnames(p.pass.arr)[[3]] <- sprintf("part.\%i", 1:p.num)
        }
        if (is.null(p.cyc.arr) == FALSE) {
            p.cyc.arr <- p.cyc.arr[, , r.idx]
            dimnames(p.cyc.arr)[[3]] <- sprintf("part.\%i", 1:p.num)
        }
        if (p.layer.sor.I == 0) {
            p.layer <- mat.sort.r(p.layer, p.num, "c.idx")
        }
        p.layer <- p.layer[, r.idx, drop = FALSE]
        p.layer["c.idx", ] <- 1:p.num
        p.layer["t", ] <- p.layer["deg.s", ] + rexp(p.num, rate = p.layer["Delta", 
            ])
        if (is.null(p.pass.arr) == FALSE) {
            p.layer["next.ev", ] <- (sign(p.layer["tau.bar", 
                ] - p.layer["t", ]) + 1)/2
            p.layer["next.ev", ] <- p.layer["next.ev", ] * p.layer["t", 
                ] + (1 - p.layer["next.ev", ]) * p.layer["tau.bar", 
                ]
        }
        if (p.layer.sor.I == 0) {
            if (is.null(p.pass.arr) == FALSE) {
                p.layer <- mat.sort.r(p.layer, p.num, "next.ev")
            }
            else {
                p.layer <- mat.sort.r(p.layer, p.num, "t")
            }
        }
        log.p.wei[] <- log(1/p.num)
        resamp.I <- 1
    }
    else {
        r.idx <- 1:p.num
        log.p.wei <- wei.control$log.p.wei
        resamp.I <- 0
    }
    list(p.idx = p.idx, r.idx = r.idx, log.p.wei = log.p.wei, 
        p.mat = p.mat, p.layer = p.layer, p.dapts = p.dapts, 
        resamp.I = resamp.I, ess = wei.control$ess, p.pass.arr = p.pass.arr, 
        p.cyc.arr = p.cyc.arr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
