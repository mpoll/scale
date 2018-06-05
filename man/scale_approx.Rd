\name{scale_approx}
\alias{scale_approx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Approximate implementation of scale algorithm%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scale_approx(p.num, t.inc, T.fin, ss.phi, ss.phiC, dimen, scale_transform, un_scale_transform, T.start = 0, x.init = NULL, ss.size = 1, ess.thresh = 0.5, resamp.method = resid.resamp, neg.wei.mech = scale_zero.wei, prev.simn = NULL, progress.check = FALSE, phi.record = FALSE, resamp.freq = p.num - 1, theta = NULL, p.path.renew = p.path.renew)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.num}{
%%     ~~Describe \code{p.num} here~~
}
  \item{t.inc}{
%%     ~~Describe \code{t.inc} here~~
}
  \item{T.fin}{
%%     ~~Describe \code{T.fin} here~~
}
  \item{ss.phi}{
%%     ~~Describe \code{ss.phi} here~~
}
  \item{ss.phiC}{
%%     ~~Describe \code{ss.phiC} here~~
}
  \item{dimen}{
%%     ~~Describe \code{dimen} here~~
}
  \item{scale.transform}{
%%     ~~Describe \code{scale.transform} here~~
}
  \item{un.scale.transform}{
%%     ~~Describe \code{un.scale.transform} here~~
}
  \item{T.start}{
%%     ~~Describe \code{T.start} here~~
}
  \item{x.init}{
%%     ~~Describe \code{x.init} here~~
}
  \item{ss.size}{
%%     ~~Describe \code{ss.size} here~~
}
  \item{ess.thresh}{
%%     ~~Describe \code{ess.thresh} here~~
}
  \item{resamp.method}{
%%     ~~Describe \code{resamp.method} here~~
}
  \item{neg.wei.mech}{
%%     ~~Describe \code{neg.wei.mech} here~~
}
  \item{prev.simn}{
%%     ~~Describe \code{prev.simn} here~~
}
  \item{progress.check}{
%%     ~~Describe \code{progress.check} here~~
}
  \item{phi.record}{
%%     ~~Describe \code{phi.record} here~~
}
  \item{resamp.freq}{
%%     ~~Describe \code{resamp.freq} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{p.path.renew}{
%%     ~~Describe \code{p.path.renew} here~~
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
function (p.num, t.inc, T.fin, ss.phi, ss.phiC, dimen, scale.transform, 
    un.scale.transform, T.start = 0, x.init = NULL, ss.size = 1, 
    ess.thresh = 0.5, resamp.method = resid.resamp, neg.wei.mech = scale_zero.wei, 
    prev.simn = NULL, progress.check = FALSE, phi.record = FALSE, 
    resamp.freq = p.num - 1, theta = NULL, p.path.renew = p.path.renew) 
{
    p.mu <- matrix(c(rep(0, dimen)), ncol = 1, dimnames = list(sprintf("dim.\%i", 
        1:dimen), "mean"))
    p.covar <- diag(dimen)
    dimnames(p.covar) <- list(sprintf("dim.\%i", 1:dimen), sprintf("dim.\%i", 
        1:dimen))
    if (is.null(prev.simn) == TRUE) {
        if (is.null(x.init) == TRUE) {
            p.mat <- matrix(0, dimen, p.num, dimnames = list(sprintf("dim.\%i", 
                1:dimen), NULL))
        }
        else {
            if (is.vector(x.init)) {
                x.init <- matrix(x.init, length(x.init), 1)
            }
            if (dim(x.init)[2] == 1) {
                x.init <- matrix(rep(x.init, p.num), dimen, p.num)
            }
            p.mat <- apply(x.init, 2, scale.transform)
            dimnames(p.mat) <- list(sprintf("dim.\%i", 1:dimen), 
                NULL)
        }
        p.idx <- 1:p.num
        log.p.wei <- rep(log(1/p.num), p.num)
        p.layer <- matrix(0, 7, p.num, dimnames = list(c("c.idx", 
            "s", "deg.s", "t", "PhiL", "PhiU", "Delta"), NULL))
        p.layer["c.idx", ] <- p.idx
        p.dapts <- array(sample(dsz, 3 * p.num * ss.size, replace = TRUE), 
            c(p.num, ss.size, 3), dimnames = list(sprintf("part.\%i", 
                1:p.num), sprintf("ss.\%i", 1:ss.size), sprintf("dapt.\%i", 
                1:3)))
    }
    else {
        T.start <- prev.simn$T.fin
        p.idx <- prev.simn$p.idx
        log.p.wei <- prev.simn$log.p.wei
        p.mat <- prev.simn$p.mat
        p.layer <- prev.simn$p.layer
        p.dapts <- prev.simn$p.dapts
    }
    p.layer[c("PhiL", "PhiU"), ] <- apply(p.mat, 2, approx.intensity)
    p.layer[c("Delta"), ] <- p.layer[c("PhiU"), ] - p.layer[c("PhiL"), 
        ]
    p.bet.degrad.r <- p.bet.degrad <- numeric(p.num)
    p.layer["s", ] <- p.layer["deg.s", ] <- T.start
    p.layer["t", ] <- p.layer["deg.s", ] + rexp(p.num, rate = p.layer["Delta", 
        ])
    p.anc.idx <- matrix(p.idx, 1, p.num, dimnames = list(NULL, 
        sprintf("pos.\%i", 1:p.num)))
    p.anc.wei <- matrix(part.normal(log.p.wei)$p.wei, 1, p.num, 
        dimnames = list(NULL, sprintf("part.\%i", 1:p.num)))
    p.anc.mat <- array(p.mat, c(dimen, p.num, 1), dimnames = list(sprintf("dim.\%i", 
        1:dimen), sprintf("part.\%i", 1:p.num), NULL))
    p.phi.hist <- matrix(0, 0, dimen + 1, dimnames = list(NULL, 
        c("ss.phi", sprintf("dim.\%i", 1:dimen))))
    p.anc.resamp <- numeric(0)
    p.anc.ess <- numeric(0)
    p.anc.neg <- matrix(0, 0, 5, dimnames = list(NULL, c("event.ch", 
        "wt-", "t", "PhiL", "PhiU")))
    p.layer <- mat.sort.r(p.layer, p.num, "t")
    anc.times <- curr.deg.time <- T.start
    t.inc.next <- T.start + t.inc
    resamp.counter <- 1
    while (curr.deg.time < T.fin) {
        if (p.layer["t", 1] >= t.inc.next) {
            p.layer <- mat.sort.r(p.layer, p.num, "c.idx")
            p.mat <- p.mat + t(sqrt(t.inc.next - p.layer["s", 
                ]) * mvrnorm(p.num, p.mu, p.covar))
            p.bet.degrad <- cbind(p.bet.degrad, rbind(p.layer[c("c.idx"), 
                ], (t.inc.next - p.layer["deg.s", ]) * p.layer["PhiL", 
                ]))
            p.wei.agg <- aggregate(x = p.bet.degrad["degrad", 
                ], by = list(p.bet.degrad["c.idx", ]), FUN = "sum")
            log.p.wei <- log.p.wei - p.wei.agg[, "x"]
            anc.times <- c(anc.times, t.inc.next)
            p.anc.idx <- rbind(p.anc.idx, p.idx)
            p.idx <- 1:p.num
            p.normalise <- part.normal(log.p.wei)
            log.p.wei <- p.normalise$log.p.wei
            p.anc.wei <- rbind(p.anc.wei, p.normalise$p.wei)
            p.anc.mat <- abind(p.anc.mat, array(p.mat, c(dimen, 
                p.num, 1)), along = 3)
            curr.deg.time <- t.inc.next
            p.layer["deg.s", ] <- p.layer["s", ] <- curr.deg.time
            t.inc.next <- min(curr.deg.time + t.inc, T.fin)
            p.resamp <- scale_resamp(p.num, p.idx, log.p.wei, 
                p.mat, p.layer, p.layer.sor.I = 1, p.dapts, ss.size, 
                resamp.method, ess.thresh)
            p.anc.ess <- c(p.anc.ess, p.resamp$ess)
            p.idx <- p.resamp$p.idx
            log.p.wei <- p.resamp$log.p.wei
            if (p.resamp$resamp.I == 1) {
                p.anc.resamp <- c(p.anc.resamp, curr.deg.time)
                p.mat <- p.resamp$p.mat
                p.layer <- p.resamp$p.layer
                p.dapts <- p.resamp$p.dapts
            }
            p.layer <- mat.sort.r(p.layer, p.num, "t")
            resamp.counter <- 1
            if (progress.check == TRUE) {
                print(curr.deg.time)
            }
        }
        else {
            p.curr <- p.layer[, 1]
            p.curr.idx <- p.curr["c.idx"]
            p.curr.dapts <- p.dapts[p.curr.idx, , , drop = FALSE]
            p.loc.next <- p.mat[, p.curr.idx] + sqrt(p.curr["t"] - 
                p.curr["s"]) * mvrnorm(1, p.mu, p.covar)
            p.phi.eval <- ss.phi(p.loc.next, p.curr.dapts, ss.size)
            p.event.ch <- (p.curr["PhiU"] - p.phi.eval)/p.curr["Delta"]
            if (phi.record == TRUE) {
                p.phi.hist <- rbind(p.phi.hist, c(p.phi.eval, 
                  p.loc.next))
            }
            if (p.event.ch < 0) {
                p.neg <- neg.wei.mech(p.event.ch, p.curr.idx, 
                  p.curr, p.anc.neg, log.p.wei, ss.phiC, ss.phiC.up)
                log.p.wei[p.curr.idx] <- -Inf
                p.anc.neg <- rbind(p.anc.neg, p.neg$p.anc.neg.append)
                ss.phiC <- p.neg$ss.phiC
            }
            else {
                log.p.wei[p.curr.idx] <- log.p.wei[p.curr.idx] + 
                  log(p.event.ch)
            }
            p.bet.degrad <- cbind(p.bet.degrad, c(p.curr.idx, 
                (p.curr["t"] - p.curr["deg.s"]) * p.curr[c("PhiL")]))
            p.mat[, p.curr.idx] <- p.loc.next
            p.dapts[p.curr.idx, , ] <- sample(dsz, 3 * ss.size, 
                replace = TRUE)
            p.phi.bds <- ss.phiC(p.loc.next)
            p.curr[c("PhiL", "PhiU")] <- c(p.phi.bds$phiL, p.phi.bds$phiU)
            p.curr["Delta"] <- p.curr["PhiU"] - p.curr["PhiL"]
            p.curr["s"] <- p.curr["deg.s"] <- p.curr["t"]
            p.curr["t"] <- p.curr["deg.s"] + rexp(1, rate = p.curr["Delta"])
            p.layer <- ins.sort.r(p.curr, p.layer[, -1], p.num - 
                1, "t")
            if (resamp.counter\%\%(resamp.freq) == 0) {
                p.bet.degrad <- cbind(p.bet.degrad, rbind(p.layer[c("c.idx"), 
                  ], (p.curr["deg.s"] - p.layer["deg.s", ]) * 
                  p.layer["PhiL", ]))
                p.wei.agg <- aggregate(x = p.bet.degrad["degrad", 
                  ], by = list(p.bet.degrad["c.idx", ]), FUN = "sum")
                log.p.wei <- log.p.wei - p.wei.agg[, "x"]
                p.layer["deg.s", ] <- curr.deg.time <- p.curr["deg.s"]
                p.bet.degrad <- p.bet.degrad.r
                p.resamp <- scale_resamp(p.num, p.idx, log.p.wei, 
                  p.mat, p.layer, p.layer.sor.I = 0, p.dapts, 
                  ss.size, resamp.method, ess.thresh)
                p.idx <- p.resamp$p.idx
                log.p.wei <- p.resamp$log.p.wei
                if (p.resamp$resamp.I == 1) {
                  p.anc.resamp <- c(p.anc.resamp, curr.deg.time)
                  p.mat <- p.resamp$p.mat
                  p.layer <- p.resamp$p.layer
                  p.dapts <- p.resamp$p.dapts
                }
            }
            resamp.counter <- resamp.counter + 1
        }
    }
    list(p.num = p.num, p.idx = p.idx, log.p.wei = log.p.wei, 
        p.mat = p.mat, p.layer = p.layer, theta = theta, p.dapts = p.dapts, 
        p.anc.idx = p.anc.idx, p.anc.wei = p.anc.wei, p.anc.mat = p.anc.mat, 
        resamp.freq = resamp.freq, resamp.method = resamp.method, 
        neg.wei.mech = neg.wei.mech, p.anc.resamp = p.anc.resamp, 
        p.anc.ess = p.anc.ess, p.anc.neg = p.anc.neg, p.phi.hist = p.phi.hist, 
        anc.times = anc.times, T.start = T.start, t.inc = t.inc, 
        T.fin = T.fin, dimen = dimen, p.mu = p.mu, ss.size = ss.size, 
        ess.thresh = ess.thresh, ss.phi = ss.phi, ss.phiC = ss.phiC, 
        p.cyc.arr = NULL, p.pass.arr = NULL, scale.transform = scale.transform, 
        un.scale.transform = un.scale.transform, progress.check = progress.check, 
        phi.record = phi.record, p.path.renew = p.path.renew)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
