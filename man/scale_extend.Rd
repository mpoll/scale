\name{scale_extend}
\alias{scale_extend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scale extend%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scale_extend(prev.simn, t.inc, T.extend)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prev.simn}{
%%     ~~Describe \code{prev.simn} here~~
}
  \item{t.inc}{
%%     ~~Describe \code{t.inc} here~~
}
  \item{T.extend}{
%%     ~~Describe \code{T.extend} here~~
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
function (prev.simn, t.inc, T.extend) 
{
    new.simn <- scale_exact(p.num = prev.simn$p.num, t.inc = t.inc, 
        T.fin = T.extend + prev.simn$T.fin, ss.phi = ss.phi, 
        ss.phiC = prev.simn$ss.phiC, dimen = prev.simn$dimen, 
        scale.transform = scale.transform, un.scale.transform = un.scale.transform, 
        T.start = prev.simn$T.fin, ss.size = prev.simn$ss.size, 
        ess.thresh = prev.simn$ess.thresh, resamp.method = prev.simn$resamp.method, 
        neg.wei.mech = prev.simn$neg.wei.mech, prev.simn = prev.simn, 
        progress.check = prev.simn$progress.check, phi.record = prev.simn$phi.record, 
        resamp.freq = prev.simn$resamp.freq, theta = prev.simn$theta, 
        p.path.renew = prev.simn$p.path.renew)
    append.len <- length(new.simn$anc.times) - 1
    if (append.len >= 1) {
        append.vec <- 2:length(new.simn$anc.times)
        p.anc.idx <- rbind(prev.simn$p.anc.idx, new.simn$p.anc.idx[append.vec, 
            , drop = FALSE])
        p.anc.wei <- rbind(prev.simn$p.anc.wei, new.simn$p.anc.wei[append.vec, 
            , drop = FALSE])
        p.anc.mat <- abind(prev.simn$p.anc.mat, new.simn$p.anc.mat[, 
            , append.vec, drop = FALSE], along = 3)
        p.anc.resamp <- c(prev.simn$p.anc.resamp, new.simn$p.anc.resamp)
        p.anc.ess <- c(prev.simn$p.anc.ess, new.simn$p.anc.ess)
        p.phi.hist <- rbind(prev.simn$p.phi.hist, new.simn$p.phi.hist)
        p.anc.neg <- rbind(prev.simn$p.anc.neg, new.simn$p.anc.neg)
        anc.times <- c(prev.simn$anc.times, new.simn$anc.times[append.vec])
    }
    else {
        p.anc.idx <- prev.simn$p.anc.idx
        p.anc.wei <- prev.simn$p.anc.wei
        p.anc.mat <- prev.simn$p.anc.mat
        p.anc.resamp <- prev.simn$p.anc.resamp
        p.anc.ess <- prev.simn$p.anc.ess
        p.phi.hist <- prev.simn$p.phi.hist
        p.anc.neg <- prev.simn$p.anc.neg
        anc.times <- prev.simn$anc.times
    }
    list(p.num = new.simn$p.num, p.idx = new.simn$p.idx, log.p.wei = new.simn$log.p.wei, 
        p.mat = new.simn$p.mat, p.layer = new.simn$p.layer, theta = new.simn$theta, 
        p.dapts = new.simn$p.dapts, p.anc.idx = p.anc.idx, p.anc.wei = p.anc.wei, 
        p.anc.mat = p.anc.mat, resamp.freq = new.simn$resamp.freq, 
        resamp.method = new.simn$resamp.method, neg.wei.mech = new.simn$neg.wei.mech, 
        p.anc.resamp = p.anc.resamp, p.anc.ess = p.anc.ess, p.anc.neg = p.anc.neg, 
        p.phi.hist = p.phi.hist, anc.times = anc.times, T.start = anc.times[1], 
        t.inc = new.simn$t.inc, T.fin = new.simn$T.fin, dimen = new.simn$dimen, 
        p.mu = new.simn$p.mu, ss.size = new.simn$ss.size, ess.thresh = new.simn$ess.thresh, 
        ss.phi = new.simn$ss.phi, ss.phiC = new.simn$ss.phiC, 
        p.cyc.arr = new.simn$p.cyc.arr, p.pass.arr = new.simn$p.pass.arr, 
        scale.transform = prev.simn$scale.transform, un.scale.transform = prev.simn$un.scale.transform, 
        progress.check = new.simn$progress.check, phi.record = new.simn$phi.record, 
        p.path.renew = new.simn$p.path.renew)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
