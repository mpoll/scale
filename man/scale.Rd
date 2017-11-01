\name{scale}
\alias{scale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scale(p.num, t.inc, T.fin, ss.phi, ss.phiC, dimen, transform, un.transform, T.start = 0, x.init = NULL, ss.size = 1, ess.thresh = 0.5, resamp.method = resid.resamp, neg.wei.mech = scale_zero.wei, prev.simn = NULL, progress.check = FALSE, phi.record = FALSE, resamp.freq = p.num - 1, theta = NULL, p.path.renew = p.path.renew, scale_version = scale_approx)
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
  \item{transform}{
%%     ~~Describe \code{transform} here~~
}
  \item{un.transform}{
%%     ~~Describe \code{un.transform} here~~
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
  \item{scale_version}{
%%     ~~Describe \code{scale_version} here~~
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
function (p.num, t.inc, T.fin, ss.phi, ss.phiC, dimen, transform, 
    un.transform, T.start = 0, x.init = NULL, ss.size = 1, ess.thresh = 0.5, 
    resamp.method = resid.resamp, neg.wei.mech = scale_zero.wei, 
    prev.simn = NULL, progress.check = FALSE, phi.record = FALSE, 
    resamp.freq = p.num - 1, theta = NULL, p.path.renew = p.path.renew, 
    scale_version = scale_approx) 
{
    scale_version(p.num, t.inc, T.fin, ss.phi, ss.phiC, dimen, 
        transform, un.transform, T.start = T.start, x.init = x.init, 
        ss.size = ss.size, ess.thresh = ess.thresh, resamp.method = resamp.method, 
        neg.wei.mech = neg.wei.mech, prev.simn = prev.simn, progress.check = progress.check, 
        phi.record = phi.record, resamp.freq = resamp.freq, theta = theta, 
        p.path.renew = p.path.renew)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
