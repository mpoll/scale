\name{p.path.renew}
\alias{p.path.renew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
p.path.renew(pass.mat, cyc.mat, curr.time, x.curr, next.time, theta, dimen)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pass.mat}{
%%     ~~Describe \code{pass.mat} here~~
}
  \item{cyc.mat}{
%%     ~~Describe \code{cyc.mat} here~~
}
  \item{curr.time}{
%%     ~~Describe \code{curr.time} here~~
}
  \item{x.curr}{
%%     ~~Describe \code{x.curr} here~~
}
  \item{next.time}{
%%     ~~Describe \code{next.time} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{dimen}{
%%     ~~Describe \code{dimen} here~~
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
function (pass.mat, cyc.mat, curr.time, x.curr, next.time, theta, 
    dimen) 
{
    if (dimen > 1) {
        x.next <- numeric(dimen)
        for (i in 1:dimen) {
            repeat {
                x.next[i] <- eabesmid(next.time, curr.time, pass.mat[i, 
                  "tau"], pass.mat[i, "tau"], x.curr[i], pass.mat[i, 
                  "y"], pass.mat[i, "y"], pass.mat[i, "minI"])$w
                if (eabesex(sV = c(curr.time, next.time), tV = c(next.time, 
                  pass.mat[i, "tau"]), xV = c(x.curr[i], x.next[i]), 
                  yV = c(x.next[i], pass.mat[i, "y"]), m = pass.mat[i, 
                    "y"], B1 = x.curr[i] + pass.mat[i, "minI"] * 
                    theta[i], B2 = x.curr[i] + pass.mat[i, "minI"] * 
                    theta[i], minI = pass.mat[i, "minI"])$accI == 
                  1) {
                  break
                }
            }
        }
        p.taus <- p.path.init(s = next.time, x = x.next, theta = theta, 
            dimen = dimen)
        pass.mat <- p.taus$pass.mat
        cyc.mat <- p.taus$cyc.mat
        list(curr.time = next.time, next.tau = cyc.mat[1, "tau"], 
            max.next = cyc.mat[1, "tau"] - next.time, pass.mat = pass.mat, 
            cyc.mat = cyc.mat, x.next = x.next)
    }
    else {
        repeat {
            x.next <- eabesmid(next.time, curr.time, pass.mat["tau"], 
                pass.mat["tau"], x.curr, pass.mat["y"], pass.mat["y"], 
                pass.mat["minI"])$w
            if (eabesex(sV = c(curr.time, next.time), tV = c(next.time, 
                pass.mat["tau"]), xV = c(x.curr, x.next), yV = c(x.next, 
                pass.mat["y"]), m = pass.mat["y"], B1 = x.curr + 
                pass.mat["minI"] * theta, B2 = x.curr + pass.mat["minI"] * 
                theta, minI = pass.mat["minI"])$accI == 1) {
                break
            }
        }
        new.fpt <- bm.pass(s = next.time, x = x.next, theta = theta, 
            Jst.t = t.opt, Jst.rat = rat.opt)
        pass.mat[] <- c(update.dimen, new.fpt$tau, new.fpt$y, 
            new.fpt$minI, x.next - theta, x.next + theta)
        list(curr.time = next.time, next.tau = pass.mat["tau"], 
            max.next = pass.mat["tau"] - next.time, pass.mat = pass.mat, 
            cyc.mat = pass.mat, x.next = x.next)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
