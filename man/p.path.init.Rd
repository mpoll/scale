\name{p.path.init}
\alias{p.path.init}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
p.path.init(s, x, theta, dimen, Jst.t = t.opt, Jst.rat = rat.opt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{dimen}{
%%     ~~Describe \code{dimen} here~~
}
  \item{Jst.t}{
%%     ~~Describe \code{Jst.t} here~~
}
  \item{Jst.rat}{
%%     ~~Describe \code{Jst.rat} here~~
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
function (s, x, theta, dimen, Jst.t = t.opt, Jst.rat = rat.opt) 
{
    passage <- matrix(0, dimen, 6)
    colnames(passage) <- c("dimen", "tau", "y", "minI", "l", 
        "u")
    for (i in 1:dimen) {
        ith.passage <- bm.pass(s, x[i], theta[i], Jst.t, Jst.rat)
        passage[i, ] <- c(i, ith.passage$tau, ith.passage$y, 
            ith.passage$minI, x[i] - theta[i], x[i] + theta[i])
    }
    list(curr.time = s, pass.mat = passage, cyc.mat = mat.sort.c(passage, 
        dimen, "tau"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
