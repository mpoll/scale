\name{dev.pr}
\alias{dev.pr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dev.pr(Jst.t = t.opt, Jst.rat = rat.opt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (Jst.t = t.opt, Jst.rat = rat.opt) 
{
    U <- runif(1, 0, 1)
    if (U < Jst.rat) {
        E <- rexp(1, 1)
        X <- Jst.t + 8 * E/(pi^2)
    }
    else {
        repeat {
            E <- rexp(2, 1)
            if (E[1]^2 <= 2 * E[2]/Jst.t) {
                X <- Jst.t/(1 + Jst.t * E[1])^2
                break
            }
        }
    }
    list(X = X, U = U, E = E)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
