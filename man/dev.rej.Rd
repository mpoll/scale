\name{dev.rej}
\alias{dev.rej}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dev.rej(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
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
function (X) 
{
    S <- pi * exp(-pi^2 * X/8)/2
    n <- 0
    Y <- runif(1, 0, 1) * S
    repeat {
        n <- n + 1
        S <- S - (-1)^n * pi * (n + 1/2) * exp(-(n + 1/2)^2 * 
            pi^2 * X/2)
        if (Y < S) {
            Acc <- 1
            break
        }
        n <- n + 1
        S <- S + (-1)^n * pi * (n + 1/2) * exp(-(n + 1/2)^2 * 
            pi^2 * X/2)
        if (Y > S) {
            Acc <- 0
            break
        }
    }
    list(Acc = Acc, S = S, n = n, Y = Y, X = X)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
