\name{eazeta}
\alias{eazeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eazeta(n, s, t, x, y, L, U)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{U}{
%%     ~~Describe \code{U} here~~
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
function (n, s, t, x, y, L, U) 
{
    if (max(x - U, y - U, L - x, L - y) >= 0) {
        1
    }
    else {
        j <- 1:(ceiling(n/2))
        P <- -2/(t - s)
        D <- U - L
        D1 <- D * j + L
        D2 <- D * j - U
        z <- y - x
        if (even(n)) {
            sum(exp(P * (D1 - x) * (D1 - y)) + exp(P * (D2 + 
                x) * (D2 + y)) - exp(P * j^2 * D^2 - P * j * 
                D * z) - exp(P * j^2 * D^2 + P * j * D * z))
        }
        else {
            sum(exp(P * (D1 - x) * (D1 - y)) + exp(P * (D2 + 
                x) * (D2 + y))) - sum(exp(P * j[1:length(j) - 
                1]^2 * D^2 - P * j[1:length(j) - 1] * D * z) + 
                exp(P * j[1:length(j) - 1]^2 * D^2 + P * j[1:length(j) - 
                  1] * D * z))
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
