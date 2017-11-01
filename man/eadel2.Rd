\name{eadel2}
\alias{eadel2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eadel2(n, s, t, m, xoy, u)
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
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{xoy}{
%%     ~~Describe \code{xoy} here~~
}
  \item{u}{
%%     ~~Describe \code{u} here~~
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
function (n, s, t, m, xoy, u) 
{
    if (max(xoy - u, m - xoy) >= 0) {
        0
    }
    else {
        if (even(n)) {
            j <- 1:(n/2)
            1 - (sum(eapsi(j, s, t, m, xoy, u) - eachi(j, s, 
                t, m, xoy, u)))/(xoy - m)
        }
        else {
            if (n > 1) {
                j <- 1:((n - 1)/2)
                1 - (sum(eapsi(j, s, t, m, xoy, u) - eachi(j, 
                  s, t, m, xoy, u)))/(xoy - m) - eapsi(max(j) + 
                  1, s, t, m, xoy, u)/(xoy - m)
            }
            else {
                1 - eapsi(1, s, t, m, xoy, u)/(xoy - m)
            }
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
