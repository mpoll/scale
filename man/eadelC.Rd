\name{eadelC}
\alias{eadelC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eadelC(mt, s, t, x, y, m, u)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mt}{
%%     ~~Describe \code{mt} here~~
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
  \item{m}{
%%     ~~Describe \code{m} here~~
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
function (mt, s, t, x, y, m, u) 
{
    if (mt >= mpfrthr) {
        pbn <- mt * mpfrpbn
        s <- mpfr(s, precBits = pbn)
        t <- mpfr(t, precBits = pbn)
        x <- mpfr(x, precBits = pbn)
        y <- mpfr(y, precBits = pbn)
        m <- mpfr(m, precBits = pbn)
        u <- mpfr(u, precBits = pbn)
        c(s1 = eadelR(mt, s, t, x, y, m, u), s2 = eadelR(mt + 
            1, s, t, x, y, m, u))
    }
    else {
        eadel_pair_cpp(mt, s, t, x, y, m, u)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
