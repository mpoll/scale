\name{ins.sort.c}
\alias{ins.sort.c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ins.sort.c(vec, sor.mat, sor.mat.sz, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
%%     ~~Describe \code{vec} here~~
}
  \item{sor.mat}{
%%     ~~Describe \code{sor.mat} here~~
}
  \item{sor.mat.sz}{
%%     ~~Describe \code{sor.mat.sz} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (vec, sor.mat, sor.mat.sz, n) 
{
    r.min <- 1
    r.max <- sor.mat.sz
    tau <- vec[n]
    if (tau > sor.mat[r.max, n]) {
        r.min <- r.max <- sor.mat.sz + 1
    }
    else {
        while (r.max > r.min) {
            r.mid <- floor(0.5 * (r.max + r.min))
            if (tau <= sor.mat[r.mid, n]) {
                r.max <- r.mid
            }
            else {
                r.min <- r.mid + 1
            }
        }
    }
    return(rbind(sor.mat[seq_len(r.min - 1), , drop = FALSE], 
        vec, sor.mat[seq_len(sor.mat.sz - r.max + 1) + r.max - 
            1, , drop = FALSE]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
