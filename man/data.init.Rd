\name{data.init}
\alias{data.init}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data initialization%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
data.init(beta.eval = beta.star, datum.idx.start = 1, length.idxs = dsz)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta.eval}{
%%     ~~Describe \code{beta.eval} here~~
}
  \item{datum.idx.start}{
%%     ~~Describe \code{datum.idx.start} here~~
}
  \item{length.idxs}{
%%     ~~Describe \code{length.idxs} here~~
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
function (beta.eval = beta.star, datum.idx.start = 1, length.idxs = dsz) 
{
    evaluate <- data.eval(beta.eval, c(datum.idx.start:(datum.idx.start + 
        length.idxs - 1)))
    alpha.cent <<- matrix(as.numeric(evaluate$grad.log.pi), nrow = 1)
    alpha.cent.sq <<- (alpha.cent) \%*\% t(alpha.cent)
    alpha.p.cent <<- sum(as.numeric(evaluate$lap.log.pi))
    phi.cent <<- (alpha.cent.sq + alpha.p.cent)/2
    list(alpha.cent = alpha.cent, alpha.cent.sq = alpha.cent.sq, 
        alpha.p.cent = alpha.p.cent, phi.cent = phi.cent)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
