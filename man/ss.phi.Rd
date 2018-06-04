\name{ss.phi}
\alias{ss.phi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SS.phi
  %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ss.phi(eta.pos, dapts, ss.size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta.pos}{
%%     ~~Describe \code{eta.pos} here~~
}
  \item{dapts}{
%%     ~~Describe \code{dapts} here~~
}
  \item{ss.size}{
%%     ~~Describe \code{ss.size} here~~
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
function (eta.pos, dapts, ss.size) 
{
    beta.pos <- un.scale.transform(eta.pos)
    factor <- dsz/ss.size
    data1 <- data.eval(beta.pos, dapts[1, , "dapt.1"], factor = factor)
    data1.star <- data.eval(beta.star, dapts[1, , "dapt.1"], 
        factor = factor)
    data2 <- data.eval(beta.pos, dapts[1, , "dapt.2"], factor = factor)
    data2.star <- data.eval(beta.star, dapts[1, , "dapt.2"], 
        factor = factor)
    data3 <- data.eval(beta.pos, dapts[1, , "dapt.3"], factor = factor)
    data3.star <- data.eval(beta.star, dapts[1, , "dapt.3"], 
        factor = factor)
    ((data1$grad.log.pi - data1.star$grad.log.pi) \%*\% t(2 * alpha.cent + 
        data2$grad.log.pi - data2.star$grad.log.pi) + sum(data3$lap.log.pi) - 
        sum(data3.star$lap.log.pi) + alpha.cent.sq + alpha.p.cent)/2
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
