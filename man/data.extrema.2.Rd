\name{data.extrema.2}
\alias{data.extrema.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Data Extema
  %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
data.extrema.2(design.min, design.max)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design.min}{
%%     ~~Describe \code{design.min} here~~
}
  \item{design.max}{
%%     ~~Describe \code{design.max} here~~
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
function (design.min, design.max) 
{
    extrema.design <- permutations(n = 2, r = dimen, v = c(0, 
        1), repeats.allowed = TRUE)
    for (i in 1:dimen) {
        extrema.design[, i] <- (1 - extrema.design[, i]) * design.min[i] + 
            extrema.design[, i] * design.max[i]
    }
    extrema.design <- rbind(extrema.design, extrema.design)
    data.extrema <- matrix(c(rep(0, dim(extrema.design)[1]/2), 
        rep(1, dim(extrema.design)[1]/2)), dim(extrema.design)[1], 
        1)
    data.precompute <- array(0, c(dim(extrema.design)[1], dimen, 
        dimen, 2))
    for (i in 1:dim(extrema.design)[1]) {
        for (dim.j in 1:dimen) {
            X.row <- extrema.design[i, ]
            data.precompute[i, 1:dimen, dim.j, 1] <- abs(-n.sigma[dim.j, 
                dim.j] * diag(n.sigma) * X.row[dim.j] * X.row)/4
            data.precompute[i, 1:dimen, dim.j, 2] <- abs(-(n.sigma[dim.j, 
                dim.j])^2 * diag(n.sigma) * (X.row[dim.j])^2 * 
                X.row) * 1/(6 * sqrt(3))
        }
    }
    dim.grad.max <- array(0, c(dimen, dimen, 2))
    for (dim.j in 1:dimen) {
        for (dim.k in 1:dimen) {
            dim.grad.max[dim.j, dim.k, 1] <- max(data.precompute[, 
                dim.j, dim.k, 1])
            dim.grad.max[dim.j, dim.k, 2] <- max(data.precompute[, 
                dim.j, dim.k, 2])
        }
    }
    dim.grad.max <<- dim.grad.max
    list(design.min = design.min, design.max = design.max, dim.grad.max = dim.grad.max)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
