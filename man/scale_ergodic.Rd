\name{scale_ergodic}
\alias{scale_ergodic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scale_ergodic(simn, retain.frac = NULL, retain.frac.range = NULL, trunc.time = NULL, trunc.time.range = NULL, symm = NULL, even.weights = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simn}{
%%     ~~Describe \code{simn} here~~
}
  \item{retain.frac}{
%%     ~~Describe \code{retain.frac} here~~
}
  \item{retain.frac.range}{
%%     ~~Describe \code{retain.frac.range} here~~
}
  \item{trunc.time}{
%%     ~~Describe \code{trunc.time} here~~
}
  \item{trunc.time.range}{
%%     ~~Describe \code{trunc.time.range} here~~
}
  \item{symm}{
%%     ~~Describe \code{symm} here~~
}
  \item{even.weights}{
%%     ~~Describe \code{even.weights} here~~
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
function (simn, retain.frac = NULL, retain.frac.range = NULL, 
    trunc.time = NULL, trunc.time.range = NULL, symm = NULL, 
    even.weights = NULL) 
{
    if (is.null(retain.frac) == TRUE) {
        frac.idxs <- 1:length(simn$anc.times)
    }
    else {
        frac.idxs <- floor(length(simn$anc.times) * retain.frac):length(simn$anc.times)
    }
    if (is.null(retain.frac.range) == FALSE) {
        frac.idxs <- floor(length(simn$anc.times) * retain.frac.range[1]):ceiling(length(simn$anc.times) * 
            retain.frac.range[2])
    }
    if (is.null(trunc.time) == TRUE) {
        time.idxs <- 1:length(simn$anc.times)
    }
    else {
        if (simn$T.fin < trunc.time) {
            time.idxs <- numeric(0)
        }
        else {
            time.idxs <- length(simn$anc.times[simn$anc.times <= 
                trunc.time]):length(simn$anc.times)
        }
    }
    if (is.null(trunc.time.range) == FALSE) {
        time.idxs <- length(simn$anc.times[simn$anc.times <= 
            trunc.time.range[1]]):length(simn$anc.times[simn$anc.times <= 
            trunc.time.range[2]])
    }
    idxs <- intersect(frac.idxs, time.idxs)
    p.num <- simn$p.num
    l.idxs <- length(idxs)
    c.idxs <- 1:(p.num * l.idxs)
    r.idxs <- numeric(0)
    for (i in 1:l.idxs) {
        r.idxs <- c(r.idxs, (0:(p.num - 1)) * l.idxs + i)
    }
    times <- simn$anc.times[idxs]
    if (length(idxs) == 0) {
        weights <- numeric(0)
        draws <- draws <- matrix(0, 0, simn$dimen, dimnames = list(NULL, 
            sprintf("dim.\%i", 1:simn$dimen)))
    }
    else {
        simn.weis <- simn$p.anc.wei[idxs, , drop = FALSE]
        weights <- simn$p.anc.wei[r.idxs]
        weights <- weights/sum(weights)
        simn.draws <- simn$p.anc.mat[, , idxs, drop = FALSE]
        draws <- matrix(0, simn$p.num * length(idxs), simn$dimen, 
            dimnames = list(NULL, sprintf("dim.\%i", 1:simn$dimen)))
        for (i in 1:simn$dimen) {
            simn.draws.dim <- simn.draws[i, , , drop = FALSE]
            draws[, i] <- simn.draws.dim[r.idxs]
        }
        if (is.null(symm) == FALSE) {
            draws <- rbind(draws, -draws)
            weights <- c(weights/2, weights/2)
        }
    }
    if (is.null(even.weights) == FALSE) {
        if (even.weights == TRUE) {
            p.idxs <- sample(length(weights), length(weights), 
                replace = TRUE, prob = weights)
            weights <- rep(1/length(weights), length(weights))
            draws <- draws[p.idxs, , drop = FALSE]
        }
    }
    list(idxs = idxs, times = times, weights = weights, draws = draws)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
