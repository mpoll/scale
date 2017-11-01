\name{eabesex}
\alias{eabesex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Murray...%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eabesex(sV, tV, xV, yV, m, B1, B2, minI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sV}{
%%     ~~Describe \code{sV} here~~
}
  \item{tV}{
%%     ~~Describe \code{tV} here~~
}
  \item{xV}{
%%     ~~Describe \code{xV} here~~
}
  \item{yV}{
%%     ~~Describe \code{yV} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{B1}{
%%     ~~Describe \code{B1} here~~
}
  \item{B2}{
%%     ~~Describe \code{B2} here~~
}
  \item{minI}{
%%     ~~Describe \code{minI} here~~
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
function (sV, tV, xV, yV, m, B1, B2, minI) 
{
    if (minI == -1) {
        xV <- -xV
        yV <- -yV
        m <- -m
        B1 <- -B1
        B2 <- -B2
    }
    u <- runif(1, 0, 1)
    mt <- ceiling(max(sqrt(max(tV - sV) + (B1 - m)^2)/(2 * (B1 - 
        m)), sqrt(max(tV - sV) + (B2 - m)^2)/(2 * (B2 - m))))
    em <- matrix(0, length(sV), 8)
    em[, 1] <- sV
    em[, 2] <- tV
    em[, 3] <- xV
    em[, 4] <- yV
    B1evI <- B2evI <- 0
    while (B1evI == 0) {
        emM <- em
        for (i in 1:dim(em)[1]) {
            if (mt >= mpfrthr) {
                em[i, 5:6] <- as.numeric(eadelC(mt, em[i, 1], 
                  em[i, 2], em[i, 3], em[i, 4], m, B1))
            }
            else {
                em[i, 5:6] <- as.numeric(eadel_pair_cpp(mt, em[i, 
                  1], em[i, 2], em[i, 3], em[i, 4], m, B1))
            }
        }
        if (u <= prod(em[, 5])) {
            B1evI <- B2evI <- 1
            con1I <- 1
            con2I <- 1
            ex1I <- 0
            ex2I <- 0
        }
        else {
            if (u > prod(em[, 6])) {
                B1evI <- 1
                con1I <- 0
                ex1I <- 1
            }
            else {
                B1evI <- 0
                con1I <- 0
                ex1I <- 0
                mt <- mt + 2
            }
        }
    }
    while (B2evI == 0) {
        for (i in 1:dim(em)[1]) {
            if (mt >= mpfrthr) {
                em[i, 7:8] <- as.numeric(eadelC(mt, em[i, 1], 
                  em[i, 2], em[i, 3], em[i, 4], m, B2))
            }
            else {
                em[i, 7:8] <- as.numeric(eadel_pair_cpp(mt, em[i, 
                  1], em[i, 2], em[i, 3], em[i, 4], m, B2))
            }
        }
        if (u <= prod(em[, 7])) {
            B2evI <- 1
            con2I <- 1
            ex1I <- 0
        }
        else {
            if (u > prod(em[, 8])) {
                B2evI <- 1
                con2I <- 0
                ex2I <- 1
            }
            else {
                B2evI <- 0
                con2I <- 0
                ex2I <- 0
                mt <- mt + 2
            }
        }
    }
    if (minI == -1) {
        em[, 3] <- -em[, 3]
        em[, 4] <- -em[, 4]
    }
    accI <- 0
    if (con1I == 1) {
        accI <- 1
    }
    else {
        if (con2I == 1) {
            if (rbinom(1, 1, 0.5) == 1) {
                accI <- 1
            }
        }
    }
    list(accI = accI, u = u, con1I = con1I, con2I = con2I, em = em)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
