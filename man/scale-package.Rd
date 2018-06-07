\name{scale-package}
\alias{scale-package}
\alias{scale}
\docType{package}
\title{
\packageTitle{scale}
}
\description{
\packageDescription{scale}
}
\details{
scale provides a reference implementation of the Scaleable Langevin
Exact Algorithm as described in Pollock et al. (2016).

Note that this implementation attempts to explicitly replicate the
mathematical formulation of the algorithm and computational efficiency
has been sacrificed to a considerable degree in order to maintain the
relationship between the mathematical presentation and the
implementation itself.
}
\section{Examples}{
This section explains how the examples presented in Pollock et
al. (2016) can be reproduced.

\itemize{
 \item Download the relevant data files from
 \url{http://www.warwick.ac.uk/mpollock/scale}
 
 \item Download the precomputed centering and preconditioning data from
 \url{http://www.warwick.ac.uk/mpollock/scale}
}
 \strong{Logistic Regression With Limited Information}
First, load this library, and download \code{uniformative.Rdata} from 
 \url{http://www.warwick.ac.uk/mpollock/scale} (this file contains 
  the data, the precomputed centering value and
  preconditioning matrix, and some helper functions.): \preformatted{    library("scale")
    load("uninformative.RData")}
  Now run the algorithm: \preformatted{ scale_logistic(fnm="uninformative-output.Rdata",p.num=1024,t.inc=0.0001,T.extend=0.0001,run.length=10,ss.size=2,ss.on=TRUE,seed.default=1,data.precompute=TRUE)}    
and after some time, the file \code{uninformative-output.Rdata} will contain
the simulation output as an object with name \code{simn}. This can be
  supplied to \code{scale_ergodic} in order to obtain quantities
  suitable for computing ergodic averages as described in Pollock et
  al. (2016).

 \strong{Logistic Regression for Airline Data}
  First, load this library, and download \code{airline.Rdata} from 
 \url{http://www.warwick.ac.uk/mpollock/scale} (this file contains
  the airline data, the precomputed centering value and
  preconditioning matrix and some helper functions.): \preformatted{    library("scale")
    load("airline.RData")}
  Now run the algorithm: \preformatted{ scale_logistic(fnm="airline-output.Rdata",p.num=1024,t.inc=0.00001,T.extend=0.00001,run.length=10,ss.size=2,ss.on=TRUE,seed.default=1,data.precompute=TRUE)}    
and after some time, the file \code{airline-output.Rdata} will contain
the simulation output as an object with name \code{simn}. This can be
  supplied to \code{scale_ergodic} in order to obtain quantities
  suitable for computing ergodic averages as described in Pollock et
  al. (2016).

 \strong{Large Logistic Regression}
 First, load this library, and download \code{large.Rdata} from 
 \url{http://www.warwick.ac.uk/mpollock/scale} (this file contains
  the precomputed centering value and preconditioning matrix and some
  )helper functions (including one to provide individual data
  values as required by the algorithm.): \preformatted{    library("scale")
    load("large.RData")}
  Now run the algorithm: \preformatted{ scale_logistic(fnm="large-output.Rdata",p.num=4096,t.inc=0.00001,T.extend=0.00001,run.length=10,ss.size=64,ss.on=TRUE,seed.default=1,data.precompute=TRUE)}    
and after some time, the file \code{large-output.Rdata} will contain
the simulation output as an object with name \code{simn}. This can be
  supplied to \code{scale_ergodic} in order to obtain quantities
  suitable for computing ergodic averages as described in Pollock et
  al. (2016).
}
\author{
\packageAuthor{scale}

Maintainer: \packageMaintainer{scale}
}
\references{
  Pollock, Murray, Paul Fearnhead, Adam M. Johansen, and Gareth O. Roberts. "The scalable Langevin exact algorithm: Bayesian inference for big data." arXiv preprint arXiv:1609.03436 (2016).

}
\keyword{ package }
\seealso{ scale_logistic, scale_ergodic
}
\examples{
}
