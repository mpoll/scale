############################################
############################################
######       Scale Algorithm         #######
######    Logistic Specification     #######
######         Exact Version         #######
######  Last Updated: 26/07/2017 MP  #######
############################################
############################################

############################################
############################################
### -1- Specification of Data Set
############################################
############################################

#if(exists("taskid")==FALSE){taskid <- 1}
#set.seed(taskid) # Set seed [**Required for datum generator**]

#############################################
#### 1.1 - Size, dimension, and parameters
#############################################

dsz             <- 10^10                    # Size of data set
dimen           <- 5                        # Dimensionality (>1)
beta.truth      <- c(1,1,-1,2,-2)   # True parameters
n.sigma         <- diag(rep(1,dimen))*dsz^(-1/2); n.sigma[1,1] <- n.sigma[1,1]/2; n.sigma.inv <- solve(n.sigma)  # Parameterise preconditioning matrix (INV FISCH INF) and specification of inverse preconditioning matrix

#############################################
#### 1.2 - Centering
#############################################

beta.star       <- matrix(beta.truth,1,dimen) # What is the choosen centering value in the original space

#############################################
#### 1.3 - Transformation specification
#############################################

transform       <- function(beta.pos){matrix(c(diag(n.sigma.inv)*(beta.pos-beta.star)),nrow=1)} # Transformation from beta to eta space
un.transform    <- function(eta.pos){matrix(c(diag(n.sigma)*eta.pos+beta.star),nrow=1)} # Transformation from eta to beta space

eta.truth       <- transform(beta.truth); eta.star <- transform(beta.star) # Specification of truth  and centering point

#############################################
#### 1.3 - Data Recaller
#############################################

library('msm')                            # Package to generate truncated Normals
design.thresh   <- 1                    # Thresholding of truncated Normals

datum     <- function(datum.idx){ # Function to re-generate data (without storage) - won't work with seeds larger than 2^53-1
    old.seed <- .Random.seed; on.exit({ .Random.seed <<- old.seed })
    seed.indexing <- datum.idx%/%2147483647
    datum.seeding <- datum.idx%%2147483647
    set.seed(datum.seeding); runif((2*dimen-1)*seed.indexing) # junk indexing by the number of random operations used to simulate data.gen
    datum.x     <- c(1,rtnorm(dimen-1,mean=0,sd=1,lower=-design.thresh,upper=design.thresh))
    datum.z     <- exp(sum(beta.truth*datum.x))
    datum.y     <- rbinom(1,1,datum.z/(1+datum.z))
    list(idx=datum.idx,x=datum.x,y=datum.y,seed=datum.seeding,seed.idx=seed.indexing)}

#############################################
#### 1.4 - Specification of derivatives (1st and 2nd)
#############################################

datum.eval   <- function(beta.pos,data.idx){
    #### CALL DATA
    datum.call      <- datum(data.idx)
    #### EVALUATE
    datum.z         <- exp(sum(datum.call$x*beta.pos))
    log.pi          <- datum.call$x*beta.pos*datum.call$y-log(1+datum.z)
    grad.log.pi     <- diag(n.sigma)*datum.call$x*(datum.call$y-datum.z/(1+datum.z))
    lap.log.pi      <- -(diag(n.sigma)*datum.call$x)^2*datum.z/(1+datum.z)^2
    list(log.pi=log.pi,grad.log.pi=grad.log.pi,lap.log.pi=lap.log.pi)}

#############################################
#### 1.5 - Function for multiple random evaluations of grad and lap log pi
#############################################

data.eval       <- function(beta.pos,data.idxs,factor=1){ # Factor is a multiple of output grad log and lap log
    grad.log.pi <- lap.log.pi <- numeric(dimen)
    for(i in 1:length(data.idxs)){
        datum.call          <- datum.eval(beta.pos,data.idxs[i])
        grad.log.pi         <- grad.log.pi + datum.call$grad.log.pi
        lap.log.pi          <- lap.log.pi + datum.call$lap.log.pi}
    list(grad.log.pi=factor*grad.log.pi,lap.log.pi=factor*lap.log.pi)}

#############################################
#### 1.6 - Function for multiple determined evaluations at centering point (ie in initialisation)
#############################################

data.init       <- function(datum.idx.start,length.idxs){data.eval(beta.star,c(datum.idx.start:(datum.idx.start+length.idxs-1)))}

############################################
############################################
### -2- Specification of Intensities
############################################
############################################

#############################################
#### 2.1 - Computation of Subsampling bounding functionals
#############################################

#### 2.1.1 - Precompute the most extreme data values

library(gtools)
extrema.design  <- rbind(cbind(1,permutations(n=2,r=dimen-1,v=c(-design.thresh,design.thresh),repeats.allowed=TRUE)),cbind(1,permutations(n=2,r=dimen-1,v=c(-design.thresh,design.thresh),repeats.allowed=TRUE)))
extrema.size <- dim(extrema.design)[1]
data.extrema <- matrix(c(rep(0,extrema.size/2),rep(1,extrema.size/2)),extrema.size,1)

#### 2.1.2 - Precompute the maxima Gradient Datum Log likelihood and Gradient Laplacian Datum Log likelihood for each datum

data.precompute <- array(0,c(extrema.size,dimen,dimen,2)); for(i in 1:extrema.size){for(dim.j in 1:dimen){
    X.row <- extrema.design[i,] # Determine row of design matrix
    data.precompute[i,1:dimen,dim.j,1] <- abs(-n.sigma[dim.j,dim.j]*diag(n.sigma)*X.row[dim.j]*X.row)/4 # maxima Gradient Datum Log likelihood
    data.precompute[i,1:dimen,dim.j,2] <- abs(-(n.sigma[dim.j,dim.j])^2*diag(n.sigma)*(X.row[dim.j])^2*X.row)*1/(6*sqrt(3))}} # Maxima Gradient Laplacian Datum Log likelihood

#### 2.1.3 - Precompute the maxima Gradient Datum Log likelihood and Laplacian Datum Log likelihood over all data

dim.grad.max  <- array(0,c(dimen,dimen,2)); for(dim.j in 1:dimen){for(dim.k in 1:dimen){
    dim.grad.max[dim.j,dim.k,1] <- max(data.precompute[,dim.j,dim.k,1])     # Maximal Gradient Log Likelihood
    dim.grad.max[dim.j,dim.k,2] <- max(data.precompute[,dim.j,dim.k,2])}}   # Maximal Laplacian Log Likelhood

#############################################
#### 2.2 - Subsampling Phi Intensity
#############################################
ss.phiC <- function(p.mat.curr,l.bound=NULL,u.bound=NULL){ # Function of current particle location, lower bounds (dimensionally), upper bounds (dimensionally)
    l.bds <- pmin(0,l.bound); u.bds <- pmax(0,u.bound) # Extend the hyper cube to the origin (i.e. eta.star)
    distance <- pmax(abs(l.bds),abs(u.bds))            # Compute the maximal distance in each dimension
    alpha.bds <- alpha.prime.bds <- numeric(dimen); for(j in 1:dimen){ # Compute the bounds for alpha and alpha prime
        alpha.bds[j]        <- dsz*sum(distance*dim.grad.max[j,,1]) # alpha bounds
        alpha.prime.bds[j]  <- dsz*sum(distance*dim.grad.max[j,,2]) # alpha prime bounds
    }
    A.bds   <- alpha.bds^2 # Compute bounds on alpha^2
    phi.bds <- (sum(2*abs(alpha.cent)*alpha.bds)+sum(A.bds)+sum(alpha.prime.bds))/2 # Compute bounds (modulo normalisation) on subsampled phi
    phiL <- phi.cent - phi.bds; phiU <- phi.cent + phi.bds; intensity <- phiU - phiL # Compute Bounding functionals and intensity
    list(distance=distance,alpha.bds=alpha.bds,alpha.prime.bds=alpha.prime.bds,A.bds=A.bds,phiL=phiL,phiU=phiU,intensity=intensity)}
