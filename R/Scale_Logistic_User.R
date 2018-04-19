############################################
############################################
######       Scale Algorithm         #######
######    Logistic Specification     #######
######         User Example          #######
######  Last Updated: 18/04/2018 MP  #######
############################################
############################################

user.logistic.example   <- function(){

############################################
############################################
### -1- Specification of Data Set
############################################
############################################

#############################################
#### 1.1 - Size, dimension, and parameters
#############################################

examp.design <- matrix(0,40,5); examp.design[,2:5] <- rnorm(160)
examp.data   <- rbinom(160,1,0.5)
beta.star    <- c(1,1,-1,2,-2) # What is the choosen centering value in the original space

if(exists("dsz")==FALSE){dsz <<- dim(examp.design)[1]}     # Size of data set
if(exists("dimen")==FALSE){dimen <<- dim(examp.design)[2]} # Dimensionality (>1)
if(exists("n.sigma")==FALSE){n.sigma <<- diag(rep(1,dimen))*dsz^(-1/2)} # Parameterise preconditioning matrix (inverse Fischers information)
if(exists("n.sigma.inv")==FALSE){n.sigma.inv <- solve(n.sigma)} # Specification of inverse preconditioning matrix





#############################################
#### 1.3 - Transformation specification
#############################################

scale.transform       <- function(beta.pos){matrix(c(diag(n.sigma.inv)*(beta.pos-beta.star)),nrow=1)} # Transformation from beta to eta space
un.scale.transform    <- function(eta.pos){matrix(c(diag(n.sigma)*eta.pos+beta.star),nrow=1)} # Transformation from eta to beta space

#############################################
#### 1.3 - Data Recaller
#############################################

datum     <- function(datum.idx){
    datum.x     <- examp.design[datum.idx,]
    datum.y     <- examp.data[datum.idx,]
    list(idx=datum.idx,x=datum.x,y=datum.y)}

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

datum.max.loglike <- function(design.row,datum)
datum.max.laplike <- function(design.row,datum)


data.extrema         <- function(){
    #### .1 - Precompute the most extreme data values
    if(exists("design.thresh")==TRUE){if(length(design.thresh)==1){design.min <- c(1,rep(-design.thresh,dimen-1)); design.max <- c(1,rep(design.thresh,dimen-1))}}
    if(exists("design.thresh")==FALSE){design.min <- apply(examp.design,2,min); design.max <- apply(examp.design,2,max)}
    extrema.design <- permutations(n=2,r=dimen,v=c(0,1),repeats.allowed=TRUE); for(i in 1:dimen){extrema.design[,i] <- (1-extrema.design[,i])*design.min[i] + extrema.design[,i]*design.max[i]}; extrema.design <- rbind(extrema.design,extrema.design)
    extrema.size <- dim(extrema.design)[1]
    data.extrema <- matrix(c(rep(0,extrema.size/2),rep(1,extrema.size/2)),extrema.size,1)
    
    #### .2 - Precompute the maxima Gradient Datum Log likelihood and Gradient Laplacian Datum Log likelihood for each datum
    data.precompute <- array(0,c(extrema.size,dimen,dimen,2)); for(i in 1:extrema.size){for(dim.j in 1:dimen){
        X.row <- extrema.design[i,] # Determine row of design matrix
        data.precompute[i,1:dimen,dim.j,1] <- abs(-n.sigma[dim.j,dim.j]*diag(n.sigma)*X.row[dim.j]*X.row)/4 # maxima Gradient Datum Log likelihood
        data.precompute[i,1:dimen,dim.j,2] <- abs(-(n.sigma[dim.j,dim.j])^2*diag(n.sigma)*(X.row[dim.j])^2*X.row)*1/(6*sqrt(3))}} # Maxima Gradient Laplacian Datum Log likelihood
    
    #### .3 - Precompute the maxima Gradient Datum Log likelihood and Laplacian Datum Log likelihood over all data
    dim.grad.max  <- array(0,c(dimen,dimen,2)); for(dim.j in 1:dimen){for(dim.k in 1:dimen){
        dim.grad.max[dim.j,dim.k,1] <- max(data.precompute[,dim.j,dim.k,1])     # Maximal Gradient Log Likelihood
        dim.grad.max[dim.j,dim.k,2] <- max(data.precompute[,dim.j,dim.k,2])}}   # Maximal Laplacian Log Likelhood
    
    #### .4 - Write to global scope dim.grad.max
    dim.grad.max <<- dim.grad.max
}



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
    phi.bds <- (sum(2*abs(scale::alpha.cent)*alpha.bds)+sum(A.bds)+sum(alpha.prime.bds))/2 # Compute bounds (modulo normalisation) on subsampled phi
    phiL <- scale::phi.cent - phi.bds; phiU <- scale::phi.cent + phi.bds; intensity <- phiU - phiL # Compute Bounding functionals and intensity
    list(distance=distance,alpha.bds=alpha.bds,alpha.prime.bds=alpha.prime.bds,A.bds=A.bds,phiL=phiL,phiU=phiU,intensity=intensity)}

############################################
############################################
### -3- End Function Specification
############################################
############################################

}

