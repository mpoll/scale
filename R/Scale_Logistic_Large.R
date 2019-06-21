########################################################################
########################################################################
######                      Scale Algorithm                      #######
######                  Logistic Specification                   #######
######                      Large Data Example                   #######
######               Last Updated: 31/05/2019 MP                 #######
########################################################################
########################################################################

large.logistic.example    <- function(){

########################################################################
#### 0.1 - Save current seed and set data seed
########################################################################

curr.seed <- .Random.seed
set.seed(1)

########################################################################
########################################################################
### -1- Specification of Data Set
########################################################################
########################################################################

########################################################################
#### 1.1 - Size, dimension, and parameters
########################################################################

dsz             <<- 2^34                    # Size of data set
dimen           <<- 5                        # Dimensionality (>1)
beta.truth      <<- c(1,1,-1,2,-2)           # True parameters

########################################################################
#### 1.2 - Centering
########################################################################

## beta.star       <<- matrix(beta.truth,1,dimen) # What is the choosen centering value in the original space (DEFAULT / CLOSE TO SIMULATED DATA)
beta.star       <<- matrix(c(0.9999943, 0.9999501, -0.9999813, 1.9999872, -1.9999819), 1, dimen) # Fit using GLM and pooling
## n.sigma         <<- diag(rep(1,dimen))*dsz^(-1/2) # (DEFAULT / CLOSE TO SIMULATED DATA)
## n.sigma[1,1]    <- n.sigma[1,1]/2 # (DEFAULT / CLOSE TO SIMULATED DATA)
n.sigma         <<- diag(c(2.166573e-05, 3.806254e-05, 3.806251e-05, 4.222116e-05, 4.222066e-05)) # Fit using GLM and pooling
n.sigma.inv     <<- solve(n.sigma)  # Parameterise preconditioning matrix (inverse Fischers information) and specification of inverse preconditioning matrix

########################################################################
#### 1.3 - Data Recaller
########################################################################

library('msm')                            # Package to generate truncated Normals
design.thresh   <<- 1                     # Thresholding of truncated Normals
design.max <<- rep(design.thresh,dimen); design.min <<- c(1,rep(-1,design.thresh))

block.length <<- 2^8; total.blocks <<- dsz/block.length

datum     <<- function(datum.idx){ # Function to re-generate data (without storage) - won't work with seeds larger than 2^53-1
    old.seed <- .Random.seed; on.exit({ .Random.seed <<- old.seed }) # Record seed prior to execution
    datum.block.idx  <- (datum.idx-1)%%block.length + 1 # Compute data position in block
    seed.indexing <- (datum.idx-1)%/%block.length # Determine seed to be used
    set.seed(seed.indexing) # Set seed for covariate generation
    datum.x     <- cbind(1,matrix(rtnorm((block.length)*(dimen-1),mean=0,sd=1,lower=-design.thresh,upper=design.thresh),nrow=block.length,ncol=dimen-1)) # Generate final row of data matrix in block to block.idx
    datum.z     <- exp(datum.x %*% beta.truth) # Compute z
    datum.y     <- rbinom(block.length,1,datum.z/(1+datum.z)) # Generate data in block to block.idx
    list(idx=datum.idx,block.idx=seed.indexing+1,datum.block.idx=datum.block.idx,x=datum.x[datum.block.idx,],y=datum.y[datum.block.idx],seed=seed.indexing)} # Output list

data.seq     <<- function(block.idx){ # Function to re-generate data (without storage) - won't work with seeds larger than 2^53-1
    set.seed((block.idx-1))
    data.idxs   <- c(((block.idx-1)*block.length+1):(block.idx*block.length))
    datum.x     <- cbind(1,matrix(rtnorm((block.length)*(dimen-1),mean=0,sd=1,lower=-design.thresh,upper=design.thresh),nrow=block.length,ncol=dimen-1)) # Generate final row of data matrix in block to block.idx
    datum.z     <- exp(datum.x %*% beta.truth) # Compute z
    datum.y     <- rbinom(block.length,1,datum.z/(1+datum.z)) # Generate data in block to block.idx
    list(idxs=data.idxs,block.idx=block.idx,x=datum.x,y=datum.y,seed=(block.idx-1))}

########################################################################
#### 1.6 - Gradient and Laplacian Evaluation at beta.star
########################################################################

grad.log.pi <<- alpha.cent <<- matrix(c(-0.08079420,-0.04208156,0.04415533,-0.10456903,0.10873715),1,5) # Centering value of grad.log.pi at beta.star computed on cluster
alpha.cent.sq <<- (alpha.cent)%*%t(alpha.cent) # alpha^2 at centering
lap.log.pi <<- matrix(c(-1.208803, -1.064498, -1.064509, -1.215643, -1.215645),1,5) # Centering value of lap.log.pi at beta.star
alpha.p.cent <<- sum(lap.log.pi) # alpha prime at centering
phi.cent <<- (alpha.cent.sq+alpha.p.cent)/2 # phi at centering

########################################################################
########################################################################
### -2- Computation of Subsampling bounding functionals
########################################################################
########################################################################

data.extrema()

########################################################################
########################################################################
### -3- End Function Specification
########################################################################
########################################################################

.Random.seed <<- curr.seed

}
