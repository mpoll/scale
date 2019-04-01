########################################################################
########################################################################
######                      Scale Algorithm                      #######
######                  Logistic Specification                   #######
######                      Large Data Example                   #######
######               Last Updated: 06/02/2019 MP                 #######
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
n.sigma         <<- diag(rep(1,dimen))*dsz^(-1/2)
n.sigma[1,1]    <- n.sigma[1,1]/2
n.sigma.inv     <<- solve(n.sigma)  # Parameterise preconditioning matrix (inverse Fischers information) and specification of inverse preconditioning matrix

########################################################################
#### 1.2 - Centering
########################################################################

beta.star       <<- matrix(beta.truth,1,dimen) # What is the choosen centering value in the original space

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

##### FUNCTIONS WHICH INDEX RANDOMNESS (suffixed .old) -- DO NOT WORK DUE TO RANDOM NUMBER OF RANDOM NUMBERS REQUIRED IN TRUNC NORM SIMULATION

datum.old     <<- function(datum.idx){ # Function to re-generate data (without storage) - won't work with seeds larger than 2^53-1
    old.seed <- .Random.seed; on.exit({ .Random.seed <<- old.seed })
    datum.seeding  <- (datum.idx-1)%%block.length +1
    seed.indexing <- (datum.idx-1)%/%block.length
    set.seed(seed.indexing); runif((2*dimen+1)*datum.seeding) # junk indexing by the number of random operations used to simulate data.gen
    datum.x     <- c(1,rtnorm(dimen-1,mean=0,sd=1,lower=-design.thresh,upper=design.thresh))
    datum.z     <- exp(sum(beta.truth*datum.x))
    datum.y     <- rbinom(1,1,datum.z/(1+datum.z))
    list(idx=datum.idx,x=datum.x,y=datum.y,seed=seed.indexing,seed.idx=datum.seeding)}

data.seq.old     <<- function(block.idx){ # Function to re-generate data (without storage) - won't work with seeds larger than 2^53-1
    data.idxs   <- c(((block.idx-1)*block.length+1):(block.idx*block.length))
    x <- matrix(0,nrow=0,ncol=dimen); y <- numeric(0)
    for(i in 1:block.length){generate.data <- datum.old(data.idxs[i]); x <- rbind(x,generate.data$x); y <- c(y,generate.data$y)}
    list(idxs=data.idxs,block.idx=block.idx,x=x,y=y)}

########################################################################
#### 1.6 - Gradient and Laplacian Evaluation at beta.star
########################################################################

alpha.cent <<- matrix(c(-0.5848635,-0.3187415,-0.2912184,0.606722,0.4416645),1,5)
alpha.cent.sq <<- 1.091649
alpha.p.cent <<- -0.2739236
phi.cent <<- 0.4088625

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
