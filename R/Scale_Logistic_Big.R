########################################################################
########################################################################
######                      Scale Algorithm                      #######
######                  Logistic Specification                   #######
######                      Big Data Example                     #######
######               Last Updated: 25/05/2018 MP                 #######
########################################################################
########################################################################

big.logistic.example    <- function(){

########################################################################
########################################################################
### -1- Specification of Data Set
########################################################################
########################################################################

########################################################################
#### 1.1 - Size, dimension, and parameters
########################################################################

dsz             <<- 10^10                    # Size of data set
dimen           <<- 5                        # Dimensionality (>1)
beta.truth      <<- c(1,1,-1,2,-2)           # True parameters
n.sigma         <<- diag(rep(1,dimen))*dsz^(-1/2); n.sigma[1,1] <<- n.sigma[1,1]/2; n.sigma.inv <<- solve(n.sigma)  # Parameterise preconditioning matrix (inverse Fischers information) and specification of inverse preconditioning matrix

########################################################################
#### 1.2 - Centering
########################################################################

beta.star       <<- matrix(beta.truth,1,dimen) # What is the choosen centering value in the original space

########################################################################
#### 1.3 - Data Recaller
########################################################################

library('msm')                            # Package to generate truncated Normals
design.thresh   <<- 1                    # Thresholding of truncated Normals

datum     <<- function(datum.idx){ # Function to re-generate data (without storage) - won't work with seeds larger than 2^53-1
    old.seed <- .Random.seed; on.exit({ .Random.seed <<- old.seed })
    seed.indexing <- datum.idx%/%2147483647
    datum.seeding <- datum.idx%%2147483647
    set.seed(datum.seeding); runif((2*dimen-1)*seed.indexing) # junk indexing by the number of random operations used to simulate data.gen
    datum.x     <- c(1,rtnorm(dimen-1,mean=0,sd=1,lower=-design.thresh,upper=design.thresh))
    datum.z     <- exp(sum(beta.truth*datum.x))
    datum.y     <- rbinom(1,1,datum.z/(1+datum.z))
    list(idx=datum.idx,x=datum.x,y=datum.y,seed=datum.seeding,seed.idx=seed.indexing)}

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

}
