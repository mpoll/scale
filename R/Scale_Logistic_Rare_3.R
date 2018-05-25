########################################################################
########################################################################
######                      Scale Algorithm                      #######
######                  Logistic Specification                   #######
######                      Rare Data Example                     #######
######               Last Updated: 30/05/2018 MP                 #######
########################################################################
########################################################################

rare.logistic.example.3    <- function(){

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
#### 1.1 - Size, dimension, and data set
########################################################################

dsz             <<- 10^7                    # Size of data set
dimen           <<- 4                        # Dimensionality (>1)
beta.truth      <<- c(0,2,-2,2)            # True parameters

########################################################################
#### 1.2 - Generate pseudo-data
########################################################################

library('msm')                            # Package to generate truncated Normals
design.thresh   <<- 1                     # Thresholding of truncated Normals

examp.design <<- cbind(rep(1,dsz),rbinom(dsz,1,10/dsz),matrix(rtnorm((dimen-2)*dsz,mean=0,sd=1,lower=-design.thresh,upper=design.thresh),dsz,dimen-2))
generate.datum <- function(x){datum.z <- exp(sum(beta.truth*x)); rbinom(1,1,datum.z/(1+datum.z))}
examp.data <<- apply(examp.design,1,generate.datum)
sum(examp.data)

########################################################################
#### 1.2 - Centering Computation
########################################################################

initialisation <- logistic.centering()

#### 1.2.1 - Centering Value

beta.star       <<- matrix(as.vector(initialisation$center),1,dimen) # What is the choosen centering value in the original

#### 1.2.2 - Preconditioning Matrix

n.sigma         <<- diag(sqrt(as.vector(diag(initialisation$precon)))); n.sigma.inv <<- solve(n.sigma) # Parameterise preconditioning matrix (inverse Fischers information) and specification of inverse preconditioning matrix

########################################################################
#### 1.3 - Data Recaller
########################################################################

datum     <<- function(datum.idx){
    datum.x     <- examp.design[datum.idx,]
    datum.y     <- examp.data[datum.idx]
    list(idx=datum.idx,x=datum.x,y=datum.y)}

########################################################################
#### 1.4 - Gradient and Laplacian Evaluation at beta.star
########################################################################

data.init()

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
