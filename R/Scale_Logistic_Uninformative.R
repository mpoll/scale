########################################################################
########################################################################
######                      Scale Algorithm                      #######
######                  Logistic Specification                   #######
######                Uninformative Data Example                 #######
######               Last Updated: 25/08/2019 MP                 #######
########################################################################
########################################################################

uninformative.logistic.example    <- function(){
    
    ########################################################################
    #### 0.1 - Save current seed and set data seed
    ########################################################################
    
    curr.seed <- .Random.seed
    set.seed(1)
    #load("data_uninformative.RData") # Data set pre-computed using the below functionals
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
    
    design.threshold.1   <<- 0.001             # Thresholding of truncated Normals column 2
    design.threshold.2   <<- 1                # Thresholding of truncated Normals remaining columns
    design.thresh        <<- c(1,design.threshold.1,rep(design.threshold.2,dimen-2))
    
    examp.design <<- matrix(1,nrow=dsz,ncol=1); for(j in 2:dimen){examp.design <<- cbind(examp.design,msm::rtnorm(dsz,mean=0,sd=1,lower=-design.thresh[j],upper=design.thresh[j]))}
    generate.datum <- function(x){datum.z <- exp(sum(beta.truth*x)); rbinom(1,1,datum.z/(1+datum.z))}
    examp.data <<- apply(examp.design,1,generate.datum)
    
    ########################################################################
    #### 1.2 - Centering Computation
    ########################################################################
    
    initialisation <- logistic.centering()
    
    #### 1.2.1 - Centering Value
    
    beta.star       <<- matrix(as.vector(initialisation$center),1,dimen) # What is the choosen centering value in the original
    
    #### 1.2.2 - Preconditioning Matrix
    
    precon <<- initialisation$precon
    precon.inv <<- initialisation$precon.inv
    
    n.sigma         <<- solve(diag(sqrt(diag(precon.inv)))); n.sigma.inv <<- solve(n.sigma) # Parameterise preconditioning matrix (inverse Fischers information) and specification of inverse preconditioning matrix
    
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
    
    intensity.initialisation <- data.init()
    
    alpha.cent <<- intensity.initialisation$alpha.cent
    alpha.cent.bds <<- intensity.initialisation$alpha.cent.bds
    alpha.cent.sq <<- intensity.initialisation$alpha.cent.sq
    alpha.p.cent <<- intensity.initialisation$alpha.p.cent
    phi.cent <<- intensity.initialisation$phi.cent
    
    H.mat <- matrix(0,dimen,dimen)
    for(i in 1:dimen){for(j in 1:dimen){H.mat[i,j] <- n.sigma[i,i]*n.sigma[j,j]*design.thresh[i]*design.thresh[j]/4}}

    Hessian.bound <<- max(eigen(H.mat)$values)
    
    ########################################################################
    ########################################################################
    ### -3- End Function Specification
    ########################################################################
    ########################################################################
    
    .Random.seed <<- curr.seed
    
}
