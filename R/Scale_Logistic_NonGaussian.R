########################################################################
########################################################################
######                      Scale Algorithm                      #######
######                  Logistic Specification                   #######
######         Non Gaussian (Extremal Small N) Example           #######
######               Last Updated: 30/08/2019 MP                 #######
########################################################################
########################################################################

small.logistic.example    <- function(){
    
    ########################################################################
    #### 0.1 - Save current seed and set data seed
    ########################################################################
    
    curr.seed <- .Random.seed
    set.seed(1)
    #load("data_small.RData") # Data set pre-computed using the below functionals
    ########################################################################
    ########################################################################
    ### -1- Specification of Data Set
    ########################################################################
    ########################################################################
    
    ########################################################################
    #### 1.1 - Size, dimension, and data set
    ########################################################################
    
    dsz             <<- 10                       # Size of data set
    dimen           <<- 2                       # Dimensionality (>1)
    
    ########################################################################
    #### 1.2 - Generate pseudo-data
    ########################################################################
    
    examp.design <<- matrix(c(rep(1,dsz),(-1)^odd(1:(dsz*(dimen-1)))/c(1:dsz)),nrow=dsz,ncol=2,byrow=FALSE)
    examp.data <<- c(rep(1,2),rep(0,8))
    
    design.thresh <<- apply(abs(examp.design),2,max)
    
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
