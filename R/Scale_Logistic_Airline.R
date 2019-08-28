########################################################################
########################################################################
######                      Scale Algorithm                      #######
######                  Logistic Specification                   #######
######                      Airline Data Example                 #######
######               Last Updated: 25/08/2019 MP                 #######
########################################################################
########################################################################

airline.logistic.example   <- function(){
    
    ########################################################################
    #### 0.1 - Save current seed and set data seed
    ########################################################################

    curr.seed <- .Random.seed
    set.seed(1)
    load("data_airline_raw.RData") # Raw data set
    #load("data_airline.RData") # Data set pre-computed using the below functionals. Note computation of the GLM fit and control variates is broken into parts subject to the constraints of the developers computer.
    
    ########################################################################
    ########################################################################
    ### -1- Specification of Data Set
    ########################################################################
    ########################################################################
    
    ########################################################################
    #### 1.1 - Size, dimension, and data set
    ########################################################################
    
    dsz             <<- 120748239                    # Size of data set
    dimen           <<- 4                            # Dimensionality (>1)

    ########################################################################
    #### 1.2 - Data
    ########################################################################
    
    ### If airline.RData is not loaded then the following steps are required
    
    ### Note the "Airline" data set ("airline_raw.RData has to be pre-processed so that examp.design (the design matrix for the example) and examp.data (the data for the example) are in the required form for execution.
    
    ### The airline data set is highly structured, and so it is advisable to permute the data before compute GLM fit, and associated control variates
    
    ### Column 1 of the airline data set corresponds to the data, the remainder the covariates
    
    airline <- airline[sample.int(dsz, dsz, replace = FALSE),]
    examp.design <- matrix(0,dsz,dimen); examp.design[,1] <- 1; for(j in 2:dimen){examp.design[,j] <- as.numeric(airline[,j])}; examp.design <<- examp.design
    examp.data <<- airline[,1]


    ########################################################################
    #### 1.2 - Centering Computation
    ########################################################################
    
    subsets <- 13 # Total Number of subsets
    subset.size <- dsz/subsets # Subset size
    data.partitions <- 1 + c(0:(subsets-1))*subset.size
    
    ### Execute seperate GLM logistic fit for specified subsets, saving output
    
    subsets.fit <- list()
    for(i in 1:subsets){
        subset.centering <- logistic.centering(data.partitions[i],subset.size)
        subsets.fit[[i]] <- list(subset=i,idx.start=data.partitions[i],idx.end=data.partitions[i]+subset.size-1,center=subset.centering$center,precon=subset.centering$precon,precon.inv=subset.centering$precon.inv)
        print(i)}
    
    #### 1.2.1 - Centering Value
    
    subset.centers <- matrix(0,subsets,dimen); for(i in 1:subsets){subset.centers[i,] <- subsets.fit[[i]]$center}
    
    beta.star       <<- matrix(as.vector(colMeans(subset.centers)),1,dimen) # What is the choosen centering value in the original
    
    #### 1.2.2 - Preconditioning Matrix
    
    precon.inv <<- matrix(0,dimen,dimen); for(i in 1:subsets){precon.inv <- precon.inv + subsets.fit[[i]]$precon.inv};
    precon <<- solve(precon.inv)
    
    n.sigma <<- solve(diag(sqrt(diag(precon.inv)))); n.sigma.inv <<- solve(n.sigma) # Parameterise preconditioning matrix (inverse Fischers information) and specification of inverse preconditioning matrix
    
    ########################################################################
    #### 1.3 - Data Recaller
    ########################################################################

    datum     <<- function(datum.idx){
        datum.x         <- as.numeric(examp.design[datum.idx,])
        datum.y         <- as.numeric(examp.data[datum.idx])
        list(idx=datum.idx,x=datum.x,y=datum.y)}
    
    ########################################################################
    #### 1.4 - Gradient and Laplacian Evaluation at beta.star
    ########################################################################
    
    # Compute alpha + alpha.p for centering point for each partition
    
    subsets.control <- list()
    for(i in 1:subsets){
        subset.init <- data.init.2(beta.star,data.partitions[i],subset.size)
        subsets.control[[i]] <- list(subset=i,idx.start=data.partitions[i],idx.end=data.partitions[i]+subset.size-1,alpha.cent=subset.init$alpha.cent,alpha.p.cent=subset.init$alpha.p.cent)
        print(i)}
    
    # Extract from list alpha + alpha.p
    
    subset.alpha.cent <- matrix(0,subsets,dimen); subset.alpha.p.cent <- numeric(subsets); for(i in 1:subsets){subset.alpha.cent[i,] <- subsets.control[[i]]$alpha.cent; subset.alpha.p.cent[i]  <- subsets.control[[i]]$alpha.p.cent}
    
    # Compute total alpha + alpha.p (summate)
    
    alpha.cent <<- matrix(as.numeric(apply(subset.alpha.cent,2,sum)),nrow=1)
    alpha.cent.bds <<- sum((2*alpha.cent)^2)^(1/2)
    alpha.cent.sq <<- (alpha.cent)%*%t(alpha.cent)
    alpha.p.cent <<- sum(subset.alpha.p.cent)
    phi.cent <<- (alpha.cent.sq+alpha.p.cent)/2
    Hessian.bound <<- sum((diag(n.sigma)^2)/4)
    
    ########################################################################
    ########################################################################
    ### -2- End Function Specification
    ########################################################################
    ########################################################################
    
    .Random.seed <<- curr.seed
    
}
