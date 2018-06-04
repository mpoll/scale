############################################
############################################
######       Scale Algorithm         #######
######    Logistic Specification     #######
######         User Example          #######
######  Last Updated: 04/06/2018 MP  #######
############################################
############################################

airline.logistic.example.2   <- function(){
    
    ########################################################################
    #### 0.1 - Save current seed and set data seed
    ########################################################################
    
    curr.seed <- .Random.seed
    set.seed(1)
    load("airline_raw.RData")
    
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
    
    airline <- airline[sample.int(dsz, dsz, replace = FALSE),]
    examp.design <- matrix(0,dsz,dimen); examp.design[,1] <- 1; for(j in 2:dimen){examp.design[,j] <- as.numeric(airline[,j])}; examp.design <<- examp.design
    examp.data <<- airline[,1]
    
    ########################################################################
    #### 1.2 - Centering Computation
    ########################################################################
    
    subsets <- 13 # Total Number of subsets
    subset.size <- dsz/subsets # Subset size
    data.partitions <- 1 + c(0:(subsets-1))*subset.size
    
    ### Execute seperate GLM logisitc fit for specified subsets, saving output
    
    subsets.fit <- list()
    for(i in 1:subsets){
        subset.centering <- logistic.centering(data.partitions[i],subset.size)
        subsets.fit[[i]] <- list(subset=i,idx.start=data.partitions[i],idx.end=data.partitions[i]+subset.size-1,center=subset.centering$center,precon=subset.centering$precon,precon.inv=subset.centering$precon.inv)
        print(i)}
    
    #### 1.2.1 - Centering Value
    
    subset.centers <- matrix(0,subsets,dimen); for(i in 1:subsets){subset.centers[i,] <- subsets.fit[[i]]$center}
    
    beta.star       <<- matrix(as.vector(colMeans(subset.centers)),1,dimen) # What is the choosen centering value in the original
    
    #### 1.2.2 - Preconditioning Matrix
    
    precon.inv <- matrix(0,dimen,dimen); for(i in 1:subsets){precon.inv <- precon.inv + subsets.fit[[i]]$precon.inv}; precon <- solve(precon.inv)
    
    n.sigma <<- diag(sqrt(as.vector(diag(precon)))); n.sigma.inv <<- solve(n.sigma) # Parameterise preconditioning matrix (inverse Fischers information) and specification of inverse preconditioning matrix
    
    #############################################
    #### 1.3 - Data Recaller
    #############################################
    
    datum     <<- function(datum.idx){
        datum.x     <- as.numeric(examp.design[datum.idx,])
        datum.y     <- as.numeric(examp.data[datum.idx])
        list(idx=datum.idx,x=datum.x,y=datum.y)}
    
    ########################################################################
    #### 1.4 - Gradient and Laplacian Evaluation at beta.star
    ########################################################################
    
    # Partition data set into executable pieces
    
    subsets <- 13 # Total Number of subsets
    subset.size <- dsz/subsets # Subset size
    data.partitions <- 1 + c(0:(subsets-1))*subset.size
    
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
    alpha.cent.sq <<- (alpha.cent)%*%t(alpha.cent)
    alpha.p.cent <<- sum(subset.alpha.p.cent)
    phi.cent <<- (alpha.cent.sq+alpha.p.cent)/2
    
    ########################################################################
    ########################################################################
    ### -2- Computation of Subsampling bounding functionals
    ########################################################################
    ########################################################################
    
    # Partition data set into executable pieces
    
    subsets <- 13 # Total Number of subsets
    subset.size <- dsz/subsets # Subset size
    data.partitions <- 1 + c(0:(subsets-1))*subset.size
    
    # Compute extremal points for each dimension for each partition
    
    subsets.extreme <- list()
    for(i in 1:subsets){
        subset.ext <- data.extrema(data.partitions[i],subset.size)
        subsets.extreme[[i]] <- list(subset=i,idx.start=data.partitions[i],idx.end=data.partitions[i]+subset.size-1,design.max=subset.ext$design.max,design.min=subset.ext$design.min)
        print(i)}
    
    # Extract from list extremal points for each dimension
    
    subset.design.min <- subset.design.max <- matrix(0,subsets,dimen); for(i in 1:subsets){subset.design.min[i,] <- subsets.extreme[[i]]$design.min; subset.design.max[i,]  <- subsets.extreme[[i]]$design.max}
    
    # Compute global extremal points for each dimension
    
    design.min <- as.numeric(apply(subset.design.min,2,min))
    design.max <- as.numeric(apply(subset.design.max,2,max))
    
    # Compute worst possible datum by considering combinations within extremal boundaries, and then compute maximal log.like and lap.like for these points (i.e. data.extrema() part II)
    
    data.extrema.2(design.min,design.max)
    
    ############################################
    ############################################
    ### -3- End Function Specification
    ############################################
    ############################################
    
    .Random.seed <<- curr.seed
    
}
