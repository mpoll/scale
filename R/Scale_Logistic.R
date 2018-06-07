########################################################################
########################################################################
######                      Scale Algorithm                      #######
######        Algorithm Specification for Logistic Data          #######
######              Last Updated: 04/06/2018 MP                  #######
########################################################################
########################################################################

########################################################################
########################################################################
### 1 - Functions transforming from Beta to Eta Space to Beta Space
########################################################################
########################################################################

scale_transform       <- function(beta.pos){matrix(c(diag(n.sigma.inv)*(beta.pos-beta.star)),nrow=1)} # Transformation from beta to eta space
un_scale_transform    <- function(eta.pos){matrix(c(diag(n.sigma)*eta.pos+beta.star),nrow=1)} # Transformation from eta to beta space

########################################################################
########################################################################
### 2 - Functions to recall and evaluate multiple data points
########################################################################
########################################################################

########################################################################
#### 2.1 - Recall and evaluate multiple data points
########################################################################

datum.eval   <- function(beta.pos,data.idx){
    #### CALL DATA
    datum.call      <- datum(data.idx)
    #### EVALUATE
    datum.z         <- exp(sum(datum.call$x*beta.pos))
    log.pi          <- datum.call$x*beta.pos*datum.call$y-log(1+datum.z)
    grad.log.pi     <- diag(n.sigma)*datum.call$x*(datum.call$y-datum.z/(1+datum.z))
    lap.log.pi      <- -(diag(n.sigma)*datum.call$x)^2*datum.z/(1+datum.z)^2
    list(log.pi=log.pi,grad.log.pi=grad.log.pi,lap.log.pi=lap.log.pi)}

########################################################################
#### 2.2 - Recall multiple data points to output specified subset of the design matrix and data
########################################################################

data.subset <- function(datum.idx.start=1,length.idxs=dsz){
    if(exists("examp.design")==FALSE){
        call.data <- datum(datum.idx.start)
        subset.design <- t(as.matrix(call.data$x)); subset.data <- call.data$y
        if(length.idxs>1){for(i in (datum.idx.start+1):(datum.idx.start+length.idxs-1)){
            call.data <- datum(i)
            subset.design <- rbind(subset.design,call.data$x); subset.data <- c(subset.data,call.data$y)    }}}
    if(exists("examp.design")==TRUE){
        subset.design  <- examp.design[datum.idx.start:(datum.idx.start+length.idxs-1),]
        subset.data <- examp.data[datum.idx.start:(datum.idx.start+length.idxs-1)]}
    list(subset.design=subset.design,subset.data=subset.data)}

########################################################################
#### 2.3 - Recall and evaluate multiple data points
########################################################################

data.eval       <- function(beta.pos,data.idxs,factor=1){ # Factor is a multiple of output grad log and lap log
    grad.log.pi <- lap.log.pi <- numeric(dimen)
    for(i in 1:length(data.idxs)){
        datum.call          <- datum.eval(beta.pos,data.idxs[i])
        grad.log.pi         <- grad.log.pi + datum.call$grad.log.pi
        lap.log.pi          <- lap.log.pi + datum.call$lap.log.pi}
    list(grad.log.pi=factor*grad.log.pi,lap.log.pi=factor*lap.log.pi)}

########################################################################
########################################################################
### 3 - Functions to initialise algorithm (centering values + phi bds)
########################################################################
########################################################################

########################################################################
#### 3.0 - Function for using GLM to compute the centering value and covariance
########################################################################

logistic.centering      <- function(datum.idx.start=1,length.idxs=dsz){
        call.data <- data.subset(datum.idx.start,length.idxs)
        subset.design.mat <- as.matrix(call.data$subset.design[,-1,drop=FALSE])
        subset.data <- call.data$subset.data
        subset.glm <- glm(formula = subset.data ~ subset.design.mat, family = binomial(link='logit'))
        center <- subset.glm$coef
        precon <- vcov(subset.glm)
        list(center=center,precon=precon,precon.inv=solve(precon),logistic.glm=subset.glm,det.metric=det(precon))
}

########################################################################
#### 3.1 - Function for evaluating data at centering point for values
########################################################################

data.init       <- function(beta.eval=beta.star,datum.idx.start=1,length.idxs=dsz){
    evaluate <- data.eval(beta.eval,c(datum.idx.start:(datum.idx.start+length.idxs-1)))
    alpha.cent <<- matrix(as.numeric(evaluate$grad.log.pi),nrow=1)
    alpha.cent.sq <<- (alpha.cent)%*%t(alpha.cent)
    alpha.p.cent <<- sum(as.numeric(evaluate$lap.log.pi))
    phi.cent <<- (alpha.cent.sq+alpha.p.cent)/2
    list(alpha.cent=alpha.cent,alpha.cent.sq=alpha.cent.sq,alpha.p.cent=alpha.p.cent,phi.cent=phi.cent)}

# More efficient initialisation in case data is fully available in form of examp.design and examp.data and is partitioned due to size

data.init.2       <- function(beta.eval=beta.star,datum.idx.start=1,length.idxs=dsz){
    call.data <- data.subset(datum.idx.start,length.idxs)
    subset.design <- call.data$subset.design
    subset.data <- call.data$subset.data
    subset.zs <- as.numeric(apply(subset.design,1,function(x) exp(sum(x*beta.eval))))
    subset.grad.log.pi <- sapply(seq_len(length.idxs),function(i) diag(n.sigma)*subset.design[i,]*(subset.data[i]-subset.zs[i]/(1+subset.zs[i])))
    subset.lap.log.pi <- sapply(seq_len(length.idxs),function(i) -(diag(n.sigma)*subset.design[i,])^2*subset.zs[i]/(1+subset.zs[i])^2)
    grad.log.pi <- as.numeric(apply(subset.grad.log.pi,1,sum))
    lap.log.pi <- as.numeric(apply(subset.lap.log.pi,1,sum))
    alpha.cent <<- matrix(grad.log.pi,nrow=1)
    alpha.cent.sq <<- (alpha.cent)%*%t(alpha.cent)
    alpha.p.cent <<- sum(lap.log.pi)
    phi.cent <<- (alpha.cent.sq+alpha.p.cent)/2
    list(alpha.cent=alpha.cent,alpha.cent.sq=alpha.cent.sq,alpha.p.cent=alpha.p.cent,phi.cent=phi.cent)}

########################################################################
#### 3.2 - Computation of Subsampling bounding functionals
########################################################################

data.extrema    <- function(datum.idx.start=1,length.idxs=dsz){
    ##### 3.2.1 - Precompute the most extreme data values
    if(exists("design.thresh")==TRUE){if(length(design.thresh)==1){design.min <- c(1,rep(-design.thresh,dimen-1)); design.max <- c(1,rep(design.thresh,dimen-1))}else{design.min <- c(1,-design.thresh[2:dimen]); design.max <- design.thresh}}
    if(exists("design.thresh")==FALSE){subset.design <- data.subset(datum.idx.start,length.idxs)$subset.design; design.min <- apply(subset.design,2,min); design.max <- apply(subset.design,2,max)}
    extrema.design <- permutations(n=2,r=dimen,v=c(0,1),repeats.allowed=TRUE); for(i in 1:dimen){extrema.design[,i] <- (1-extrema.design[,i])*design.min[i] + extrema.design[,i]*design.max[i]}; extrema.design <- rbind(extrema.design,extrema.design)
    data.extrema <- matrix(c(rep(0,dim(extrema.design)[1]/2),rep(1,dim(extrema.design)[1]/2)),dim(extrema.design)[1],1)
    
    #### 3.2.2 - Precompute the maxima Gradient Datum Log likelihood and Gradient Laplacian Datum Log likelihood for each datum
    data.precompute <- array(0,c(dim(extrema.design)[1],dimen,dimen,2)); for(i in 1:dim(extrema.design)[1]){for(dim.j in 1:dimen){
            X.row <- extrema.design[i,] # Determine row of design matrix
            data.precompute[i,1:dimen,dim.j,1] <- abs(-n.sigma[dim.j,dim.j]*diag(n.sigma)*X.row[dim.j]*X.row)/4 # maxima Gradient Datum Log likelihood
            data.precompute[i,1:dimen,dim.j,2] <- abs(-(n.sigma[dim.j,dim.j])^2*diag(n.sigma)*(X.row[dim.j])^2*X.row)*1/(6*sqrt(3))}} # Maxima Gradient Laplacian Datum Log likelihood

    #### 3.2.3 - Precompute the maxima Gradient Datum Log likelihood and Laplacian Datum Log likelihood over all data
    dim.grad.max  <- array(0,c(dimen,dimen,2)); for(dim.j in 1:dimen){for(dim.k in 1:dimen){
        dim.grad.max[dim.j,dim.k,1] <- max(data.precompute[,dim.j,dim.k,1])     # Maximal Gradient Log Likelihood
        dim.grad.max[dim.j,dim.k,2] <- max(data.precompute[,dim.j,dim.k,2])}}   # Maximal Laplacian Log Likelhood

    #### 3.2.4 - Write to global scope dim.grad.max and output
    dim.grad.max <<- dim.grad.max
    list(design.min=design.min,design.max=design.max,dim.grad.max=dim.grad.max)}

# Functional to compute dim.grad.max when design.min and design.max specfied

data.extrema.2     <- function(design.min,design.max){
    ##### 3.2.b.1 - Precompute the most extreme data values
    extrema.design <- permutations(n=2,r=dimen,v=c(0,1),repeats.allowed=TRUE); for(i in 1:dimen){extrema.design[,i] <- (1-extrema.design[,i])*design.min[i] + extrema.design[,i]*design.max[i]}; extrema.design <- rbind(extrema.design,extrema.design)
    data.extrema <- matrix(c(rep(0,dim(extrema.design)[1]/2),rep(1,dim(extrema.design)[1]/2)),dim(extrema.design)[1],1)
   
   #### 3.2.b.2 - Precompute the maxima Gradient Datum Log likelihood and Gradient Laplacian Datum Log likelihood for each datum
   data.precompute <- array(0,c(dim(extrema.design)[1],dimen,dimen,2)); for(i in 1:dim(extrema.design)[1]){for(dim.j in 1:dimen){
       X.row <- extrema.design[i,] # Determine row of design matrix
       data.precompute[i,1:dimen,dim.j,1] <- abs(-n.sigma[dim.j,dim.j]*diag(n.sigma)*X.row[dim.j]*X.row)/4 # maxima Gradient Datum Log likelihood
       data.precompute[i,1:dimen,dim.j,2] <- abs(-(n.sigma[dim.j,dim.j])^2*diag(n.sigma)*(X.row[dim.j])^2*X.row)*1/(6*sqrt(3))}} # Maxima Gradient Laplacian Datum Log likelihood
   
   #### 3.2.b.3 - Precompute the maxima Gradient Datum Log likelihood and Laplacian Datum Log likelihood over all data
   dim.grad.max  <- array(0,c(dimen,dimen,2)); for(dim.j in 1:dimen){for(dim.k in 1:dimen){
       dim.grad.max[dim.j,dim.k,1] <- max(data.precompute[,dim.j,dim.k,1])     # Maximal Gradient Log Likelihood
       dim.grad.max[dim.j,dim.k,2] <- max(data.precompute[,dim.j,dim.k,2])}}   # Maximal Laplacian Log Likelhood
   
   #### 3.2.4 - Write to global scope dim.grad.max and output
   dim.grad.max <<- dim.grad.max
    list(design.min=design.min,design.max=design.max,dim.grad.max=dim.grad.max)}

########################################################################
########################################################################
### 4 - Functional to execute exact version of ScaLE for logistic data
########################################################################
########################################################################

scale_logistic <- function(fnm="logistic_default.RData",p.num=2^10,t.inc=0.01,T.extend=0.01,run.length=25,ss.size=10,ss.on=TRUE,seed.default=1,data.precompute=TRUE){
    # Set seed to default seed
    set.seed(seed.default)
    # Data precomputed by default. If not, compute using the following files
    if(data.precompute=="big.logistic.example"){big.logistic.example()}
    if(data.precompute=="rare.logistic.example"){rare.logistic.example()}
    if(data.precompute=="airline.logistic.example"){airline.logistic.example()}
    # Save current seed and current data
    curr.seed <- .Random.seed
    save(file=fnm,curr.seed=curr.seed)
    # Run initial Scale run
    time.elapse <- system.time(simn <<- scale_exact(p.num=p.num,t.inc=t.inc,T.fin=T.extend,ss.phi=ss.phi,ss.phiC=ss.phiC,dimen,scale_transform,un_scale_transform,T.start=0,x.init=NULL,ss.size=ss.size,ess.thresh=0.5,resamp.method=resid.resamp,neg.wei.mech=scale_zero.wei,prev.simn=NULL,progress.check=FALSE,phi.record=FALSE,resamp.freq=p.num-1,theta=NULL,p.path.renew=p.path.renew))[3]
    # Save current seed and current data
    curr.seed <- .Random.seed
    save(file=fnm,simn=simn,curr.seed=curr.seed,time.elapse=time.elapse); print(time.elapse)
    # Run Scale outputting current seed and current simulation
    while(simn$T.fin < run.length){time.inc <- system.time(simn <<- scale_extend(simn,t.inc,T.extend))[3]
        time.elapse <- time.elapse + time.inc
        curr.seed <- .Random.seed
        save(file=fnm,simn=simn,curr.seed=curr.seed,time.elapse=time.elapse); print(time.elapse)}}

#############################################
#############################################
#### 5 - Functions to combine multiply split data for inititialisation
#############################################
#############################################

#############################################
#### 5.1 - Pass 1: Computing centering point and preconditioning matrix combination
#############################################

scale_fileapply_pass_1 <- function(fnm="combined_data_1.RData",pathformstart="scale_precomp_1_",directory=getwd()){
    file.names <- dir(path = directory, pattern = pathformstart, full.names=TRUE); combined.data_1 <- list() # Initialise combiner
    for(i in 1:length(file.names)){ # Loop over all file names
        flist <- NULL; flist <- ls() # Surmise current data
        load(file.names[i], envir = .GlobalEnv) # Load data
        run <- as.integer(substring(file.names[i],nchar(directory)+nchar(pathformstart)+2,nchar(file.names[i])-nchar(".RData"))) # Find run number
        combined.data_1[[i]]    <- list(run=run,subset=subset,idx.start=idx.start,idx.end=idx.end,center=center,precon=precon,precon.inv=precon.inv) # Apply data function
        rm(list = ls()[ls() %in% flist == FALSE])} # Remove any loaded/generated data
    save(file=fnm,combined.data_1=combined.data_1)}

#############################################
#### 5.2 - Pass 2: Computing centering values and data extrema at centering point
#############################################

scale_fileapply_pass_2  <- function(fnm="combined_data_2.RData",pathformstart="scale_precomp_2_",directory=getwd()){
    file.names <- dir(path = directory, pattern = pathformstart, full.names=TRUE); combined.data_2 <- list() # Initialise combiner
    for(i in 1:length(file.names)){ # Loop over all file names
        flist <- NULL; flist <- ls() # Surmise current data
        load(file.names[i], envir = .GlobalEnv) # Load data
        run <- as.integer(substring(file.names[i],nchar(directory)+nchar(pathformstart)+2,nchar(file.names[i])-nchar(".RData"))) # Find run number
        combined.data_2[[i]]    <- list(run=run,subset=i,idx.start=idx.start,idx.end=idx.end,alpha.cent=alpha.cent,alpha.p.cent=alpha.p.cent,design.max=design.max,design.min=design.min) # Apply data function
        rm(list = ls()[ls() %in% flist == FALSE])} # Remove any loaded/generated data
    save(file=fnm,combined.data_2=combined.data_2)}
