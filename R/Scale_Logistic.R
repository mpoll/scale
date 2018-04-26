########################################################################
########################################################################
######                      Scale Algorithm                      #######
######        Algorithm Specification for Logistic Data          #######
######              Last Updated: 22/04/2018 MP                  #######
########################################################################
########################################################################

########################################################################
########################################################################
### 1 - Functions transforming from Beta to Eta Space to Beta Space
########################################################################
########################################################################

scale.transform       <- function(beta.pos){matrix(c(diag(n.sigma.inv)*(beta.pos-beta.star)),nrow=1)} # Transformation from beta to eta space
un.scale.transform    <- function(eta.pos){matrix(c(diag(n.sigma)*eta.pos+beta.star),nrow=1)} # Transformation from eta to beta space

########################################################################
########################################################################
### 2 - Functions to recall and evaluate multiple data points
########################################################################
########################################################################

########################################################################
#### 2.1 - Recall and evaluate a single datum
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
#### 2.1 - Recall and evaluate multiple data points
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
#### 3.1 - Function for evaluating data at centering point for values
########################################################################

data.init       <- function(beta.eval=beta.star,datum.idx.start=1,length.idxs=dsz){
    evaluate <- data.eval(beta.eval,c(datum.idx.start:(datum.idx.start+length.idxs-1)))
    alpha.cent <<- matrix(evaluate$grad.log.pi,nrow=1)
    alpha.cent.sq <<- (alpha.cent)%*%t(alpha.cent)
    alpha.p.cent <<- sum(evaluate$lap.log.pi)
    phi.cent <<- (alpha.cent.sq+alpha.p.cent)/2
    list(alpha.cent=alpha.cent,alpha.cent.sq=alpha.cent.sq,alpha.p.cent=alpha.p.cent,phi.cent=phi.cent)}

########################################################################
#### 3.2 - Computation of Subsampling bounding functionals
########################################################################

data.extrema    <- function(){
    ##### 3.2.1 - Precompute the most extreme data values
    if(exists("design.thresh")==TRUE){if(length(design.thresh)==1){design.min <- c(1,rep(-design.thresh,dimen-1)); design.max <- c(1,rep(design.thresh,dimen-1))}}
    if(exists("design.thresh")==FALSE){design.min <- apply(examp.design,2,min); design.max <- apply(examp.design,2,max)}
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
    list(dim.grad.max=dim.grad.max)}

########################################################################
#### 3.3 - Initialise ScaLE algorithm
########################################################################

scale_logistic_init <- function(fnm="logistic_default.RData"){data.init(); data.extrema(); save(file=fnm,dim.grad.max=dim.grad.max,alpha.cent=alpha.cent,alpha.cent.sq=alpha.cent.sq,alpha.p.cent=alpha.p.cent,phi.cent=phi.cent)}

########################################################################
########################################################################
### 4 - Functional to execute exact version of ScaLE for logistic data
########################################################################
########################################################################

scale_logistic <- function(fnm="logistic_default.RData",p.num=2^10,t.inc=0.1,T.extend=0.1,run.length=100000,ss.size=10,ss.on=TRUE,seed.default=1,default.data=TRUE){
    set.seed(seed.default)
    if(default.data==TRUE){big.logistic.example()}else{user.logistic.example()}
    save(file=fnm,dim.grad.max=dim.grad.max,alpha.cent=alpha.cent,alpha.cent.sq=alpha.cent.sq,alpha.p.cent=alpha.p.cent,phi.cent=phi.cent,curr.seed=.Random.seed)
    time.elapse <- system.time(simn <<- scale_exact(p.num=p.num,t.inc=t.inc,T.fin=T.extend,ss.phi=ss.phi,ss.phiC=ss.phiC,dimen,scale.transform,un.scale.transform,T.start=0,x.init=NULL,ss.size=ss.size,ess.thresh=0.5,resamp.method=resid.resamp,neg.wei.mech=scale_zero.wei,prev.simn=NULL,progress.check=FALSE,phi.record=FALSE,resamp.freq=p.num-1,theta=NULL,p.path.renew=p.path.renew))[3]
    save(file=fnm,simn=simn,dim.grad.max=dim.grad.max,alpha.cent=alpha.cent,alpha.cent.sq=alpha.cent.sq,alpha.p.cent=alpha.p.cent,phi.cent=phi.cent,curr.seed=.Random.seed,time.elapse=time.elapse); print(time.elapse)
    while(simn$T.fin < run.length){time.inc <- system.time(simn <<- scale_extend(simn,t.inc,T.extend))[3]
    time.elapse <- time.elapse + time.inc
        save(file=fnm,simn,time.elapse); print(time.elapse)}}

