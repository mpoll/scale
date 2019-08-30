########################################################################
########################################################################
######                      Scale Algorithm                      #######
######              Last Updated: 25/08/2019 MP                  #######
########################################################################
########################################################################

########################################################################
########################################################################
### 0 - Libraries, Global Functions & Euler
########################################################################
########################################################################

mpfrthr <- 9; mpfrpbn <- 10 # mpfr - high  precision kick in levels for alternating series
mat.sort.c <- function(mat,mat.sz,n) {mat[rank(mat[,n],ties.method="first"),] <- mat[1:mat.sz,]; return(mat)} # Define matrix sort (ranking along column)
mat.sort.r <- function(mat,mat.sz,n) {mat[,rank(mat[n,],ties.method="first")] <- mat[,1:mat.sz]; return(mat)} # Define matrix sort (ranking along a row)
ins.sort.c <- function(vec,sor.mat,sor.mat.sz,n){ # Insert a vector into a sorted matrix (ranking along a column)
    r.min <- 1; r.max <- sor.mat.sz; tau <- vec[n]
    if(tau > sor.mat[r.max,n]){r.min <- r.max <- sor.mat.sz+1}else{
        while(r.max > r.min){
            r.mid <- floor(0.5*(r.max+r.min))
            if(tau <= sor.mat[r.mid,n]){r.max <- r.mid}else{r.min <- r.mid+1}}}
    return(rbind(sor.mat[seq_len(r.min-1),,drop=FALSE],vec,sor.mat[seq_len(sor.mat.sz-r.max+1)+r.max-1,,drop=FALSE]))}

ins.sort.r <- function(vec,sor.mat,sor.mat.sz,n){ # Insert a vector into a sorted matrix (ranking along a row)
    r.min <- 1; r.max <- sor.mat.sz; tau <- vec[n]
    if(tau > sor.mat[n,r.max]){r.min <- r.max <- sor.mat.sz+1}else{
        while(r.max > r.min){
            r.mid <- floor(0.5*(r.max+r.min))
            if(tau <= sor.mat[n,r.mid]){r.max <- r.mid}else{r.min <- r.mid+1}}}
    return(cbind(sor.mat[,seq_len(r.min-1),drop=FALSE],vec,sor.mat[,seq_len(sor.mat.sz-r.max+1)+r.max-1,drop=FALSE]))}

########################################################################
########################################################################
### 1 - Resampling Steps
########################################################################
########################################################################

#############################################
#### 1.1 - Particle Normalisation
#############################################

part.control    <- function(log.p.wei){log.p.wei - max(log.p.wei)} # Control floating point precision of weights
part.normal     <- function(log.p.wei){p.wei <- exp(log.p.wei)/sum(exp(log.p.wei)); list(log.p.wei=log(p.wei),p.wei=p.wei)} # Normalise log particle weights
part.ess        <- function(log.p.wei,p.num){normalised <- part.normal(log.p.wei);list(log.p.wei=normalised$log.p.wei,p.wei=normalised$p.wei,ess=1/(p.num*sum(normalised$p.wei^2)))} # Normalise log particle weights and compute ESS

#############################################
#### 1.2 - Resampling Algorithms
#############################################
##### 1.2.1 - Multinomial Resampling
#############################################

multi.resamp     <- function(p.wei,n=length(p.wei)){ # Multinomial Resampling
    if((sum(p.wei)>0)&(n>0)){ # Check whether resampling possible
        p.idx        <- sample(1:length(p.wei),n,replace=TRUE,prob=p.wei) # Sampled Index
    }else{p.idx     <- numeric(0)} # If resampling not possible return empty vector
    list(p.idx=p.idx)} # Return particle indices

#############################################
##### 1.2.2 - Systematic Resampling
#############################################

system.resamp    <- function(p.wei,n=length(p.wei)){ # Systematic Resampling
    if((sum(p.wei)>0)&(n>0)){ # Check whether resampling possible
        cum.wei        <- cumsum(p.wei)/sum(p.wei) # Normalise and calculate cumulative weights
        samp.vec    <- seq(runif(1,0,1/n),1,1/n) # Vector of cdf samples
        p.idx        <- numeric(0); for(i in 1:n){p.idx <- c(p.idx,length(cum.wei[samp.vec[i]>=cum.wei])+1)} # Sample from cdf
    }else{p.idx     <- numeric(0)} # If resampling not possible return empty vector
    list(p.idx=p.idx)} # Return particle indices

#############################################
##### 1.2.3 - Stratified Resampling
#############################################

strat.resamp    <- function(p.wei,n=length(p.wei)){ # Stratified Resampling
    vec.length  <- length(p.wei) # Calculate the length of the input vector
    cum.wei        <- cumsum(p.wei) # Calculate cumulative weights
    if(cum.wei[vec.length]>0){ # Check whether resampling possible
        cum.wei     <- cum.wei/cum.wei[vec.length]  # Normalise cumulative weights
        samp.vec    <- seq(0,(n-1)/n,1/n) + runif(n,0,1/n) # Vector of cdf samples
        p.idx <- findInterval(samp.vec,cum.wei)+1 # Sample from cdf
    }else{p.idx     <- numeric(0)} # If resampling not possible return empty vector
    list(p.idx=p.idx)} # Return particle indices

#############################################
##### 1.2.4 - Residual Resampling
#############################################

resid.resamp    <- function(p.wei,n=length(p.wei),nest.resamp=strat.resamp){ # Residual Resampling
    if((sum(p.wei)>0)&(n>0)){ # Check whether resampling possible
        det.resamp    <- floor(n*p.wei/sum(p.wei))
        p.idx        <- c(rep(1:length(p.wei),det.resamp),nest.resamp(p.wei-det.resamp*(sum(p.wei)/n),n=n-sum(det.resamp))$p.idx)
    }else{p.idx     <- numeric(0)} # If resampling not possible return empty vector
    list(p.idx=p.idx)} # Return particle indices

########################################################################
########################################################################
#### 2 - Exact Algorithm Constructions
########################################################################
########################################################################

#############################################
#### 2.1 - Elementary EA Functions
#############################################

#############################################
##### 2.1.1 - Bessel Functionals - See Pollock et al. 2015
#############################################

eazeta        <- function(n,s,t,x,y,L,U){if(max(x-U,y-U,L-x,L-y)>=0){1}else{j<-1:(ceiling(n/2));P<--2/(t-s);D<-U-L;D1<-D*j+L;D2<-D*j-U;z<-y-x;if(even(n)){sum(exp(P*(D1-x)*(D1-y))+exp(P*(D2+x)*(D2+y))-exp(P*j^2*D^2-P*j*D*z)-exp(P*j^2*D^2+P*j*D*z))}else{sum(exp(P*(D1-x)*(D1-y))+exp(P*(D2+x)*(D2+y)))-sum(exp(P*j[1:length(j)-1]^2*D^2-P*j[1:length(j)-1]*D*z)+exp(P*j[1:length(j)-1]^2*D^2+P*j[1:length(j)-1]*D*z))}}} # Zeta Functional
eagamma        <- function(n,s,t,x,y,L,U){1-eazeta(n,s,t,x,y,L,U)} # Gamma Functional
eapsi        <- function(j,s,t,m,xoy,u){P<--2*abs(u-m)*j/(t-s);(2*abs(u-m)*j-(xoy-m))*exp(P*(abs(u-m)*j-(xoy-m)))} # Psi Functional
eachi        <- function(j,s,t,m,xoy,u){P<--2*abs(u-m)*j/(t-s);(2*abs(u-m)*j+(xoy-m))*exp(P*(abs(u-m)*j+(xoy-m)))} # Chi Functional
eadel2      <- function(n,s,t,m,xoy,u){if(max(xoy-u,m-xoy)>=0){0}else{if(even(n)){j<-1:(n/2);1-(sum(eapsi(j,s,t,m,xoy,u)-eachi(j,s,t,m,xoy,u)))/(xoy-m)}else{if(n>1){j<-1:((n-1)/2);1-(sum(eapsi(j,s,t,m,xoy,u)-eachi(j,s,t,m,xoy,u)))/(xoy-m)-eapsi(max(j)+1,s,t,m,xoy,u)/(xoy-m)}else{1-eapsi(1,s,t,m,xoy,u)/(xoy-m)}}}}
eadelR        <- function(n,s,t,x,y,m,u){if(x==m){xI<-1}else{xI<-0};if(y==m){yI<-1}else{yI<-0};if(max(xI,yI)==1){delT<-2}else{delT<-1};if(m>min(x,y)){x<--x;y<--y;m<--m;u<--u};if(max(x-u,y-u,m-x,m-y)>=0){out<-0};if(delT==1){out<-eagamma(n,s,t,x,y,m,u)/(1-exp(-2*(x-m)*(y-m)/(t-s)))};if(delT==2){if(xI*yI==0){xoy<-max(x,y);out<-eadel2(n,s,t,m,xoy,u)}else{out<-0}};if(out<0){out<-0};if(out>1){out<-1};if((t-s)==0){out <- 1}; out}
eadelC         <- function(mt,s,t,x,y,m,u){if(mt>=mpfrthr){pbn<-mt*mpfrpbn;s<-mpfr(s,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);y<-mpfr(y,precBits=pbn);m<-mpfr(m,precBits=pbn);u<-mpfr(u,precBits=pbn);c(s1=eadelR(mt,s,t,x,y,m,u),s2=eadelR(mt+1,s,t,x,y,m,u))}else{eadel_pair_cpp(mt,s,t,x,y,m,u)}}

#############################################
#### 2.1.2 - Simulate Bessel 3 Process Mid Points along with Bessel Exceedance Evaluation - See Pollock (PhD Thesis, 2013)
#############################################

scale.midpt <- function(q, s, tau, x, y, minI, bdry){
    
    ## Repeat until acceptance
    repeat{
        ### Simulation of Bessel proposal
        c.scaling     <- c(tau-s,minI*(x-y)*(tau-q)/(tau-s)^1.5,sqrt((tau-q)*(q-s))/(tau-s))
        bb.sims          <- rnorm(3,0,sd=c.scaling[3])
        w                 <- y + minI*sqrt(c.scaling[1]*((c.scaling[2]+bb.sims[1])^2+bb.sims[2]^2+bb.sims[3]^2))
        
        ### Determine acceptance / rejection of proposal
        u               <- runif(1,0,1) # Uniform RV to make accept / reject decision
        counter     <- ceiling(sqrt((tau-s)+(bdry-y)^2)/(2*abs(bdry-y))); if(even(counter)==1){counter <- counter + 1} # Determining the minimum threshold for computing the boundary
        repeat{
            bounds     <- eadelC(counter,s,q,x,w,y,bdry)*eadelC(counter,q,tau,w,y,y,bdry) # Determine crossing probability
            if(u <= bounds[1]){accI <- 1; break} # Determine whether to accept
            if(u > bounds[2]){accI <- 0; break} # Determine whether to reject
            counter <- counter + 2} # Increase counter
        
        ### If accept then break loop
        if(accI == 1){break}}
    ### Output midpt
    list(w=w)}

#############################################
#### 2.2 - Simulate first hitting times of univariate Brownian motion with symmetric boundary (Devroye 2009)
#############################################

#############################################
##### 2.2.1 - Optimal rejection sampler specification
#############################################

# t.opt <- 0.64; p.opt <- 0.57810262346829443; q.opt <- 0.422599094; rat.opt <- p.opt/(p.opt+q.opt); rat.opt <- 0.5776972

#############################################
##### 2.2.2 - Proposal and rejection sampler
#############################################

dev.pr        <- function(Jst.t=0.64,Jst.rat=0.5776972){ # First Hitting Time (J*) Proposal Function
    U <- runif(1,0,1) # Simulate uniform
    if(U < Jst.rat){ # Case 1, U <- p/(p+q)
        E <- rexp(1,1); X <- Jst.t + 8*E/(pi^2)
    }else{ # Case 2, U >= p/(p+q)
        repeat{E <- rexp(2,1); if(E[1]^2 <= 2*E[2]/Jst.t){X <- Jst.t/(1+Jst.t*E[1])^2; break}}}
    list(X=X,U=U,E=E)}

dev.rej         <- function(X,Jst.t=0.64){ # First Hitting Time (J*) Rejection Sampler given a proposed point (X)
    if(X <= Jst.t){ # Series sampler 1, used around critical point 0.64 (t*<= 0.64, see Appendix C)
        S <- exp(-1/(2*X))/2; n <- -1; Y <- runif(1,0,1)*S # Initialise alternating sequence and simulate appropriate Uniform
        repeat{
            n <- n+2 # Indexation of sequence for subsequent odd / even terms
            S <- S - (n+0.5)*exp(-2*(n+0.5)^2/X); if(Y <= S){Acc <- 1; n <- n-1; break} # Odd n
            S <- S + (n+1.5)*exp(-2*(n+1.5)^2/X); if(Y >= S){Acc <- 0; break}} # Even n
    }else{ # Series sampler 2, used around critical point 0.64 (t*> 0.64, see Appendix C)
        S <- exp(-pi^2*X/8)/2; n <- -1; Y <- runif(1,0,1)*S # Initialise alternating sequence and simulate appropriate Uniform
        repeat{
            n <- n+2 # Indexation of sequence for subsequent odd / even terms
            S <- S - (n+0.5)*exp(-(n+0.5)^2*pi^2*X/2); if(Y <= S){Acc <- 1; n <- n-1; break} # Odd n
            S <- S + (n+1.5)*exp(-(n+1.5)^2*pi^2*X/2); if(Y >= S){Acc <- 0; break}}} # Even n
    list(Acc=Acc,S=S,n=n,Y=Y,X=X)} # Output A/R. Note S will not match target density as constants have been removed (see Appendix C)

#############################################
##### 2.2.3 - First Hitting Time (J*) Simulation
#############################################

bm.pass        <- function(s=0,x=0,theta=1,Jst.t=0.64,Jst.rat=0.5776972){ # theta denotes passage level
    repeat{sim <- dev.rej(dev.pr()$X); if(sim$Acc==1){break}}
    tau <- s + theta^2*sim$X; minI <- 2*rbinom(1,1,0.5)-1; y <- x - theta*minI
    list(tau=tau,y=y,minI=minI)}

#############################################
#### 2.3 - Simulate first hitting times of multivariate Brownian motion with symmetric boundary
#############################################

#############################################
##### 2.3.1 - Initialisation of first passage times
#############################################

p.path.init     <- function(s,x,theta,dimen,Jst.t=0.64,Jst.rat=0.5776972){ # x and theta are d-dimensional vectors
    passage <- matrix(0,dimen,6); colnames(passage) <- c("dimen","tau","y","minI","l","u") # d-dimensional storage vector
    for(i in 1:dimen){ith.passage <- bm.pass(s,x[i],theta[i],Jst.t,Jst.rat); passage[i,] <- c(i,ith.passage$tau,ith.passage$y,ith.passage$minI,x[i]-theta[i],x[i]+theta[i])} # Simulate first hitting times
    list(curr.time=s,pass.mat=passage,cyc.mat=mat.sort.c(passage,dimen,"tau"))}

#############################################
##### 2.3.2 - Given a current point in space iterate to future time (with accept / reject for opposing layer and renewing the path layer)
#############################################

#### Updates the particle to a single point in space

p.path.update <- function(pass.mat,cyc.mat,curr.time,x.curr,next.time,theta,dimen){
    if(dimen>1){
        x.next <- numeric(dimen) # Storage vector for x.next
        for(i in 1:dimen){x.next[i] <- scale.midpt(q=next.time,s=curr.time,tau=pass.mat[i,"tau"],x=x.curr[i],y=pass.mat[i,"y"],minI=pass.mat[i,"minI"],bdry=(pass.mat[i,"minI"]*(pass.mat[i,"u"]-pass.mat[i,"l"])+pass.mat[i,"l"]+pass.mat[i,"u"])/2)$w} # Update path to intermediate point
        list(curr.time=next.time,next.tau=cyc.mat[1,"tau"],pass.mat=pass.mat,cyc.mat=cyc.mat,x.next=x.next)}
    else{
        x.next <- scale.midpt(q=next.time,s=curr.time,tau=pass.mat["tau"],x=x.curr,y=pass.mat["y"],minI=pass.mat["minI"],bdry=(pass.mat["minI"]*(pass.mat["u"]-pass.mat["l"])+pass.mat["l"]+pass.mat["u"])/2)$w # Update path to intermediate point
        list(curr.time=next.time,next.tau=cyc.mat["tau"],pass.mat=pass.mat,cyc.mat=cyc.mat,x.next=x.next)}}


# Updates a particle to a single point in space, corresponding with the fpt of a dimension and so also renewing the fpt accordingly

p.path.renew.fpt <- function(pass.mat,cyc.mat,curr.time,x.curr,next.time,theta,dimen){
    if(dimen>1){
        x.next                     <- p.path.update(pass.mat,cyc.mat,curr.time,x.curr,next.time,theta,dimen)$x.next # Update path to intermediate point
        dimen.update             <- cyc.mat[1,"dimen"] # Select the dimension that has to be updated
        fpt.update                  <- bm.pass(s=next.time,x=x.next[dimen.update],theta=theta[dimen.update]) # Find the next FPT for given dimension
        pass.mat[dimen.update,] <- entry.update <- c(dimen.update, fpt.update$tau, fpt.update$y, fpt.update$minI, x.next[dimen.update] - theta[dimen.update], x.next[dimen.update] + theta[dimen.update]) # Determine entry for new FPT within dimension passage matrices and update
        cyc.mat                 <- ins.sort.c(entry.update,cyc.mat[-1,,drop=FALSE],dimen-1,2) # Update tau sorted passage matrix
        list(curr.time=next.time,next.tau=cyc.mat[1,"tau"],pass.mat=pass.mat,cyc.mat=cyc.mat,x.next=x.next)}
    else{
        x.next                     <- p.path.update(pass.mat,cyc.mat,curr.time,x.curr,next.time,theta,dimen)$x.next # Update path to intermediate point
        fpt.update                 <- bm.pass(s=next.time,x=x.next,theta=theta)
        pass.mat                   <- c(1, fpt.update$tau, fpt.update$y, fpt.update$minI, x.next - theta, x.next + theta) # Determine entry for new FPT within dimension passage matrices and update
        list(curr.time=next.time,next.tau=pass.mat["tau"],pass.mat=pass.mat,cyc.mat=pass.mat,x.next=x.next)}}

## Update path selecting amongst above two appropriately

p.path.renew <- function(p.vec,x.curr,pass.mat,cyc.mat,theta,dimen){
    if(p.vec["t"]<= p.vec["tau.bar"]){fpt.event <- 0; oth.event <- 1}else{fpt.event <- 1; oth.event <- 0}
    if(fpt.event==1){path.update <-     p.path.renew.fpt(pass.mat,cyc.mat,p.vec["s"],x.curr,p.vec["next.ev"],theta,dimen)}else{path.update <- p.path.update(pass.mat,cyc.mat,p.vec["s"],x.curr,p.vec["next.ev"],theta,dimen)}
    list(curr.time=path.update$curr.time,next.tau=path.update$next.tau,pass.mat=path.update$pass.mat,cyc.mat=path.update$cyc.mat,x.next=path.update$x.next,fpt.event=fpt.event,oth.event=oth.event)}

########################################################################
########################################################################
### 3 - Scale Subfunctions - Event Rotational Scale Implementation
########################################################################
########################################################################

#############################################
#### 3.1 - Actual Phi Specification
#############################################

ac.phi          <- function(eta.pos,dapts){
    beta.pos <- un_scale_transform(eta.pos); data.evaluation <- data.eval(beta.pos,1:dsz) # Specify beta position and evaluate data
    (data.evaluation$grad.log.pi%*%t(data.evaluation$grad.log.pi) + sum(data.evaluation$lap.log.pi))/2}

#############################################
#### 3.2 - Actual Subsample Phi Specification
#############################################

##### 3.2.1 - Subsample Phi evaluation

ss.phi          <- function(eta.pos,dapts,ss.size){
    beta.pos <- un_scale_transform(eta.pos); factor <- dsz/ss.size # Specify beta position and factor
    scale.data.counter <<- scale.data.counter + ss.size*2 # Index data access counter
    data1 <- data.eval(beta.pos,dapts[1,,"dapt.1"],factor=factor); data1.star <- data.eval(beta.star,dapts[1,,"dapt.1"],factor=factor) # Evaluate data grad log and lap log (ss 1)
    data2 <- data.eval(beta.pos,dapts[1,,"dapt.2"],factor=factor); data2.star <- data.eval(beta.star,dapts[1,,"dapt.2"],factor=factor) # Evaluate data grad log and lap log (ss 2)
    ((data1$grad.log.pi - data1.star$grad.log.pi)%*%t(2*alpha.cent + data2$grad.log.pi - data2.star$grad.log.pi) + sum(data1$lap.log.pi))/2}

##### 3.2.2 - Subsampling Phi intensity bounds

ss.phiC         <- function(l.bound,u.bound){ # Function of current particle location, lower bounds (dimensionally), upper bounds (dimensionally)
    distance <- sum(pmax(abs(l.bound),abs(u.bound))^2)^(1/2) # Compute the maximal distance from the centering value
    phiL <- phiL.fn(distance) # Computation of lower bound
    phiU <- phiU.fn(distance) # Computation of upper bound
    intensity <- phiU - phiL # Compute Bounding functionals and intensity
    list(distance=distance,phiL=phiL,phiU=phiU,intensity=intensity)}

#############################################
#### 3.3 - Functionals only loaded if not otherwise defined
#############################################

.onLoad <- function(libname,pkgname){
    ### DATA ACCESS Counter
    if(exists("scale.data.counter")==FALSE){scale.data.counter <<- 0} # Set the data access counter to zero if it doesn't already exist
    
    ### TURN OFF THE SUB SAMPLER (REQUIRES SETTING ss.on <- FALSE)
    if(exists("ss.on")==TRUE){if(ss.on==FALSE){ss.phi <- function(eta.pos,dapts,ss.size){ac.phi(eta.pos,dapts)}}}
}

#############################################
#### 3.4 - Approximate Algorithm Apply Functionals
#############################################

diffn.scaling <- 3 # Approximate absolute scaling of diffuion within unit time
approx.intensity  <- function(eta.pos,time=t.inc){position.bound <- abs(eta.pos)+diffn.scaling*sqrt(time); phiCfn <- ss.phiC(-position.bound, position.bound); c(phiCfn$phiL,phiCfn$phiU)}

#############################################
#### 3.5 - Resampling Step - Event Rotational Implementation
#############################################

scale_resamp <- function(p.num,p.idx,log.p.wei,p.mat,p.layer,p.layer.sor.I,ss.size,resamp.method,ess.thresh,p.pass.arr=NULL,p.cyc.arr=NULL){
    ###### 1 ###### Perform weight control function - Subtract maximum log weight to control floating point precision, normalise and calculate ESS
    wei.control <- part.ess(log.p.wei,p.num)
    ###### 2 ###### Subtract maximum log weight to control floating point precision
    if(wei.control$ess <= ess.thresh){ # If ESS < Threshold
        ### Determine resampling indices
        r.idx <- resamp.method(wei.control$p.wei)$p.idx # Update Resampled Particle Indices
        ### Determine ancestoral lineages
        p.idx <- p.idx[r.idx] # Update particle ancestory
        ### Update p.mat (ordered by c.idx)
        p.mat <- p.mat[,r.idx,drop=FALSE] # Update particle matrix
        if(is.null(p.pass.arr)==FALSE){p.pass.arr <- p.pass.arr[,,r.idx]; dimnames(p.pass.arr)[[3]] <- sprintf("part.%i",1:p.num)} # Update p.pass.arr
        if(is.null(p.cyc.arr)==FALSE){p.cyc.arr <- p.cyc.arr[,,r.idx]; dimnames(p.cyc.arr)[[3]] <- sprintf("part.%i",1:p.num)} # Update p.cyc.arr
        ### Update p.layer (potential unordered and requiring partial resimulation
        if(p.layer.sor.I==0){p.layer <- mat.sort.r(p.layer,p.num,"c.idx")} # Sort matrix by c.idx if currently sorted by "t"
        p.layer <- p.layer[,r.idx,drop=FALSE]; p.layer["c.idx",] <- 1:p.num # Resample p.layers and and assign c.idx as 1:p.num
        p.layer["t",] <- p.layer["deg.s",] + rexp(p.num,rate=p.layer["Delta",]) # Resimulate next event times for p.layer
        if(is.null(p.pass.arr)==FALSE){p.layer["next.ev",] <- (sign(p.layer["tau.bar",]-p.layer["t",])+1)/2; p.layer["next.ev",] <- p.layer["next.ev",]*p.layer["t",] + (1-p.layer["next.ev",])*p.layer["tau.bar",]} # Update "next.ev" if this exists
        if(p.layer.sor.I==0){if(is.null(p.pass.arr)==FALSE){p.layer <- mat.sort.r(p.layer,p.num,"next.ev")}else{p.layer <- mat.sort.r(p.layer,p.num,"t")}} # Resort matrix to "t" or "next.ev" form if it has been sorted
        ### Update log weights
        log.p.wei[] <- log(1/p.num) # Update Particle Weights (even weights)
        ### Update resampling indicator
        resamp.I <- 1 # Set resampling indicator to 1
    }else{ # If ESS > Thresh
        r.idx <- 1:p.num # Retain Particle Indices
        log.p.wei <- wei.control$log.p.wei # Update Particle Weights
        resamp.I <- 0} # Set resampling indicator to 0
    list(p.idx=p.idx,r.idx=r.idx,log.p.wei=log.p.wei,p.mat=p.mat,p.layer=p.layer,resamp.I=resamp.I,ess=wei.control$ess,p.pass.arr=p.pass.arr,p.cyc.arr=p.cyc.arr)}

#############################################
#### 3.6 - Negative Weight Handling Sub-Algorithms (Approx Alg.) == Set weights to zero
#############################################

scale_zero.wei  <- function(p.event.ch,p.curr.idx,p.curr,p.anc.neg,log.p.wei,ss.phiC){
    p.anc.neg.append <- c(p.event.ch,exp(log.p.wei[p.curr.idx]),p.curr[c("t","PhiL","PhiU")]) # Ancestral Negative Weight Record
    list(p.anc.neg.append=p.anc.neg.append,ss.phiC=ss.phiC)}

########################################################################
########################################################################
### 4 - Scale Algorithm
########################################################################
########################################################################

#############################################
#### 4.1 - Scale Algorithm Approximate Version
#############################################

scale_approx <- function(p.num,t.inc,T.fin,ss.phi,ss.phiC,dimen,scale_transform,un_scale_transform,T.start=0,x.init=NULL,x.init.random=TRUE,ss.size=1,ess.thresh=0.5,resamp.method=resid.resamp,neg.wei.mech=scale_zero.wei,prev.simn=NULL,progress.check=FALSE,phi.record=FALSE,resamp.freq=5*p.num,theta=NULL,p.path.renew=p.path.renew){
    #################
    #### (0)1 #### Initalise Algorithm
    #################
    ##### (0)1.1 #### Define dimension and other variables not pre-defined in input
    p.mu <- matrix(c(rep(0,dimen)),ncol=1,dimnames=list(sprintf("dim.%i",1:dimen),"mean")); p.covar <- diag(dimen); dimnames(p.covar) <- list(sprintf("dim.%i",1:dimen),sprintf("dim.%i",1:dimen)) # Specify algorithm dimension and Normal covariance matrix
    ##### (0)1.2 Initialise particle storage framework
    if(is.null(prev.simn)==TRUE){ #i.e. starting algorithm from 0
        ### Algorithm initialisation point
        if(is.null(x.init)==TRUE){ # Determine if user has specified particle initialisation
            if(x.init.random==TRUE){ # Check if particles have to be randomly initialised
                if(exists("precon")==FALSE){p.mat <- matrix(rnorm(dimen*p.num,0,1),dimen,p.num,dimnames=list(sprintf("dim.%i",1:dimen),NULL))}else{p.mat <- matrix(apply(t(mvrnorm(n=p.num,mu=beta.star,Sigma=as.matrix(precon))),2,scale_transform),dimen,p.num,dimnames=list(sprintf("dim.%i",1:dimen),NULL))} # Initialise from N(0,1) marginals in transformed space if there exists no preconditioning matrix, otherwise simulate from preconditioning matrix and transform.
            }else{p.mat <- matrix(0,dimen,p.num,dimnames=list(sprintf("dim.%i",1:dimen),NULL))}
        }else{
            if(is.vector(x.init)){x.init <- matrix(x.init,length(x.init),1)} # If x.init is a vector transform to matrix
            if(dim(x.init)[2]==1){x.init <- matrix(rep(x.init,p.num),dimen,p.num)} # If x.init is a point then repeat
            p.mat <- apply(x.init,2,scale_transform); dimnames(p.mat)<-list(sprintf("dim.%i",1:dimen),NULL) # Transform x.init and define as p.mat
        } # Setting of p.mat depending on x.init NULL, d dimn vector, dxp dimn vector, 1xd matrix or pxd matrix
        ### Define storage matrices
        p.idx <- 1:p.num; log.p.wei <- rep(log(1/p.num),p.num) # Set idxs and weights
        p.layer <- matrix(0,7,p.num,dimnames=list(c("c.idx","s","deg.s","t","PhiL","PhiU","Delta"),NULL)) # Establish particle layer matrix
        p.layer["c.idx",] <- p.idx # Set particle idxs in layer matrix
    }else{ #i.e. restarting algorithm from some point
        T.start <- prev.simn$T.fin; p.idx <- prev.simn$p.idx; log.p.wei <- prev.simn$log.p.wei; p.mat <- prev.simn$p.mat; p.layer <- prev.simn$p.layer} # Restablish Algorithm
    ##### (0)1.3 Initialise particle system
    p.layer[c("PhiL","PhiU"),] <- apply(p.mat,2,approx.intensity); p.layer[c("Delta"),] <- p.layer[c("PhiU"),] - p.layer[c("PhiL"),] # Set p.layer current particle intensity
    p.bet.degrad.r <- p.bet.degrad <- matrix(0,nrow=2,ncol=0,dimnames=list(c("c.idx","degrad"))) # Initialise between resampling log weight degration
    p.layer["s",] <- p.layer["deg.s",] <- T.start # Initialise algorithms start point
    p.layer["t",] <- p.layer["deg.s",] + rexp(p.num,rate=p.layer["Delta",]) # Set p.layer next event times
    ##### (0)1.4 Initialise particle ancestral storage framework
    p.anc.idx <- matrix(p.idx,1,p.num,dimnames=list(NULL,sprintf("pos.%i",1:p.num))); p.anc.wei <- matrix(part.normal(log.p.wei)$p.wei,1,p.num,dimnames=list(NULL,sprintf("part.%i",1:p.num)));  p.anc.mat <- array(p.mat,c(dimen,p.num,1),dimnames=list(sprintf("dim.%i",1:dimen),sprintf("part.%i",1:p.num),NULL)) # Ancestral recording
    p.phi.hist <- matrix(0,0,dimen+1,dimnames=list(NULL,c("ss.phi",sprintf("dim.%i",1:dimen)))); p.anc.resamp <- numeric(0); p.anc.ess <- numeric(0); p.anc.neg <- matrix(0,0,5,dimnames=list(NULL,c("event.ch","wt-","t","PhiL","PhiU"))) # Diagnostic ancestral recording
    ##### (0)1.5 Sort p.layer matrix
    p.layer <- mat.sort.r(p.layer,p.num,"t") # Sort p.layer matrix by next event time
    ######################
    #### (t.inc) #### LOOP
    ######################
    ##### (t.inc)2 #### Define loop stopping criteria
    ######################
    if(exists("scale.data.counter")==FALSE){scale.data.counter <<- 0} # Initialise data access counter if it doesnt exist
    anc.times <- curr.deg.time <- T.start; t.inc.next <- T.start + t.inc; resamp.counter <- 1 # Initialise current degradation time, time of first increment and resampling counter
    while(curr.deg.time < T.fin){ # Iterate until final time is reached
        ##### (t.inc)2.1 ##### Increment particle set to next event time
        ######################
        if(p.layer["t",1] >= t.inc.next){ # If we have reached the next increment point
            ###### (t.inc)2.2.a.1 ##### Compute particles at intermediate time point
            p.layer <- mat.sort.r(p.layer,p.num,"c.idx") # Sort p.layer matrix by current particle index
            p.mat <- p.mat + t(sqrt(t.inc.next-p.layer["s",])*mvrnorm(p.num,p.mu,p.covar)) # Update particle locations
            ###### (t.inc)2.2.a.2 ##### Update particle weights
            p.bet.degrad <- cbind(p.bet.degrad,rbind(p.layer[c("c.idx"),],(t.inc.next-p.layer["deg.s",])*p.layer["PhiL",])) # Find full particle between resampling degradation matrix
            p.wei.agg <- aggregate(x=p.bet.degrad["degrad",],by=list(p.bet.degrad["c.idx",]),FUN="sum") # Aggregate weight degradation
            log.p.wei <- log.p.wei - p.wei.agg[,"x"] # Update weight of particles due to weight degradation to mesh point
            ###### (t.inc)2.2.a.3 ##### Append increment to ancesotral recording
            anc.times <- c(anc.times,t.inc.next) # Record time of resampling
            p.anc.idx <- rbind(p.anc.idx,p.idx); p.idx <- 1:p.num # Append indices and resent current indices
            p.normalise <- part.normal(log.p.wei); log.p.wei <- p.normalise$log.p.wei; p.anc.wei <- rbind(p.anc.wei,p.normalise$p.wei) # Control, normalise log weights and append normalised weights
            p.anc.mat <- abind(p.anc.mat,array(p.mat,c(dimen,p.num,1)),along=3) # Append increment to ancestoral recording
            ###### (t.inc)2.2.a.4 ##### Update particle system
            curr.deg.time <- t.inc.next # Index left hand time point and current time
            p.layer["deg.s",] <- p.layer["s",] <- curr.deg.time  # Index left hand time point and current time
            t.inc.next <- min(curr.deg.time + t.inc, T.fin) # Index increment time
            ###### (t.inc)2.2.a.5 ##### Resample
            p.resamp <- scale_resamp(p.num,p.idx,log.p.wei,p.mat,p.layer,p.layer.sor.I=1,ss.size,resamp.method,ess.thresh) # Resample particles
            p.anc.ess <- c(p.anc.ess,p.resamp$ess) # Record ancestoral ESS at mesh point
            p.idx <- p.resamp$p.idx; log.p.wei <- p.resamp$log.p.wei # Update particle indicies and weights (normalised)
            if(p.resamp$resamp.I==1){p.anc.resamp <- c(p.anc.resamp,curr.deg.time); p.mat <- p.resamp$p.mat; p.layer <- p.resamp$p.layer} # If there has been resampling, update ancestral resampling times, layer
            p.layer <- mat.sort.r(p.layer,p.num,"t") # Resort p.layer by time
            resamp.counter <- 1 # Reset resample counter
            if(progress.check==TRUE){print(curr.deg.time)} # Output current time if checking progress
        }else{ # If we haven't reached the next increment point
            ###### (t.inc)2.2.b.1 ##### Select particle to be updated
            p.curr <- p.layer[,1] # Next particle to be updated
            p.curr.idx <- p.curr["c.idx"] # Define current index
            p.curr.dapts <- array(sample(dsz,2*ss.size,replace=TRUE),c(1,ss.size,2),dimnames=list("",sprintf("ss.%i",1:ss.size),sprintf("dapt.%i",1:2))) # Find particle subsampled data points
            ###### (t.inc)2.2.b.2 ##### Update weight
            p.loc.next <- p.mat[,p.curr.idx] + sqrt(p.curr["t"]-p.curr["s"])*mvrnorm(1,p.mu,p.covar) # Particle location at t
            p.phi.eval <- ss.phi(p.loc.next,p.curr.dapts,ss.size)
            p.event.ch <- (p.curr["PhiU"] - p.phi.eval)/p.curr["Delta"]; if(phi.record==TRUE){p.phi.hist <- rbind(p.phi.hist,c(p.phi.eval,p.loc.next))} # Raw Poisson event increment
            if(p.event.ch < 0){p.neg <- neg.wei.mech(p.event.ch,p.curr.idx,p.curr,p.anc.neg,log.p.wei,ss.phiC); log.p.wei[p.curr.idx] <- -Inf; p.anc.neg <- rbind(p.anc.neg,p.neg$p.anc.neg.append); ss.phiC <- p.neg$ss.phiC}else{log.p.wei[p.curr.idx] <- log.p.wei[p.curr.idx] + log(p.event.ch)} # Apply chosen negative weight mechanism
            p.bet.degrad <- cbind(p.bet.degrad,c(p.curr.idx,(p.curr["t"]-p.curr["deg.s"])*p.curr[c("PhiL")])) # Record individual particle degradation (since last resampling point)
            ###### (t.inc)2.2.b.3 ##### Update particle to new point
            p.mat[,p.curr.idx] <- p.loc.next # Update particle location
            p.phi.bds <- approx.intensity(p.loc.next,t.inc); p.curr[c("PhiL","PhiU")] <- c(p.phi.bds[1],p.phi.bds[2]); p.curr["Delta"] <- p.curr["PhiU"] - p.curr["PhiL"] # Update phi bounds and intensity
            p.curr["s"] <- p.curr["deg.s"] <- p.curr["t"] # Update current particle degradation time
            p.curr["t"] <- p.curr["deg.s"] + rexp(1,rate=p.curr["Delta"]) # Update current particle event time
            p.layer <- ins.sort.r(p.curr,p.layer[,-1,drop=FALSE],p.num-1,"t") # Re-insert particle into p.layer
            ###### (t.inc)2.2.b.4 ##### Resample
            if(resamp.counter%%(resamp.freq)==0){ # Resample if resample counter sufficient
                ### Update particle degradation matrix and aggregate weight degradation
                p.bet.degrad <- cbind(p.bet.degrad,rbind(p.layer[c("c.idx"),],(p.curr["deg.s"]-p.layer["deg.s",])*p.layer["PhiL",])) # Find full particle between resampling degradation matrix
                p.wei.agg <- aggregate(x=p.bet.degrad["degrad",],by=list(p.bet.degrad["c.idx",]),FUN="sum") # Aggregate weight degradation
                log.p.wei <- log.p.wei - p.wei.agg[,"x"] # Update weight of particle set due to weight degradation to event time
                ### Update current degradation time, update p.layer and resent particle degradation matrix
                p.layer["deg.s",] <- curr.deg.time <- p.curr["deg.s"] # Index current degradation time
                p.bet.degrad <- p.bet.degrad.r  # Reset particle between resampling degradation matrix
                ### Resample
                p.resamp <- scale_resamp(p.num,p.idx,log.p.wei,p.mat,p.layer,p.layer.sor.I=0,ss.size,resamp.method,ess.thresh) # Resample particles
                p.idx <- p.resamp$p.idx; log.p.wei <- p.resamp$log.p.wei # Update particle indicies and weights (normalised)
                if(p.resamp$resamp.I==1){p.anc.resamp <- c(p.anc.resamp,curr.deg.time); p.mat <- p.resamp$p.mat; p.layer <- p.resamp$p.layer} # If there has been resampling, update ancestral resampling times, layer and dapts
            } # Close reasampling if statement
            ###### (t.inc)2.2.b.4 ##### Resample Counter Index
            resamp.counter <- resamp.counter + 1 # Index resample counter
            ###### (t.inc)2.2.c ##### Close loop
        } # Close particle set increment loop
        ###### (t.inc)2.3 ##### Close loop
        #######################
    } # Close time increment loop
    ########################
    #### (T.fin) #### Output
    ########################
    list(p.num=p.num,p.idx=p.idx,log.p.wei=log.p.wei,p.mat=p.mat,p.layer=p.layer,theta=theta,p.anc.idx=p.anc.idx,p.anc.wei=p.anc.wei,p.anc.mat=p.anc.mat,resamp.freq=resamp.freq,resamp.method=resamp.method,neg.wei.mech=neg.wei.mech,p.anc.resamp=p.anc.resamp,p.anc.ess=p.anc.ess,p.anc.neg=p.anc.neg,p.phi.hist=p.phi.hist,anc.times=anc.times,T.start=T.start,t.inc=t.inc,T.fin=T.fin,dimen=dimen,p.mu=p.mu,ss.size=ss.size,ess.thresh=ess.thresh,ss.phi=ss.phi,ss.phiC=ss.phiC,p.cyc.arr=NULL,p.pass.arr=NULL,scale_transform=scale_transform,un_scale_transform=un_scale_transform,progress.check=progress.check,phi.record=phi.record,p.path.renew=p.path.renew)} # Output list

#############################################
#### 4.2 - Scale Algorithm Exact Version
#############################################

scale_exact <- function(p.num,t.inc,T.fin,ss.phi,ss.phiC,dimen,scale_transform,un_scale_transform,T.start=0,x.init=NULL,x.init.random=TRUE,ss.size=1,ess.thresh=0.5,resamp.method=resid.resamp,neg.wei.mech=scale_zero.wei,prev.simn=NULL,progress.check=FALSE,phi.record=FALSE,resamp.freq=5*p.num,theta=NULL,p.path.renew=p.path.renew){
    #################
    #### (0)1 #### Initalise Algorithm
    #################
    ##### (0)1.1 #### Define dimension and other variables not pre-defined in input
    p.mu <- matrix(c(rep(0,dimen)),ncol=1,dimnames=list(sprintf("dim.%i",1:dimen),"mean")); p.covar <- diag(dimen); dimnames(p.covar) <- list(sprintf("dim.%i",1:dimen),sprintf("dim.%i",1:dimen)) # Specify algorithm dimension and Normal covariance matrix
    if(is.null(theta)){theta <- rep(dimen^(1/4)*sqrt(t.inc),dimen)}
    ##### (0)1.2 Initialise particle storage framework
    if(is.null(prev.simn)==TRUE){ #i.e. starting algorithm from 0
        ### Algorithm initialisation point
        if(is.null(x.init)==TRUE){ # Determine if user has specified particle initialisation
            if(x.init.random==TRUE){ # Check if particles have to be randomly initialised
                if(exists("precon")==FALSE){p.mat <- matrix(rnorm(dimen*p.num,0,1),dimen,p.num,dimnames=list(sprintf("dim.%i",1:dimen),NULL))}else{p.mat <- matrix(apply(t(mvrnorm(n=p.num,mu=beta.star,Sigma=as.matrix(precon))),2,scale_transform),dimen,p.num,dimnames=list(sprintf("dim.%i",1:dimen),NULL))} # Initialise from N(0,1) marginals in transformed space if there exists no preconditioning matrix, otherwise simulate from preconditioning matrix and transform.
            }else{p.mat <- matrix(0,dimen,p.num,dimnames=list(sprintf("dim.%i",1:dimen),NULL))} # If simulation is not at random, then initialise at origin.
        }else{
            if(is.vector(x.init)){x.init <- matrix(x.init,length(x.init),1)} # If x.init is a vector transform to matrix
            if(dim(x.init)[2]==1){x.init <- matrix(rep(x.init,p.num),dimen,p.num)} # If x.init is a point then repeat
            p.mat <- apply(x.init,2,scale_transform); dimnames(p.mat)<-list(sprintf("dim.%i",1:dimen),NULL) # Transform x.init and define as p.mat
        } # Setting of p.mat depending on x.init NULL, d dimn vector, dxp dimn vector, 1xd matrix or pxd matrix
        ### Define storage matrices
        p.idx <- 1:p.num; log.p.wei <- rep(log(1/p.num),p.num) # Set idxs and weights
        p.layer <- matrix(0,9,p.num,dimnames=list(c("c.idx","s","deg.s","t","tau.bar","next.ev","PhiL","PhiU","Delta"),NULL)) # Establish particle layer matrix
        p.layer["c.idx",] <- p.idx # Set particle idxs in layer matrix
        p.pass.arr <- p.cyc.arr <- array(0,c(dimen,6,p.num),dimnames=list(NULL,c("dimen","tau","y","minI","l","u"),sprintf("part.%i",1:p.num))) # Define particle passage and cyclic arrays
        for(i in 1:p.num){p.taus <- p.path.init(s=T.start,x=p.mat[,i],theta=theta,dimen=dimen); p.pass.arr[,,i] <- p.taus$pass.mat; p.cyc.arr[,,i] <- p.taus$cyc.mat} # Initialise particle first passage times
        p.layer["tau.bar",] <- p.cyc.arr[1,"tau",] # Determine each particles tau.bar
    }else{ #i.e. restarting algorithm from some point
        T.start <- prev.simn$T.fin; p.idx <- prev.simn$p.idx; log.p.wei <- prev.simn$log.p.wei; p.mat <- prev.simn$p.mat; p.layer <- prev.simn$p.layer; p.pass.arr <- prev.simn$p.pass.arr; p.cyc.arr <- prev.simn$p.cyc.arr} # Restablish Algorithm
    ##### (0)1.3 Initialise particle system
    for(i in 1:p.num){p.phi.bds <- ss.phiC(p.pass.arr[,"l",i],p.pass.arr[,"u",i]); p.layer[c("PhiL","PhiU"),i]<-c(p.phi.bds$phiL,p.phi.bds$phiU)}; p.layer[c("Delta"),] <- p.layer["PhiU",] - p.layer["PhiL",] # Set p.layer current particle intensity
    p.bet.degrad.r <- p.bet.degrad <- numeric(p.num) # Initialise between resampling log weight degration
    p.layer["s",] <- p.layer["deg.s",] <- T.start # Initialise algorithms start point
    p.layer["t",] <- p.layer["deg.s",] + rexp(p.num,rate=p.layer["Delta",]) # Set p.layer next event times
    p.layer["next.ev",] <- (sign(p.layer["tau.bar",]-p.layer["t",])+1)/2; p.layer["next.ev",] <- p.layer["next.ev",]*p.layer["t",] + (1-p.layer["next.ev",])*p.layer["tau.bar",] # For each particle determine whether a tau.bar or poisson event time occurs next and at what time
    ##### (0)1.4 Initialise particle ancestral storage framework
    p.anc.idx <- matrix(p.idx,1,p.num,dimnames=list(NULL,sprintf("pos.%i",1:p.num))); p.anc.wei <- matrix(part.normal(log.p.wei)$p.wei,1,p.num,dimnames=list(NULL,sprintf("part.%i",1:p.num)));  p.anc.mat <- array(p.mat,c(dimen,p.num,1),dimnames=list(sprintf("dim.%i",1:dimen),sprintf("part.%i",1:p.num),NULL)) # Ancestral recording
    p.phi.hist <- matrix(0,0,dimen+1,dimnames=list(NULL,c("ss.phi",sprintf("dim.%i",1:dimen)))); p.anc.resamp <- numeric(0); p.anc.ess <- numeric(0); p.anc.neg <- matrix(0,0,5,dimnames=list(NULL,c("event.ch","wt-","t","PhiL","PhiU"))) # Diagnostic ancestral recording
    ##### (0)1.5 Sort p.layer matrix
    p.layer <- mat.sort.r(p.layer,p.num,"next.ev") # Sort p.layer matrix by next event time
    ######################
    #### (t.inc) #### LOOP
    ######################
    ##### (t.inc)2 #### Define loop stopping criteria
    ######################
    if(exists("scale.data.counter")==FALSE){scale.data.counter <<- 0} # Initialise data access counter if it doesnt exist
    anc.times <- curr.deg.time <- T.start; t.inc.next <- T.start + t.inc; resamp.counter <- 1 # Initialise current degradation time, time of first increment and resampling counter
    while(curr.deg.time < T.fin){ # Iterate until final time is reached
        ##### (t.inc)2.1 ##### Increment particle set to next event time
        ######################
        if(p.layer["next.ev",1] >= t.inc.next){ # If we have reached the next increment point
            ###### (t.inc)2.2.a.1 ##### Compute particles at intermediate time point
            p.layer <- mat.sort.r(p.layer,p.num,"c.idx") # Sort p.layer matrix by current particle index
            for(loc.up in 1:p.num){p.mat[,loc.up]<-p.path.update(pass.mat=p.pass.arr[,,loc.up],cyc.mat=p.cyc.arr[,,loc.up],curr.time=p.layer["s",loc.up],x.curr=p.mat[,loc.up],next.time=t.inc.next,theta,dimen)$x.next}
            ###### (t.inc)2.2.a.2 ##### Update particle weights
            p.bet.degrad[p.layer[c("c.idx"),]] <- p.bet.degrad[p.layer[c("c.idx"),]] + (p.curr["deg.s"]-p.layer["deg.s",])*p.layer["PhiL",] # Find full particle between resampling degradation matrix
            log.p.wei <- log.p.wei - p.bet.degrad # Update weight of particles due to weight degradation to mesh point
            ###### (t.inc)2.2.a.3 ##### Append increment to ancesotral recording
            anc.times <- c(anc.times,t.inc.next) # Record time of resampling
            p.anc.idx <- rbind(p.anc.idx,p.idx); p.idx <- 1:p.num # Append indices and resent current indices
            p.normalise <- part.normal(log.p.wei); log.p.wei <- p.normalise$log.p.wei; p.anc.wei <- rbind(p.anc.wei,p.normalise$p.wei) # Control, normalise log weights and append normalised weights
            p.anc.mat <- abind(p.anc.mat,array(p.mat,c(dimen,p.num,1)),along=3) # Append increment to ancestoral recording
            ###### (t.inc)2.2.a.4 ##### Update particle system
            curr.deg.time <- t.inc.next # Index left hand time point and current time
            p.layer["deg.s",] <- p.layer["s",] <- curr.deg.time  # Index left hand time point and current time
            t.inc.next <- min(curr.deg.time + t.inc, T.fin) # Index increment time
            ###### (t.inc)2.2.a.5 ##### Resample
            p.resamp <- scale_resamp(p.num,p.idx,log.p.wei,p.mat,p.layer,p.layer.sor.I=1,ss.size,resamp.method,ess.thresh,p.pass.arr=p.pass.arr,p.cyc.arr=p.cyc.arr) # Resample particles
            p.anc.ess <- c(p.anc.ess,p.resamp$ess) # Record ancestoral ESS at mesh point
            p.idx <- p.resamp$p.idx; log.p.wei <- p.resamp$log.p.wei # Update particle indicies and weights (normalised)
            if(p.resamp$resamp.I==1){p.anc.resamp <- c(p.anc.resamp,curr.deg.time); p.mat <- p.resamp$p.mat; p.layer <- p.resamp$p.layer; p.pass.arr <- p.resamp$p.pass.arr; p.cyc.arr <- p.resamp$p.cyc.arr} # If there has been resampling, update ancestral resampling times, layer, dapts and cycling arrays
            p.layer <- mat.sort.r(p.layer,p.num,"next.ev") # Redefine storage vectors II
            resamp.counter <- 1 # Reset resample counter
            if(progress.check==TRUE){print(curr.deg.time)} # Output current time if checking progress
        }else{ # If we haven't reached the next increment point
            ###### (t.inc)2.2.b.1 ##### Select particle to be updated
            p.curr <- p.layer[,1] # Next particle to be updated
            p.curr.idx <- p.curr["c.idx"] # Define current index
            p.curr.dapts <- array(sample(dsz,2*ss.size,replace=TRUE),c(1,ss.size,2),dimnames=list("",sprintf("ss.%i",1:ss.size),sprintf("dapt.%i",1:2))) # Find particle subsampled data points
            ###### (t.inc)2.2.b.2 ##### Update particle trajectory
            traj.update <- p.path.renew(p.vec=p.curr,x.curr=p.mat[,p.curr.idx],pass.mat=p.pass.arr[,,p.curr.idx],cyc.mat=p.cyc.arr[,,p.curr.idx],theta,dimen) # Iterate passage matrix
            p.mat[,p.curr.idx] <- p.loc.next <- traj.update$x.next # Particle location at next.ev
            ###### (t.inc)2.2.b.3 ##### Update particle degradation matrix
            p.bet.degrad[p.curr.idx] <- p.bet.degrad[p.curr.idx] + (p.curr["next.ev"]-p.curr["deg.s"])*p.curr[c("PhiL")] # Record individual particle degradation (since last resampling point)
            ###### (t.inc)2.2.b.4 ##### If particle to be updated has an event before the next tau.bar, then update to Poisson event time
            if(traj.update$oth.event==1){
                ###### (t.inc)2.2.b.4.2 ##### Update Poisson weight
                p.phi.eval <- ss.phi(p.loc.next,p.curr.dapts,ss.size)
                p.event.ch <- (p.curr["PhiU"] - p.phi.eval)/p.curr["Delta"]; if(phi.record==TRUE){p.phi.hist <- rbind(p.phi.hist,c(p.phi.eval,p.loc.next))} # Raw Poisson event increment
                if(p.event.ch < 0){p.neg <- neg.wei.mech(p.event.ch,p.curr.idx,p.curr,p.anc.neg,log.p.wei,ss.phiC); log.p.wei[p.curr.idx] <- -Inf; p.anc.neg <- rbind(p.anc.neg,p.neg$p.anc.neg.append); ss.phiC <- p.neg$ss.phiC}else{log.p.wei[p.curr.idx] <- log.p.wei[p.curr.idx] + log(p.event.ch)} # Apply chosen negative weight mechanism
            }else{
                ###### (t.inc)2.2.b.5 ##### If particle to be updated has no event before the next tau.bar, then update passage trajectory
                p.pass.arr[,,p.curr.idx] <- traj.update$pass.mat; p.cyc.arr[,,p.curr.idx] <- traj.update$cyc.mat # Update pass.mat and cyc.mat
                p.curr["tau.bar"] <- traj.update$next.tau # Determine next particle passage time
                ###### (t.inc)2.2.b.6 ##### Update particle Poisson intensity
                p.phi.bds <- ss.phiC(p.pass.arr[,"l",p.curr.idx],p.pass.arr[,"u",p.curr.idx]); p.curr[c("PhiL","PhiU")] <- c(p.phi.bds$phiL,p.phi.bds$phiU); p.curr["Delta"] <- p.curr["PhiU"] - p.curr["PhiL"] # Update phi bounds and intensity
            }
            ###### (t.inc)2.2.b.7 ##### Update particle event times and layer set
            p.curr["s"] <- p.curr["deg.s"] <- p.curr["next.ev"] # Update current particle degradation time
            p.curr["t"] <- p.curr["next.ev"] + rexp(1,rate=p.curr["Delta"]) # Update current particle Poisson event time
            p.curr["next.ev"] <- min(p.curr["tau.bar"],p.curr["t"]) # Determine next particle event time
            p.layer <- ins.sort.r(p.curr,p.layer[,-1,drop=FALSE],p.num-1,"next.ev") # Re-insert particle into p.layer
            ###### (t.inc)2.2.b.8 ##### Resample
            if(resamp.counter%%(resamp.freq)==0){ # Resample if resample counter sufficient
                ### Update particle degradation matrix and aggregate weight degradation
                p.bet.degrad[p.layer[c("c.idx"),]] <- p.bet.degrad[p.layer[c("c.idx"),]] + (p.curr["deg.s"]-p.layer["deg.s",])*p.layer["PhiL",] # Find full particle between resampling degradation matrix
                log.p.wei <- log.p.wei - p.bet.degrad # Update weight of particle set due to weight degradation to event time
                ### Update current degradation time, update p.layer and resent particle degradation matrix
                p.layer["deg.s",] <- curr.deg.time <- p.curr["deg.s"] # Index current degradation time
                p.bet.degrad <- p.bet.degrad.r  # Reset particle between resampling degradation matrix
                ### Resample
                p.resamp <- scale_resamp(p.num,p.idx,log.p.wei,p.mat,p.layer,p.layer.sor.I=0,ss.size,resamp.method,ess.thresh,p.pass.arr=p.pass.arr,p.cyc.arr=p.cyc.arr) # Resample particles
                p.idx <- p.resamp$p.idx; log.p.wei <- p.resamp$log.p.wei # Update particle indicies and weights (normalised)
                if(p.resamp$resamp.I==1){p.anc.resamp <- c(p.anc.resamp,curr.deg.time); p.mat <- p.resamp$p.mat; p.layer <- p.resamp$p.layer; p.pass.arr <- p.resamp$p.pass.arr; p.cyc.arr <- p.resamp$p.cyc.arr} # If there has been resampling, update ancestral resampling times, layer
            } # Close reasampling if statement
            ###### (t.inc)2.2.b.4 ##### Resample Counter Index
            resamp.counter <- resamp.counter + 1 # Index resample counter
            ###### (t.inc)2.2.c ##### Close loop
        } # Close particle set increment loop
        ###### (t.inc)2.3 ##### Close loop
        #######################
    } # Close time increment loop
    ########################
    #### (T.fin) #### Output
    ########################
    list(p.num=p.num,p.idx=p.idx,log.p.wei=log.p.wei,p.mat=p.mat,p.layer=p.layer,theta=theta,p.anc.idx=p.anc.idx,p.anc.wei=p.anc.wei,p.anc.mat=p.anc.mat,resamp.freq=resamp.freq,resamp.method=resamp.method,neg.wei.mech=neg.wei.mech,p.anc.resamp=p.anc.resamp,p.anc.ess=p.anc.ess,p.anc.neg=p.anc.neg,p.phi.hist=p.phi.hist,anc.times=anc.times,T.start=T.start,t.inc=t.inc,T.fin=T.fin,dimen=dimen,p.mu=p.mu,ss.size=ss.size,ess.thresh=ess.thresh,ss.phi=ss.phi,ss.phiC=ss.phiC,p.cyc.arr=p.cyc.arr,p.pass.arr=p.pass.arr,scale_transform=scale_transform,un_scale_transform=un_scale_transform,progress.check=progress.check,phi.record=phi.record,p.path.renew=p.path.renew)} # Output list

########################################################################
########################################################################
#### 5 - Extend existing scale algorithm
########################################################################
########################################################################

scale_extend <- function(prev.simn,t.inc,T.extend,scale_method=scale_exact){
    ### Extend current simulation
    new.simn <- scale_method(p.num=prev.simn$p.num,t.inc=t.inc,T.fin=T.extend+prev.simn$T.fin,ss.phi=ss.phi,ss.phiC=prev.simn$ss.phiC,dimen=prev.simn$dimen,scale_transform=scale_transform,un_scale_transform=un_scale_transform,T.start=prev.simn$T.fin,ss.size=prev.simn$ss.size,ess.thresh=prev.simn$ess.thresh,resamp.method=prev.simn$resamp.method,neg.wei.mech=prev.simn$neg.wei.mech,prev.simn=prev.simn,progress.check=prev.simn$progress.check,phi.record=prev.simn$phi.record,resamp.freq=prev.simn$resamp.freq,theta=prev.simn$theta,p.path.renew=prev.simn$p.path.renew)
    ### Append simulation
    append.len          <- length(new.simn$anc.times)-1 # Compute the number of additional mesh points to be added and
    if(append.len>=1){ # append if there are one or more
        append.vec      <- 2:length(new.simn$anc.times)
        p.anc.idx       <- rbind(prev.simn$p.anc.idx,new.simn$p.anc.idx[append.vec,,drop=FALSE]) # Append new ancestral particle indices
        p.anc.wei       <- rbind(prev.simn$p.anc.wei,new.simn$p.anc.wei[append.vec,,drop=FALSE]) # Append new ancestral weights
        p.anc.mat       <- abind(prev.simn$p.anc.mat,new.simn$p.anc.mat[,,append.vec,drop=FALSE],along=3) # Append new ancestral particle matrix
        p.anc.resamp    <- c(prev.simn$p.anc.resamp,new.simn$p.anc.resamp) # Append new ancestral resampling times
        p.anc.ess       <- c(prev.simn$p.anc.ess,new.simn$p.anc.ess) # Append new ancestral resampling times
        p.phi.hist      <- rbind(prev.simn$p.phi.hist,new.simn$p.phi.hist) # Append new phi values
        p.anc.neg       <- rbind(prev.simn$p.anc.neg,new.simn$p.anc.neg) # Append new negative weight time details and count number of negative weights
        anc.times       <- c(prev.simn$anc.times,new.simn$anc.times[append.vec]) # Append new ancestral mesh times
    }else{ # Case where no simulation takes place - output prev.simn information
        p.anc.idx <- prev.simn$p.anc.idx;p.anc.wei <- prev.simn$p.anc.wei; p.anc.mat <- prev.simn$p.anc.mat; p.anc.resamp <- prev.simn$p.anc.resamp; p.anc.ess <- prev.simn$p.anc.ess; p.phi.hist <- prev.simn$p.phi.hist; p.anc.neg <- prev.simn$p.anc.neg; anc.times <- prev.simn$anc.times}
    list(p.num=new.simn$p.num,p.idx=new.simn$p.idx,log.p.wei=new.simn$log.p.wei,p.mat=new.simn$p.mat,p.layer=new.simn$p.layer,theta=new.simn$theta,p.anc.idx=p.anc.idx,p.anc.wei=p.anc.wei,p.anc.mat=p.anc.mat,resamp.freq=new.simn$resamp.freq,resamp.method=new.simn$resamp.method,neg.wei.mech=new.simn$neg.wei.mech,p.anc.resamp=p.anc.resamp,p.anc.ess=p.anc.ess,p.anc.neg=p.anc.neg,p.phi.hist=p.phi.hist,anc.times=anc.times,T.start=anc.times[1],t.inc=new.simn$t.inc,T.fin=new.simn$T.fin,dimen=new.simn$dimen,p.mu=new.simn$p.mu,ss.size=new.simn$ss.size,ess.thresh=new.simn$ess.thresh,ss.phi=new.simn$ss.phi,ss.phiC=new.simn$ss.phiC,p.cyc.arr=new.simn$p.cyc.arr,p.pass.arr=new.simn$p.pass.arr,scale_transform=prev.simn$scale_transform,un_scale_transform=prev.simn$un_scale_transform,progress.check=new.simn$progress.check,phi.record=new.simn$phi.record,p.path.renew=new.simn$p.path.renew)} # Output list

########################################################################
########################################################################
#### 6 - Scale Ergodic Extraction
########################################################################
########################################################################

# Extract from simulation particle locations and weights appropriately normalised and according to user specified resolution
scale_ergodic <- function(simn,retain.frac=NULL,retain.frac.range=NULL,trunc.time=NULL,trunc.time.range=NULL,even.weights=TRUE,untransform=FALSE,retain.every=1){
    #### 1.1 #### Assign frac.idxs depending on user assignation of retain frac of simulation
    if(is.null(retain.frac)==TRUE){
        frac.idxs <- 1:length(simn$anc.times) # If no specification include all
    }else{
        frac.idxs <- floor(length(simn$anc.times)*(1-retain.frac)):length(simn$anc.times)} # If specified include specified number
    if(is.null(retain.frac.range)==FALSE){
        frac.idxs <- floor(length(simn$anc.times)*retain.frac.range[1]):ceiling(length(simn$anc.times)*retain.frac.range[2])} # If a range of idxs to be included has been specified then choose specified range
    #### 1.2 #### Assign time.idxs depending on user assignation of trunc.time of simulation
    if(is.null(trunc.time)==TRUE){
        time.idxs <- 1:length(simn$anc.times) # If no specification include all
    }else{
        if(simn$T.fin < trunc.time){time.idxs <- numeric(0)}else{time.idxs <- length(simn$anc.times[simn$anc.times <= trunc.time]):length(simn$anc.times)}} # If specified include specified number
    if(is.null(trunc.time.range)==FALSE){
        time.idxs <- length(simn$anc.times[simn$anc.times <= trunc.time.range[1]]):length(simn$anc.times[simn$anc.times <= trunc.time.range[2]])} # If a range of idxs to be included has been specified then choose specified range
    #### 2 #### Find intersection of ranges and define storage matrices
    idxs <- intersect(frac.idxs,time.idxs); p.num <- simn$p.num; l.idxs <- length(idxs); c.idxs <- 1:(p.num*l.idxs); r.idxs <- numeric(0); for(i in 1:l.idxs){r.idxs <- c(r.idxs,(0:(p.num-1))*l.idxs+i)} # Define various index extractions
    idxs <- intersect(idxs,c(1:length(simn$anc.times))*retain.every) # Redefine indexes to be extracted based on frequency determined by user
    times <- simn$anc.times[idxs] # Define times according to indicies
    simn.draws <- simn$p.anc.mat[,,idxs,drop=FALSE] # Extract particle locations corresponding to indicies
    simn.weights <- simn$p.anc.wei[idxs,,drop=FALSE] # Extract particle weights according to indicies
    if(untransform==TRUE){for(j in 1:length(idxs)){simn.draws[,,j] <- apply(simn.draws[,,j],2,un_scale_transform)}} # If required, transform particles from transformed to original space
    if(even.weights==TRUE){ # If particles have to be of even weight
        for(j in 1:length(idxs)){ # If particles have to be of even weight, resample particles
            simn.draws[,,j] <- simn.draws[,sample(simn$p.num,simn$p.num,replace=TRUE,prob=simn.weights[j,]),j]
            simn.weights[j,] <- 1/simn$p.num}
    }else{for(j in 1:length(idxs)){simn.weights[j,] <- 1/sum(simn.weights[j,])}} # If particles are not of even weight, ensure each time slice has even weighting
    list(idxs=idxs,times=times,weights=simn.weights,draws=simn.draws)}

########################################################################
########################################################################
#### 7 - Ergodic Output Effective Sample Size Calculation
########################################################################
########################################################################

scale_ess <- function(ergodic_output){
    number.times <- length(ergodic_output$times) # Determine the number of time slices included within the ergodic output
    number.dimen <- dim(ergodic_output$draws)[1] # Determine the number of dimensions
    p.means <- p.sds <- matrix(0, number.dimen,number.times) # Storage matrices for particle means and sds
    for(t in 1:number.times){for(d in 1:number.dimen){p.means[d,t] <- mean(ergodic_output$draws[d,,t]); p.sds[d,t] <- sd(ergodic_output$draws[d,,t])}} # Determine mean and sds for each time slice and dimension
    ess <- numeric(number.dimen) # Storage vector for ESS proxy
    for(d in 1:number.dimen){ # Computation of ESS by dimension
        dimen.pos.mean <- mean(p.means[d,]) # Mean of time slice means
        dimen.mean.var <- (p.means[d,]-dimen.pos.mean)^2 # Variance of time slice means
        dimen.pos.var <- mean((p.sds[d,]^2)+dimen.mean.var) # Variance of time slice
        dimen.mar.ess <- dimen.pos.var/mean(dimen.mean.var) # Raw ESS
        dimen.act <- (1+acf(p.means[d,])[[1]][2])/(1-acf(p.means[d,])[[1]][2]) # Dimension autocorrelation time
        ess[d] <- dimen.mar.ess/dimen.act} # ESS accounting for autocorrelation
    list(ess=ess,p.means=p.means,p.sds=p.sds)}
