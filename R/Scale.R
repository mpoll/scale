############################################
############################################
######       Scale Algorithm         #######
######  Last Updated: 30/05/2018 MP  #######
############################################
############################################

############################################
############################################
### -1- C++ preliminaries
############################################
############################################

library(Rcpp); #sourceCpp('../src/scale_rcpp.cpp') # Call CPP code

#############################################
#############################################
### 0 - Libraries, Global Functions & Euler
#############################################
#############################################

library(gtools); library(Rmpfr); library(statmod); library(abind); library(MASS)

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

#############################################
#############################################
### 1 - Resampling Steps
#############################################
#############################################

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

multi.resamp 	<- function(p.wei,n=length(p.wei)){ # Multinomial Resampling
    if((sum(p.wei)>0)&(n>0)){ # Check whether resampling possible
        p.idx		<- sample(1:length(p.wei),n,replace=TRUE,prob=p.wei) # Sampled Index
    }else{p.idx 	<- numeric(0)} # If resampling not possible return empty vector
    list(p.idx=p.idx)} # Return particle indices

#############################################
##### 1.2.2 - Systematic Resampling
#############################################

system.resamp	<- function(p.wei,n=length(p.wei)){ # Systematic Resampling
    if((sum(p.wei)>0)&(n>0)){ # Check whether resampling possible
        cum.wei		<- cumsum(p.wei)/sum(p.wei) # Normalise and calculate cumulative weights
        samp.vec	<- seq(runif(1,0,1/n),1,1/n) # Vector of cdf samples
        p.idx		<- numeric(0); for(i in 1:n){p.idx <- c(p.idx,length(cum.wei[samp.vec[i]>=cum.wei])+1)} # Sample from cdf
    }else{p.idx 	<- numeric(0)} # If resampling not possible return empty vector
    list(p.idx=p.idx)} # Return particle indices

#############################################
##### 1.2.3 - Stratified Resampling
#############################################

strat.resamp	<- function(p.wei,n=length(p.wei)){ # Stratified Resampling
    vec.length  <- length(p.wei) # Calculate the length of the input vector
    cum.wei		<- cumsum(p.wei) # Calculate cumulative weights
    if(cum.wei[vec.length]>0){ # Check whether resampling possible
        cum.wei     <- cum.wei/cum.wei[vec.length]  # Normalise cumulative weights
        samp.vec	<- seq(0,(n-1)/n,1/n) + runif(n,0,1/n) # Vector of cdf samples
        p.idx <- findInterval(samp.vec,cum.wei)+1 # Sample from cdf
    }else{p.idx 	<- numeric(0)} # If resampling not possible return empty vector
    list(p.idx=p.idx)} # Return particle indices

#############################################
##### 1.2.4 - Residual Resampling
#############################################

resid.resamp	<- function(p.wei,n=length(p.wei),nest.resamp=strat.resamp){ # Residual Resampling
    if((sum(p.wei)>0)&(n>0)){ # Check whether resampling possible
        det.resamp	<- floor(n*p.wei/sum(p.wei))
        p.idx		<- c(rep(1:length(p.wei),det.resamp),nest.resamp(p.wei-det.resamp*(sum(p.wei)/n),n=n-sum(det.resamp))$p.idx)
    }else{p.idx 	<- numeric(0)} # If resampling not possible return empty vector
    list(p.idx=p.idx)} # Return particle indices

#############################################
#############################################
#### 2 - Exact Algorithm Constructions
#############################################
#############################################

#############################################
#### 2.1 - Elementary EA Functions
#############################################

#############################################
##### 2.1.1 - Bessel Functionals - See Pollock et al. 2015
#############################################

eazeta		<- function(n,s,t,x,y,L,U){if(max(x-U,y-U,L-x,L-y)>=0){1}else{j<-1:(ceiling(n/2));P<--2/(t-s);D<-U-L;D1<-D*j+L;D2<-D*j-U;z<-y-x;if(even(n)){sum(exp(P*(D1-x)*(D1-y))+exp(P*(D2+x)*(D2+y))-exp(P*j^2*D^2-P*j*D*z)-exp(P*j^2*D^2+P*j*D*z))}else{sum(exp(P*(D1-x)*(D1-y))+exp(P*(D2+x)*(D2+y)))-sum(exp(P*j[1:length(j)-1]^2*D^2-P*j[1:length(j)-1]*D*z)+exp(P*j[1:length(j)-1]^2*D^2+P*j[1:length(j)-1]*D*z))}}} # Zeta Functional
eagamma		<- function(n,s,t,x,y,L,U){1-eazeta(n,s,t,x,y,L,U)} # Gamma Functional
eapsi		<- function(j,s,t,m,xoy,u){P<--2*abs(u-m)*j/(t-s);(2*abs(u-m)*j-(xoy-m))*exp(P*(abs(u-m)*j-(xoy-m)))} # Psi Functional
eachi		<- function(j,s,t,m,xoy,u){P<--2*abs(u-m)*j/(t-s);(2*abs(u-m)*j+(xoy-m))*exp(P*(abs(u-m)*j+(xoy-m)))} # Chi Functional
eadel2  	<- function(n,s,t,m,xoy,u){if(max(xoy-u,m-xoy)>=0){0}else{if(even(n)){j<-1:(n/2);1-(sum(eapsi(j,s,t,m,xoy,u)-eachi(j,s,t,m,xoy,u)))/(xoy-m)}else{if(n>1){j<-1:((n-1)/2);1-(sum(eapsi(j,s,t,m,xoy,u)-eachi(j,s,t,m,xoy,u)))/(xoy-m)-eapsi(max(j)+1,s,t,m,xoy,u)/(xoy-m)}else{1-eapsi(1,s,t,m,xoy,u)/(xoy-m)}}}}
eadel		<- function(n,s,t,x,y,m,u){if(x==m){xI<-1}else{xI<-0};if(y==m){yI<-1}else{yI<-0};if(max(xI,yI)==1){delT<-2}else{delT<-1};if(m>max(x,y)){x<--x;y<--y;m<--m;u<--u};if(max(x-u,y-u,m-x,m-y)>=0){out<-0};if(delT==1){out<-eagamma_cpp(n,s,t,x,y,m,u)/(1-exp(-2*(x-m)*(y-m)/(t-s)))};if(delT==2){if(xI*yI==0){xoy<-max(x,y);out<-eadel2_cpp(n,s,t,m,xoy,u)}else{out<-0}};if(out<0){out<-0};if(out>1){out<-1};if((t-s)==0){out <- 1}; out}
eadelR		<- function(n,s,t,x,y,m,u){if(x==m){xI<-1}else{xI<-0};if(y==m){yI<-1}else{yI<-0};if(max(xI,yI)==1){delT<-2}else{delT<-1};if(m>max(x,y)){x<--x;y<--y;m<--m;u<--u};if(max(x-u,y-u,m-x,m-y)>=0){out<-0};if(delT==1){out<-eagamma(n,s,t,x,y,m,u)/(1-exp(-2*(x-m)*(y-m)/(t-s)))};if(delT==2){if(xI*yI==0){xoy<-max(x,y);out<-eadel2(n,s,t,m,xoy,u)}else{out<-0}};if(out<0){out<-0};if(out>1){out<-1};if((t-s)==0){out <- 1}; out}
eadelC 		<- function(mt,s,t,x,y,m,u){if(mt>=mpfrthr){pbn<-mt*mpfrpbn;s<-mpfr(s,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);y<-mpfr(y,precBits=pbn);m<-mpfr(m,precBits=pbn);u<-mpfr(u,precBits=pbn);c(s1=eadelR(mt,s,t,x,y,m,u),s2=eadelR(mt+1,s,t,x,y,m,u))}else{eadel_pair_cpp(mt,s,t,x,y,m,u);}}

#############################################
#### 2.1.2 - Simulate Bessel 3 Process Mid Points - See Pollock (PhD Thesis, 2013)
#############################################

eabesmid	<- function(q,s,tau,t,x,m,y,minI){
	if(minI==-1){x<--x;y<--y;m<--m} # Reflection
    bI<-0;if(q==s){bI<-1;w<-x};if(q==tau){bI<-1;w<-m};if(q==t){bI<-1;w<-y} # Boundary
    t<-t-s;tau<-tau-s;q<-q-s # Rescale time
    if(bI==0){if(q<tau){Ra1<-sqrt(tau);Ra2<-(x-m)*(tau-q)/((tau)^(3/2));Ra3<-(tau-q)/(tau)}else{Ra1<-sqrt(t-tau);Ra2<-(y-m)*(q-tau)/((t-tau)^(3/2));Ra3<-(q-tau)/(t-tau)};BB3<-rnorm(3,0,sqrt(Ra3*(1-Ra3)));w<-m+Ra1*sqrt((Ra2+BB3[1])^2+(BB3[2])^2+(BB3[3])^2)}
	list(w=minI*w)}

#############################################
##### 2.1.3 - Bessel Exceedance Evaluation - See Pollock (PhD Thesis, 2013)
#############################################

eabesex		<- function(sV,tV,xV,yV,m,B1,B2,minI){ # Vectors of equal length
    if(minI==-1){xV<--xV;yV<--yV;m<--m;B1<--B1;B2<--B2}
    u <- runif(1,0,1); mt<-ceiling(max(sqrt(max(tV-sV)+(B1-m)^2)/(2*(B1-m)),sqrt(max(tV-sV)+(B2-m)^2)/(2*(B2-m))))
    em<-matrix(0,length(sV),8); em[,1]<-sV; em[,2]<-tV; em[,3]<-xV; em[,4]<-yV
    B1evI<-B2evI<-0; while(B1evI==0){
        emM <- em; for(i in 1:dim(em)[1]){if(mt>=mpfrthr) {em[i,5:6]<-as.numeric(eadelC(mt,em[i,1],em[i,2],em[i,3],em[i,4],m,B1))}else{em[i,5:6] <- as.numeric(eadel_pair_cpp(mt,em[i,1],em[i,2],em[i,3],em[i,4],m,B1))}}
        if(u<=prod(em[,5])){B1evI<-B2evI<-1;con1I<-1;con2I<-1;ex1I<-0;ex2I<-0}else{if(u>prod(em[,6])){B1evI<-1;con1I<-0;ex1I<-1}else{B1evI<-0;con1I<-0;ex1I<-0;mt<-mt+2}}}
    while(B2evI==0){for(i in 1:dim(em)[1]){
        if(mt >= mpfrthr){   em[i,7:8]<-as.numeric(eadelC(mt,em[i,1],em[i,2],em[i,3],em[i,4],m,B2))} else { em[i,7:8]<-as.numeric(eadel_pair_cpp(mt,em[i,1],em[i,2],em[i,3],em[i,4],m,B2))}
    };if(u<=prod(em[,7])){B2evI<-1;con2I<-1;ex1I<-0}else{if(u>prod(em[,8])){B2evI<-1;con2I<-0;ex2I<-1}else{B2evI<-0;con2I<-0;ex2I<-0;mt<-mt+2}}}
    if(minI==-1){em[,3]<--em[,3];em[,4]<--em[,4]}; accI<-0; if(con1I==1){accI<-1}else{if(con2I==1){if(rbinom(1,1,0.5)==1){accI<-1}}}
    list(accI=accI,u=u,con1I=con1I,con2I=con2I,em=em)}

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

dev.pr		<- function(Jst.t=0.64,Jst.rat=0.5776972){ # First Hitting Time (J*) Proposal Function
	U <- runif(1,0,1) # Simulate uniform
    if(U < Jst.rat){ # Case 1, U <- p/(p+q)
         E <- rexp(1,1); X <- Jst.t + 8*E/(pi^2)
    }else{ # Case 2, U >= p/(p+q)
    repeat{E <- rexp(2,1); if(E[1]^2 <= 2*E[2]/Jst.t){X <- Jst.t/(1+Jst.t*E[1])^2; break}}}
    list(X=X,U=U,E=E)}

dev.rej		<- function(X){ # First Hitting Time (J*) Rejection Sampler
	S <- pi*exp(-pi^2*X/8)/2; n <- 0; Y <- runif(1,0,1)*S # Initialise sequence and simulate Uniform
    repeat{ # Loop
        n <- n + 1; S <- S - (-1)^n*pi*(n+1/2)*exp(-(n+1/2)^2*pi^2*X/2); if(Y < S){Acc <- 1; break} # Odd n
        n <- n + 1; S <- S + (-1)^n*pi*(n+1/2)*exp(-(n+1/2)^2*pi^2*X/2); if(Y > S){Acc <- 0; break}} # Even n
     list(Acc=Acc,S=S,n=n,Y=Y,X=X)}

#############################################
##### 2.2.3 - First Hitting Time (J*) Simulation
#############################################

bm.pass		<- function(s=0,x=0,theta=1,Jst.t=0.64,Jst.rat=0.5776972){ # theta denotes passage level
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

p.path.renew <- function(pass.mat,cyc.mat,curr.time,x.curr,next.time,theta,dimen){ #***next.time <= tau.bar***
    if(dimen>1){
        x.next <- numeric(dimen) # Storage vector for x.next
        for(i in 1:dimen){repeat{
            x.next[i] <- eabesmid(next.time,curr.time,pass.mat[i,"tau"],pass.mat[i,"tau"],x.curr[i],pass.mat[i,"y"],pass.mat[i,"y"],pass.mat[i,"minI"])$w # Update sample path in given dimension to next.time
            if(eabesex(sV=c(curr.time,next.time),tV=c(next.time,pass.mat[i,"tau"]),xV=c(x.curr[i],x.next[i]),yV=c(x.next[i],pass.mat[i,"y"]),m=pass.mat[i,"y"],B1=x.curr[i]+pass.mat[i,"minI"]*theta[i],B2=x.curr[i]+pass.mat[i,"minI"]*theta[i],minI=pass.mat[i,"minI"])$accI==1){break}}} # If sample path retains alternate barrier then accept else reject
        p.taus <- p.path.init(s=next.time,x=x.next,theta=theta,dimen=dimen); pass.mat <- p.taus$pass.mat; cyc.mat <- p.taus$cyc.mat # Initialise particle first passage times
        list(curr.time=next.time,next.tau=cyc.mat[1,"tau"],max.next=cyc.mat[1,"tau"]-next.time,pass.mat=pass.mat,cyc.mat=cyc.mat,x.next=x.next)}
    else{
        repeat{x.next <- eabesmid(next.time,curr.time,pass.mat["tau"],pass.mat["tau"],x.curr,pass.mat["y"],pass.mat["y"],pass.mat["minI"])$w; if(eabesex(sV=c(curr.time,next.time),tV=c(next.time,pass.mat["tau"]),xV=c(x.curr,x.next),yV=c(x.next,pass.mat["y"]),m=pass.mat["y"],B1=x.curr+pass.mat["minI"]*theta,B2=x.curr+pass.mat["minI"]*theta,minI=pass.mat["minI"])$accI==1){break}} # If sample path retains alternate barrier then accept else reject
        new.fpt <- bm.pass(s=next.time,x=x.next,theta=theta,Jst.t=0.64,Jst.rat=0.5776972); pass.mat[] <- c(1,new.fpt$tau,new.fpt$y,new.fpt$minI,x.next-theta,x.next+theta)
        list(curr.time=next.time,next.tau=pass.mat["tau"],max.next=pass.mat["tau"]-next.time,pass.mat=pass.mat,cyc.mat=pass.mat,x.next=x.next)}}

#############################################
#############################################
### 3 - Scale Subfunctions - Event Rotational Scale Implementation
#############################################
#############################################

#############################################
#### 3.1 - Actual Phi Specification
#############################################

ac.phi          <- function(eta.pos,dapts){
    beta.pos <- un.scale.transform(eta.pos); data.evaluation <- data.eval(beta.pos,1:dsz) # Specify beta position and evaluate data
    (data.evaluation$grad.log.pi%*%t(data.evaluation$grad.log.pi) + sum(data.evaluation$lap.log.pi))/2}

#############################################
#### 3.2 - Actual Subsample Phi Specification
#############################################

##### 3.2.1 - Subsample Phi evaluation

ss.phi          <- function(eta.pos,dapts,ss.size){
    beta.pos <- un.scale.transform(eta.pos); factor <- dsz/ss.size # Specify beta position and factor
    data1 <- data.eval(beta.pos,dapts[1,,"dapt.1"],factor=factor); data1.star <- data.eval(beta.star,dapts[1,,"dapt.1"],factor=factor) # Evaluate data grad log and lap log (ss 1)
    data2 <- data.eval(beta.pos,dapts[1,,"dapt.2"],factor=factor); data2.star <- data.eval(beta.star,dapts[1,,"dapt.2"],factor=factor) # Evaluate data grad log and lap log (ss 2)
    data3 <- data.eval(beta.pos,dapts[1,,"dapt.3"],factor=factor); data3.star <- data.eval(beta.star,dapts[1,,"dapt.3"],factor=factor) # Evaluate data grad log and lap log (ss 3)
    ((data1$grad.log.pi - data1.star$grad.log.pi)%*%t(2*alpha.cent + data2$grad.log.pi - data2.star$grad.log.pi) + sum(data3$lap.log.pi)-sum(data3.star$lap.log.pi) + alpha.cent.sq + alpha.p.cent)/2}

##### 3.2.2 - Subsampling Phi intensity bounds

ss.phiC         <- function(p.mat.curr,l.bound=NULL,u.bound=NULL){ # Function of current particle location, lower bounds (dimensionally), upper bounds (dimensionally)
    l.bds <- pmin(0,l.bound); u.bds <- pmax(0,u.bound) # Extend the hyper cube to the origin (i.e. eta.star)
    distance <- pmax(abs(l.bds),abs(u.bds))            # Compute the maximal distance in each dimension
    alpha.bds <- alpha.prime.bds <- numeric(dimen); for(j in 1:dimen){ # Compute the bounds for alpha and alpha prime
        alpha.bds[j]        <- dsz*sum(distance*dim.grad.max[j,,1]) # alpha bounds
        alpha.prime.bds[j]  <- dsz*sum(distance*dim.grad.max[j,,2]) # alpha prime bounds
    }
    A.bds   <- alpha.bds^2 # Compute bounds on alpha^2
    phi.bds <- (sum(2*abs(alpha.cent)*alpha.bds)+sum(A.bds)+sum(alpha.prime.bds))/2 # Compute bounds (modulo normalisation) on subsampled phi
    phiL <- phi.cent - phi.bds; phiU <- phi.cent + phi.bds; intensity <- phiU - phiL # Compute Bounding functionals and intensity
    list(distance=distance,alpha.bds=alpha.bds,alpha.prime.bds=alpha.prime.bds,A.bds=A.bds,phiL=phiL,phiU=phiU,intensity=intensity)}


#############################################
#### 3.3 - TURN OFF THE SUB SAMPLER (REQUIRES SETTING ss.on <- FALSE)
#############################################

if(exists("ss.on")==TRUE){if(ss.on==FALSE){ss.phi <- function(eta.pos,dapts,ss.size){ac.phi(eta.pos,dapts)}}}

#############################################
#### 3.4 - Approximate Algorithm Apply Functionals
#############################################

approx.intensity  <- function(eta.pos){phiCfn <- ss.phiC(eta.pos); c(phiCfn$phiL,phiCfn$phiU)}

#############################################
#### 3.5 - Resampling Step - Event Rotational Implementation
#############################################

scale_resamp <- function(p.num,p.idx,log.p.wei,p.mat,p.layer,p.layer.sor.I,p.dapts,ss.size,resamp.method,ess.thresh,p.pass.arr=NULL,p.cyc.arr=NULL){
    ###### 1 ###### Perform weight control function - Subtract maximum log weight to control floating point precision, normalise and calculate ESS
    wei.control <- part.ess(log.p.wei,p.num)
    ###### 2 ###### Subtract maximum log weight to control floating point precision
    if(wei.control$ess <= ess.thresh){ # If ESS < Threshold
        ### Determine resampling indices
        r.idx <- resamp.method(wei.control$p.wei)$p.idx # Update Resampled Particle Indices
        ### Determine ancestoral lineages
        p.idx <- p.idx[r.idx] # Update particle ancestory
        ### Update p.mat and p.dapts (ordered by c.idx)
        p.mat <- p.mat[,r.idx,drop=FALSE] # Update particle matrix
        p.dapts[] <- sample(dsz,3*p.num*ss.size,replace=TRUE) # Update particle dapt matrix *** IE REFRESH DATA POINTS ***
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
    list(p.idx=p.idx,r.idx=r.idx,log.p.wei=log.p.wei,p.mat=p.mat,p.layer=p.layer,p.dapts=p.dapts,resamp.I=resamp.I,ess=wei.control$ess,p.pass.arr=p.pass.arr,p.cyc.arr=p.cyc.arr)}

#############################################
#### 3.6 - Negative Weight Handling Sub-Algorithms (Approx Alg.)
#############################################

#############################################
######## 3.3.1 - Set weights to zero
#############################################

scale_zero.wei  <- function(p.event.ch,p.curr.idx,p.curr,p.anc.neg,log.p.wei,ss.phiC,ss.phiC.up){
    p.anc.neg.append <- c(p.event.ch,exp(log.p.wei[p.curr.idx]),p.curr[c("t","PhiL","PhiU")]) # Ancestral Negative Weight Record
    list(p.anc.neg.append=p.anc.neg.append,ss.phiC=ss.phiC)}

#############################################
######## 3.3.3 - Set weights and index PhiC
#############################################

scale_zero.inc  <- function(p.event.ch,p.curr.idx,p.curr,p.anc.neg,log.p.wei,ss.phiC,ss.phiC.up){
    p.anc.neg.append <- c(p.event.ch,exp(log.p.wei[p.curr.idx]),p.curr[c("t","PhiL","PhiU")]) # Ancestral Negative Weight Record
    p.phi <- p.curr["PhiU"]-p.event.ch*p.curr["Delta"] # Compute what phi was
    print(sprintf("Old PhiL:%f, Old PhiU:%f",p.curr["PhiL"],p.curr["PhiU"]))
    ss.phiC <- ss.phiC.up(p.phi,p.curr["PhiL"],p.curr["PhiU"]) # Index intensity
    print(sprintf("New PhiL:%f, New PhiU:%f",ss.phiC(1)$phiL,ss.phiC(1)$phiU))
    list(p.anc.neg.append=p.anc.neg.append,ss.phiC=ss.phiC)}

#############################################
######## 3.3.3.1 - ss.phiC.up default index for Lipschitz likelihood
#############################################
ss.phiC.up      <- function(p.phi,p.phiL,p.phiU){function(p.mat.curr){list(phiL=p.phiL-2*(p.phiU-p.phiL),phiU=p.phiU+2*(p.phiU-p.phiL))}}

#############################################
#############################################
### 4 - Scale Algorithm
#############################################
#############################################

#############################################
#### 4.1 - Scale Algorithm Approximate Version
#############################################

scale_approx <- function(p.num,t.inc,T.fin,ss.phi,ss.phiC,dimen,scale.transform,un.scale.transform,T.start=0,x.init=NULL,ss.size=1,ess.thresh=0.5,resamp.method=resid.resamp,neg.wei.mech=scale_zero.wei,prev.simn=NULL,progress.check=FALSE,phi.record=FALSE,resamp.freq=p.num-1,theta=NULL,p.path.renew=p.path.renew){
    #################
    #### (0)1 #### Initalise Algorithm
    #################
    ##### (0)1.1 #### Define dimension and other variables not pre-defined in input
    p.mu <- matrix(c(rep(0,dimen)),ncol=1,dimnames=list(sprintf("dim.%i",1:dimen),"mean")); p.covar <- diag(dimen); dimnames(p.covar) <- list(sprintf("dim.%i",1:dimen),sprintf("dim.%i",1:dimen)) # Specify algorithm dimension and Normal covariance matrix
    ##### (0)1.2 Initialise particle storage framework
    if(is.null(prev.simn)==TRUE){ #i.e. starting algorithm from 0
        ### Algorithm initialisation point
        if(is.null(x.init)==TRUE){p.mat <- matrix(0,dimen,p.num,dimnames=list(sprintf("dim.%i",1:dimen),NULL))}else{
            if(is.vector(x.init)){x.init <- matrix(x.init,length(x.init),1)} # If x.init is a vector transform to matrix
            if(dim(x.init)[2]==1){x.init <- matrix(rep(x.init,p.num),dimen,p.num)} # If x.init is a point then repeat
            p.mat <- apply(x.init,2,scale.transform); dimnames(p.mat)<-list(sprintf("dim.%i",1:dimen),NULL) # Transform x.init and define as p.mat
        } # Setting of p.mat depending on x.init NULL, d dimn vector, dxp dimn vector, 1xd matrix or pxd matrix
        ### Define storage matrices
        p.idx <- 1:p.num; log.p.wei <- rep(log(1/p.num),p.num) # Set idxs and weights
        p.layer <- matrix(0,7,p.num,dimnames=list(c("c.idx","s","deg.s","t","PhiL","PhiU","Delta"),NULL)) # Establish particle layer matrix
        p.layer["c.idx",] <- p.idx # Set particle idxs in layer matrix
        p.dapts <- array(sample(dsz,3*p.num*ss.size,replace=TRUE),c(p.num,ss.size,3),dimnames=list(sprintf("part.%i",1:p.num),sprintf("ss.%i",1:ss.size),sprintf("dapt.%i",1:3))) # Set p.layer data points
    }else{ #i.e. restarting algorithm from some point
        T.start <- prev.simn$T.fin; p.idx <- prev.simn$p.idx; log.p.wei <- prev.simn$log.p.wei; p.mat <- prev.simn$p.mat; p.layer <- prev.simn$p.layer; p.dapts <- prev.simn$p.dapts} # Restablish Algorithm
    ##### (0)1.3 Initialise particle system
    p.layer[c("PhiL","PhiU"),] <- apply(p.mat,2,approx.intensity); p.layer[c("Delta"),] <- p.layer[c("PhiU"),] - p.layer[c("PhiL"),] # Set p.layer current particle intensity
    p.bet.degrad.r <- p.bet.degrad <- numeric(p.num) # Initialise between resampling log weight degration
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
            p.resamp <- scale_resamp(p.num,p.idx,log.p.wei,p.mat,p.layer,p.layer.sor.I=1,p.dapts,ss.size,resamp.method,ess.thresh) # Resample particles
            p.anc.ess <- c(p.anc.ess,p.resamp$ess) # Record ancestoral ESS at mesh point
            p.idx <- p.resamp$p.idx; log.p.wei <- p.resamp$log.p.wei # Update particle indicies and weights (normalised)
            if(p.resamp$resamp.I==1){p.anc.resamp <- c(p.anc.resamp,curr.deg.time); p.mat <- p.resamp$p.mat; p.layer <- p.resamp$p.layer; p.dapts <- p.resamp$p.dapts} # If there has been resampling, update ancestral resampling times, layer and dapts
            p.layer <- mat.sort.r(p.layer,p.num,"t") # Resort p.layer by time
            resamp.counter <- 1 # Reset resample counter
            if(progress.check==TRUE){print(curr.deg.time)} # Output current time if checking progress
        }else{ # If we haven't reached the next increment point
            ###### (t.inc)2.2.b.1 ##### Select particle to be updated
            p.curr <- p.layer[,1] # Next particle to be updated
            p.curr.idx <- p.curr["c.idx"] # Define current index
            p.curr.dapts <- p.dapts[p.curr.idx,,,drop=FALSE] # Find particle subsampled data points
            ###### (t.inc)2.2.b.2 ##### Update weight
            p.loc.next <- p.mat[,p.curr.idx] + sqrt(p.curr["t"]-p.curr["s"])*mvrnorm(1,p.mu,p.covar) # Particle location at t
            p.phi.eval <- ss.phi(p.loc.next,p.curr.dapts,ss.size)
            p.event.ch <- (p.curr["PhiU"] - p.phi.eval)/p.curr["Delta"]; if(phi.record==TRUE){p.phi.hist <- rbind(p.phi.hist,c(p.phi.eval,p.loc.next))} # Raw Poisson event increment
            if(p.event.ch < 0){p.neg <- neg.wei.mech(p.event.ch,p.curr.idx,p.curr,p.anc.neg,log.p.wei,ss.phiC,ss.phiC.up); log.p.wei[p.curr.idx] <- -Inf; p.anc.neg <- rbind(p.anc.neg,p.neg$p.anc.neg.append); ss.phiC <- p.neg$ss.phiC}else{log.p.wei[p.curr.idx] <- log.p.wei[p.curr.idx] + log(p.event.ch)} # Apply chosen negative weight mechanism
            p.bet.degrad <- cbind(p.bet.degrad,c(p.curr.idx,(p.curr["t"]-p.curr["deg.s"])*p.curr[c("PhiL")])) # Record individual particle degradation (since last resampling point)
            ###### (t.inc)2.2.b.3 ##### Update particle to new point
            p.mat[,p.curr.idx] <- p.loc.next # Update particle location
            p.dapts[p.curr.idx,,] <- sample(dsz,3*ss.size,replace=TRUE) # Update particle data points
            p.phi.bds <- ss.phiC(p.loc.next); p.curr[c("PhiL","PhiU")] <- c(p.phi.bds$phiL,p.phi.bds$phiU); p.curr["Delta"] <- p.curr["PhiU"] - p.curr["PhiL"] # Update phi bounds and intensity
            p.curr["s"] <- p.curr["deg.s"] <- p.curr["t"] # Update current particle degradation time
            p.curr["t"] <- p.curr["deg.s"] + rexp(1,rate=p.curr["Delta"]) # Update current particle event time
            p.layer <- ins.sort.r(p.curr,p.layer[,-1],p.num-1,"t") # Re-insert particle into p.layer
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
                p.resamp <- scale_resamp(p.num,p.idx,log.p.wei,p.mat,p.layer,p.layer.sor.I=0,p.dapts,ss.size,resamp.method,ess.thresh) # Resample particles
                p.idx <- p.resamp$p.idx; log.p.wei <- p.resamp$log.p.wei # Update particle indicies and weights (normalised)
                if(p.resamp$resamp.I==1){p.anc.resamp <- c(p.anc.resamp,curr.deg.time); p.mat <- p.resamp$p.mat; p.layer <- p.resamp$p.layer; p.dapts <- p.resamp$p.dapts} # If there has been resampling, update ancestral resampling times, layer and dapts
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
    list(p.num=p.num,p.idx=p.idx,log.p.wei=log.p.wei,p.mat=p.mat,p.layer=p.layer,theta=theta,p.dapts=p.dapts,p.anc.idx=p.anc.idx,p.anc.wei=p.anc.wei,p.anc.mat=p.anc.mat,resamp.freq=resamp.freq,resamp.method=resamp.method,neg.wei.mech=neg.wei.mech,p.anc.resamp=p.anc.resamp,p.anc.ess=p.anc.ess,p.anc.neg=p.anc.neg,p.phi.hist=p.phi.hist,anc.times=anc.times,T.start=T.start,t.inc=t.inc,T.fin=T.fin,dimen=dimen,p.mu=p.mu,ss.size=ss.size,ess.thresh=ess.thresh,ss.phi=ss.phi,ss.phiC=ss.phiC,p.cyc.arr=NULL,p.pass.arr=NULL,scale.transform=scale.transform,un.scale.transform=un.scale.transform,progress.check=progress.check,phi.record=phi.record,p.path.renew=p.path.renew)} # Output list

#############################################
#### 4.2 - Scale Algorithm Exact Version
#############################################

### theta <- NULL; prev.simn <- NULL; x.init <- NULL; T.start <- 0; theta <- NULL; phi.record <- FALSE; ess.thresh <- 0.5

scale_exact <- function(p.num,t.inc,T.fin,ss.phi,ss.phiC,dimen,scale.transform,un.scale.transform,T.start=0,x.init=NULL,ss.size=1,ess.thresh=0.5,resamp.method=resid.resamp,neg.wei.mech=scale_zero.wei,prev.simn=NULL,progress.check=FALSE,phi.record=FALSE,resamp.freq=p.num-1,theta=NULL,p.path.renew=p.path.renew){
    #################
    #### (0)1 #### Initalise Algorithm
    #################
    ##### (0)1.1 #### Define dimension and other variables not pre-defined in input
    p.mu <- matrix(c(rep(0,dimen)),ncol=1,dimnames=list(sprintf("dim.%i",1:dimen),"mean")); p.covar <- diag(dimen); dimnames(p.covar) <- list(sprintf("dim.%i",1:dimen),sprintf("dim.%i",1:dimen)) # Specify algorithm dimension and Normal covariance matrix
    if(is.null(theta)){theta <- rep(signif(dimen^(1/4)*sqrt(t.inc),digits=1),dimen)} # Specify theta
    ##### (0)1.2 Initialise particle storage framework
    if(is.null(prev.simn)==TRUE){ #i.e. starting algorithm from 0
        ### Algorithm initialisation point
        if(is.null(x.init)==TRUE){p.mat <- matrix(0,dimen,p.num,dimnames=list(sprintf("dim.%i",1:dimen),NULL))}else{
            if(is.vector(x.init)){x.init <- matrix(x.init,length(x.init),1)} # If x.init is a vector transform to matrix
            if(dim(x.init)[2]==1){x.init <- matrix(rep(x.init,p.num),dimen,p.num)} # If x.init is a point then repeat
            p.mat <- apply(x.init,2,scale.transform); dimnames(p.mat)<-list(sprintf("dim.%i",1:dimen),NULL) # Transform x.init and define as p.mat
        } # Setting of p.mat depending on x.init NULL, d dimn vector, dxp dimn vector, 1xd matrix or pxd matrix
        ### Define storage matrices
        p.idx <- 1:p.num; log.p.wei <- rep(log(1/p.num),p.num) # Set idxs and weights
        p.layer <- matrix(0,9,p.num,dimnames=list(c("c.idx","s","deg.s","t","tau.bar","next.ev","PhiL","PhiU","Delta"),NULL)) # Establish particle layer matrix
        p.layer["c.idx",] <- p.idx # Set particle idxs in layer matrix
        p.dapts <- array(sample(dsz,3*p.num*ss.size,replace=TRUE),c(p.num,ss.size,3),dimnames=list(sprintf("part.%i",1:p.num),sprintf("ss.%i",1:ss.size),sprintf("dapt.%i",1:3))) # Set p.layer data points
        p.pass.arr <- p.cyc.arr <- array(0,c(dimen,6,p.num),dimnames=list(NULL,c("dimen","tau","y","minI","l","u"),sprintf("part.%i",1:p.num))) # Define particle passage and cyclic arrays
        for(i in 1:p.num){p.taus <- p.path.init(s=T.start,x=p.mat[,i],theta=theta,dimen=dimen); p.pass.arr[,,i] <- p.taus$pass.mat; p.cyc.arr[,,i] <- p.taus$cyc.mat} # Initialise particle first passage times
        p.layer["tau.bar",] <- p.cyc.arr[1,"tau",] # Determine each particles tau.bar
    }else{ #i.e. restarting algorithm from some point
        T.start <- prev.simn$T.fin; p.idx <- prev.simn$p.idx; log.p.wei <- prev.simn$log.p.wei; p.mat <- prev.simn$p.mat; p.layer <- prev.simn$p.layer; p.dapts <- prev.simn$p.dapts; p.pass.arr <- prev.simn$p.pass.arr; p.cyc.arr <- prev.simn$p.cyc.arr} # Restablish Algorithm
    ##### (0)1.3 Initialise particle system
    for(i in 1:p.num){p.phi.bds <- ss.phiC(p.mat[,i],p.pass.arr[,"l",i],p.pass.arr[,"u",i]); p.layer[c("PhiL","PhiU"),i]<-c(p.phi.bds$phiL,p.phi.bds$phiU)}; p.layer[c("Delta"),] <- p.layer["PhiU",] - p.layer["PhiL",] # Set p.layer current particle intensity
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
    anc.times <- curr.deg.time <- T.start; t.inc.next <- T.start + t.inc; resamp.counter <- 1 # Initialise current degradation time, time of first increment and resampling counter
    while(curr.deg.time < T.fin){ # Iterate until final time is reached
        ##### (t.inc)2.1 ##### Increment particle set to next event time
        ######################
        if(p.layer["next.ev",1] >= t.inc.next){ # If we have reached the next increment point
            ###### (t.inc)2.2.a.1 ##### Compute particles at intermediate time point
            p.layer <- mat.sort.r(p.layer,p.num,"c.idx") # Sort p.layer matrix by current particle index
            for(loc.up in 1:p.num){p.mat[,loc.up]<-p.path.renew(pass.mat=p.pass.arr[,,loc.up],cyc.mat=p.cyc.arr[,,loc.up],curr.time=p.layer["s",loc.up],x.curr=p.mat[,loc.up],next.time=t.inc.next,theta,dimen)$x.next}
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
            p.resamp <- scale_resamp(p.num,p.idx,log.p.wei,p.mat,p.layer,p.layer.sor.I=1,p.dapts,ss.size,resamp.method,ess.thresh,p.pass.arr=p.pass.arr,p.cyc.arr=p.cyc.arr) # Resample particles
            p.anc.ess <- c(p.anc.ess,p.resamp$ess) # Record ancestoral ESS at mesh point
            p.idx <- p.resamp$p.idx; log.p.wei <- p.resamp$log.p.wei # Update particle indicies and weights (normalised)
            if(p.resamp$resamp.I==1){p.anc.resamp <- c(p.anc.resamp,curr.deg.time); p.mat <- p.resamp$p.mat; p.layer <- p.resamp$p.layer; p.dapts <- p.resamp$p.dapts; p.pass.arr <- p.resamp$p.pass.arr; p.cyc.arr <- p.resamp$p.cyc.arr} # If there has been resampling, update ancestral resampling times, layer, dapts and cycling arrays
            p.layer <- mat.sort.r(p.layer,p.num,"next.ev") # Redefine storage vectors II
            resamp.counter <- 1 # Reset resample counter
            if(progress.check==TRUE){print(curr.deg.time)} # Output current time if checking progress
        }else{ # If we haven't reached the next increment point
            ###### (t.inc)2.2.b.1 ##### Select particle to be updated
            p.curr <- p.layer[,1] # Next particle to be updated
            p.curr.idx <- p.curr["c.idx"] # Define current index
            p.curr.dapts <- p.dapts[p.curr.idx,,,drop=FALSE] # Find particle subsampled data points
            ###### (t.inc)2.2.b.2 ##### Update particle trajectory
            traj.update <- p.path.renew(pass.mat=p.pass.arr[,,p.curr.idx],cyc.mat=p.cyc.arr[,,p.curr.idx],curr.time=p.curr["s"],x.curr=p.mat[,p.curr.idx],next.time=p.curr["next.ev"],theta,dimen) # Iterate passage matrix
            p.mat[,p.curr.idx] <- p.loc.next <- traj.update$x.next # Particle location at next.ev
            ###### (t.inc)2.2.b.3 ##### If particle to be updated has an event before the next tau.bar, then update to Poisson event time
            if(p.curr["t"]<=p.curr["tau.bar"]){
                ###### (t.inc)2.2.b.3.1 ##### Update data set
                p.dapts[p.curr.idx,,] <- sample(dsz,3*ss.size,replace=TRUE) # Update particle data points
                ###### (t.inc)2.2.b.3.2 ##### Update Poisson weight
                p.phi.eval <- ss.phi(p.loc.next,p.curr.dapts,ss.size)
                p.event.ch <- (p.curr["PhiU"] - p.phi.eval)/p.curr["Delta"]; if(phi.record==TRUE){p.phi.hist <- rbind(p.phi.hist,c(p.phi.eval,p.loc.next))} # Raw Poisson event increment
                if(p.event.ch < 0){p.neg <- neg.wei.mech(p.event.ch,p.curr.idx,p.curr,p.anc.neg,log.p.wei,ss.phiC,ss.phiC.up); log.p.wei[p.curr.idx] <- -Inf; p.anc.neg <- rbind(p.anc.neg,p.neg$p.anc.neg.append); ss.phiC <- p.neg$ss.phiC}else{log.p.wei[p.curr.idx] <- log.p.wei[p.curr.idx] + log(p.event.ch)} # Apply chosen negative weight mechanism
            }else{
                ###### (t.inc)2.2.b.4 ##### If particle to be updated has no event before the next tau.bar, then update passage trajectory
                p.pass.arr[,,p.curr.idx] <- traj.update$pass.mat; p.cyc.arr[,,p.curr.idx] <- traj.update$cyc.mat # Update pass.mat and cyc.mat
                p.curr["tau.bar"] <- traj.update$next.tau # Determine next particle passage time
            }
            ###### (t.inc)2.2.b.5 ##### Update particle degradation matrix
            p.bet.degrad[p.curr.idx] <- p.bet.degrad[p.curr.idx] + (p.curr["next.ev"]-p.curr["deg.s"])*p.curr[c("PhiL")] # Record individual particle degradation (since last resampling point)
            ###### (t.inc)2.2.b.6 ##### Update particle Poisson intensity
            p.phi.bds <- ss.phiC(p.loc.next,p.pass.arr[,"l",i],p.pass.arr[,"u",i]); p.curr[c("PhiL","PhiU")] <- c(p.phi.bds$phiL,p.phi.bds$phiU); p.curr["Delta"] <- p.curr["PhiU"] - p.curr["PhiL"] # Update phi bounds and intensity
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
                p.resamp <- scale_resamp(p.num,p.idx,log.p.wei,p.mat,p.layer,p.layer.sor.I=0,p.dapts,ss.size,resamp.method,ess.thresh,p.pass.arr=p.pass.arr,p.cyc.arr=p.cyc.arr) # Resample particles
                p.idx <- p.resamp$p.idx; log.p.wei <- p.resamp$log.p.wei # Update particle indicies and weights (normalised)
                if(p.resamp$resamp.I==1){p.anc.resamp <- c(p.anc.resamp,curr.deg.time); p.mat <- p.resamp$p.mat; p.layer <- p.resamp$p.layer; p.dapts <- p.resamp$p.dapts; p.pass.arr <- p.resamp$p.pass.arr; p.cyc.arr <- p.resamp$p.cyc.arr} # If there has been resampling, update ancestral resampling times, layer and dapts
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
    list(p.num=p.num,p.idx=p.idx,log.p.wei=log.p.wei,p.mat=p.mat,p.layer=p.layer,theta=theta,p.dapts=p.dapts,p.anc.idx=p.anc.idx,p.anc.wei=p.anc.wei,p.anc.mat=p.anc.mat,resamp.freq=resamp.freq,resamp.method=resamp.method,neg.wei.mech=neg.wei.mech,p.anc.resamp=p.anc.resamp,p.anc.ess=p.anc.ess,p.anc.neg=p.anc.neg,p.phi.hist=p.phi.hist,anc.times=anc.times,T.start=T.start,t.inc=t.inc,T.fin=T.fin,dimen=dimen,p.mu=p.mu,ss.size=ss.size,ess.thresh=ess.thresh,ss.phi=ss.phi,ss.phiC=ss.phiC,p.cyc.arr=p.cyc.arr,p.pass.arr=p.pass.arr,scale.transform=scale.transform,un.scale.transform=un.scale.transform,progress.check=progress.check,phi.record=phi.record,p.path.renew=p.path.renew)} # Output list

### Modified version of Scale with Theta specification using global assignment

scale_exact <- function(p.num,t.inc,T.fin,ss.phi,ss.phiC,dimen,scale.transform,un.scale.transform,T.start=0,x.init=NULL,ss.size=1,ess.thresh=0.5,resamp.method=resid.resamp,neg.wei.mech=scale_zero.wei,prev.simn=NULL,progress.check=FALSE,phi.record=FALSE,resamp.freq=p.num-1,theta=NULL,p.path.renew=p.path.renew){
    #################
    #### (0)1 #### Initalise Algorithm
    #################
    ##### (0)1.1 #### Define dimension and other variables not pre-defined in input
    p.mu <- matrix(c(rep(0,dimen)),ncol=1,dimnames=list(sprintf("dim.%i",1:dimen),"mean")); p.covar <- diag(dimen); dimnames(p.covar) <- list(sprintf("dim.%i",1:dimen),sprintf("dim.%i",1:dimen)) # Specify algorithm dimension and Normal covariance matrix
    if(is.null(theta)){theta <- rep(signif(dimen^(1/4)*sqrt(t.inc)*(colMeans(dim.grad.max[,,1])^(-1)/max(colMeans(dim.grad.max[,,1])^(-1))),digits=1),dimen)} # Specify theta
    ##### (0)1.2 Initialise particle storage framework
    if(is.null(prev.simn)==TRUE){ #i.e. starting algorithm from 0
        ### Algorithm initialisation point
        if(is.null(x.init)==TRUE){p.mat <- matrix(0,dimen,p.num,dimnames=list(sprintf("dim.%i",1:dimen),NULL))}else{
            if(is.vector(x.init)){x.init <- matrix(x.init,length(x.init),1)} # If x.init is a vector transform to matrix
            if(dim(x.init)[2]==1){x.init <- matrix(rep(x.init,p.num),dimen,p.num)} # If x.init is a point then repeat
            p.mat <- apply(x.init,2,scale.transform); dimnames(p.mat)<-list(sprintf("dim.%i",1:dimen),NULL) # Transform x.init and define as p.mat
        } # Setting of p.mat depending on x.init NULL, d dimn vector, dxp dimn vector, 1xd matrix or pxd matrix
        ### Define storage matrices
        p.idx <- 1:p.num; log.p.wei <- rep(log(1/p.num),p.num) # Set idxs and weights
        p.layer <- matrix(0,9,p.num,dimnames=list(c("c.idx","s","deg.s","t","tau.bar","next.ev","PhiL","PhiU","Delta"),NULL)) # Establish particle layer matrix
        p.layer["c.idx",] <- p.idx # Set particle idxs in layer matrix
        p.dapts <- array(sample(dsz,3*p.num*ss.size,replace=TRUE),c(p.num,ss.size,3),dimnames=list(sprintf("part.%i",1:p.num),sprintf("ss.%i",1:ss.size),sprintf("dapt.%i",1:3))) # Set p.layer data points
        p.pass.arr <- p.cyc.arr <- array(0,c(dimen,6,p.num),dimnames=list(NULL,c("dimen","tau","y","minI","l","u"),sprintf("part.%i",1:p.num))) # Define particle passage and cyclic arrays
        for(i in 1:p.num){p.taus <- p.path.init(s=T.start,x=p.mat[,i],theta=theta,dimen=dimen); p.pass.arr[,,i] <- p.taus$pass.mat; p.cyc.arr[,,i] <- p.taus$cyc.mat} # Initialise particle first passage times
        p.layer["tau.bar",] <- p.cyc.arr[1,"tau",] # Determine each particles tau.bar
    }else{ #i.e. restarting algorithm from some point
        T.start <- prev.simn$T.fin; p.idx <- prev.simn$p.idx; log.p.wei <- prev.simn$log.p.wei; p.mat <- prev.simn$p.mat; p.layer <- prev.simn$p.layer; p.dapts <- prev.simn$p.dapts; p.pass.arr <- prev.simn$p.pass.arr; p.cyc.arr <- prev.simn$p.cyc.arr} # Restablish Algorithm
    ##### (0)1.3 Initialise particle system
    for(i in 1:p.num){p.phi.bds <- ss.phiC(p.mat[,i],p.pass.arr[,"l",i],p.pass.arr[,"u",i]); p.layer[c("PhiL","PhiU"),i]<-c(p.phi.bds$phiL,p.phi.bds$phiU)}; p.layer[c("Delta"),] <- p.layer["PhiU",] - p.layer["PhiL",] # Set p.layer current particle intensity
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
    anc.times <- curr.deg.time <- T.start; t.inc.next <- T.start + t.inc; resamp.counter <- 1 # Initialise current degradation time, time of first increment and resampling counter
    while(curr.deg.time < T.fin){ # Iterate until final time is reached
        ##### (t.inc)2.1 ##### Increment particle set to next event time
        ######################
        if(p.layer["next.ev",1] >= t.inc.next){ # If we have reached the next increment point
            ###### (t.inc)2.2.a.1 ##### Compute particles at intermediate time point
            p.layer <- mat.sort.r(p.layer,p.num,"c.idx") # Sort p.layer matrix by current particle index
            for(loc.up in 1:p.num){p.mat[,loc.up]<-p.path.renew(pass.mat=p.pass.arr[,,loc.up],cyc.mat=p.cyc.arr[,,loc.up],curr.time=p.layer["s",loc.up],x.curr=p.mat[,loc.up],next.time=t.inc.next,theta,dimen)$x.next}
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
            p.resamp <- scale_resamp(p.num,p.idx,log.p.wei,p.mat,p.layer,p.layer.sor.I=1,p.dapts,ss.size,resamp.method,ess.thresh,p.pass.arr=p.pass.arr,p.cyc.arr=p.cyc.arr) # Resample particles
            p.anc.ess <- c(p.anc.ess,p.resamp$ess) # Record ancestoral ESS at mesh point
            p.idx <- p.resamp$p.idx; log.p.wei <- p.resamp$log.p.wei # Update particle indicies and weights (normalised)
            if(p.resamp$resamp.I==1){p.anc.resamp <- c(p.anc.resamp,curr.deg.time); p.mat <- p.resamp$p.mat; p.layer <- p.resamp$p.layer; p.dapts <- p.resamp$p.dapts; p.pass.arr <- p.resamp$p.pass.arr; p.cyc.arr <- p.resamp$p.cyc.arr} # If there has been resampling, update ancestral resampling times, layer, dapts and cycling arrays
            p.layer <- mat.sort.r(p.layer,p.num,"next.ev") # Redefine storage vectors II
            resamp.counter <- 1 # Reset resample counter
            if(progress.check==TRUE){print(curr.deg.time)} # Output current time if checking progress
        }else{ # If we haven't reached the next increment point
            ###### (t.inc)2.2.b.1 ##### Select particle to be updated
            p.curr <- p.layer[,1] # Next particle to be updated
            p.curr.idx <- p.curr["c.idx"] # Define current index
            p.curr.dapts <- p.dapts[p.curr.idx,,,drop=FALSE] # Find particle subsampled data points
            ###### (t.inc)2.2.b.2 ##### Update particle trajectory
            traj.update <- p.path.renew(pass.mat=p.pass.arr[,,p.curr.idx],cyc.mat=p.cyc.arr[,,p.curr.idx],curr.time=p.curr["s"],x.curr=p.mat[,p.curr.idx],next.time=p.curr["next.ev"],theta,dimen) # Iterate passage matrix
            p.mat[,p.curr.idx] <- p.loc.next <- traj.update$x.next # Particle location at next.ev
            ###### (t.inc)2.2.b.3 ##### If particle to be updated has an event before the next tau.bar, then update to Poisson event time
            if(p.curr["t"]<=p.curr["tau.bar"]){
                ###### (t.inc)2.2.b.3.1 ##### Update data set
                p.dapts[p.curr.idx,,] <- sample(dsz,3*ss.size,replace=TRUE) # Update particle data points
                ###### (t.inc)2.2.b.3.2 ##### Update Poisson weight
                p.phi.eval <- ss.phi(p.loc.next,p.curr.dapts,ss.size)
                p.event.ch <- (p.curr["PhiU"] - p.phi.eval)/p.curr["Delta"]; if(phi.record==TRUE){p.phi.hist <- rbind(p.phi.hist,c(p.phi.eval,p.loc.next))} # Raw Poisson event increment
                if(p.event.ch < 0){p.neg <- neg.wei.mech(p.event.ch,p.curr.idx,p.curr,p.anc.neg,log.p.wei,ss.phiC,ss.phiC.up); log.p.wei[p.curr.idx] <- -Inf; p.anc.neg <- rbind(p.anc.neg,p.neg$p.anc.neg.append); ss.phiC <- p.neg$ss.phiC}else{log.p.wei[p.curr.idx] <- log.p.wei[p.curr.idx] + log(p.event.ch)} # Apply chosen negative weight mechanism
            }else{
                ###### (t.inc)2.2.b.4 ##### If particle to be updated has no event before the next tau.bar, then update passage trajectory
                p.pass.arr[,,p.curr.idx] <- traj.update$pass.mat; p.cyc.arr[,,p.curr.idx] <- traj.update$cyc.mat # Update pass.mat and cyc.mat
                p.curr["tau.bar"] <- traj.update$next.tau # Determine next particle passage time
            }
            ###### (t.inc)2.2.b.5 ##### Update particle degradation matrix
            p.bet.degrad[p.curr.idx] <- p.bet.degrad[p.curr.idx] + (p.curr["next.ev"]-p.curr["deg.s"])*p.curr[c("PhiL")] # Record individual particle degradation (since last resampling point)
            ###### (t.inc)2.2.b.6 ##### Update particle Poisson intensity
            p.phi.bds <- ss.phiC(p.loc.next,p.pass.arr[,"l",i],p.pass.arr[,"u",i]); p.curr[c("PhiL","PhiU")] <- c(p.phi.bds$phiL,p.phi.bds$phiU); p.curr["Delta"] <- p.curr["PhiU"] - p.curr["PhiL"] # Update phi bounds and intensity
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
                p.resamp <- scale_resamp(p.num,p.idx,log.p.wei,p.mat,p.layer,p.layer.sor.I=0,p.dapts,ss.size,resamp.method,ess.thresh,p.pass.arr=p.pass.arr,p.cyc.arr=p.cyc.arr) # Resample particles
                p.idx <- p.resamp$p.idx; log.p.wei <- p.resamp$log.p.wei # Update particle indicies and weights (normalised)
                if(p.resamp$resamp.I==1){p.anc.resamp <- c(p.anc.resamp,curr.deg.time); p.mat <- p.resamp$p.mat; p.layer <- p.resamp$p.layer; p.dapts <- p.resamp$p.dapts; p.pass.arr <- p.resamp$p.pass.arr; p.cyc.arr <- p.resamp$p.cyc.arr} # If there has been resampling, update ancestral resampling times, layer and dapts
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
    list(p.num=p.num,p.idx=p.idx,log.p.wei=log.p.wei,p.mat=p.mat,p.layer=p.layer,theta=theta,p.dapts=p.dapts,p.anc.idx=p.anc.idx,p.anc.wei=p.anc.wei,p.anc.mat=p.anc.mat,resamp.freq=resamp.freq,resamp.method=resamp.method,neg.wei.mech=neg.wei.mech,p.anc.resamp=p.anc.resamp,p.anc.ess=p.anc.ess,p.anc.neg=p.anc.neg,p.phi.hist=p.phi.hist,anc.times=anc.times,T.start=T.start,t.inc=t.inc,T.fin=T.fin,dimen=dimen,p.mu=p.mu,ss.size=ss.size,ess.thresh=ess.thresh,ss.phi=ss.phi,ss.phiC=ss.phiC,p.cyc.arr=p.cyc.arr,p.pass.arr=p.pass.arr,scale.transform=scale.transform,un.scale.transform=un.scale.transform,progress.check=progress.check,phi.record=phi.record,p.path.renew=p.path.renew)} # Output list


#############################################
#############################################
#### 5 - Scale Ergodic Extraction
#############################################
#############################################
scale_ergodic <- function(simn,retain.frac=NULL,retain.frac.range=NULL,trunc.time=NULL,trunc.time.range=NULL,symm=NULL,even.weights=NULL){
    #### 1.1 #### Assign frac.idxs depending on user assignation of retain frac of simulation
    if(is.null(retain.frac)==TRUE){
        frac.idxs <- 1:length(simn$anc.times) # If no specification include all
    }else{
        frac.idxs <- floor(length(simn$anc.times)*retain.frac):length(simn$anc.times)} # If specified include specified number
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
    times <- simn$anc.times[idxs] # Define times
    if(length(idxs)==0){weights <- numeric(0); draws <- draws <- matrix(0,0,simn$dimen,dimnames=list(NULL,sprintf("dim.%i",1:simn$dimen)))}else{
        #### 3 #### Extract particle weights
        simn.weis <- simn$p.anc.wei[idxs,,drop=FALSE] # Define subset of weights
        weights <- simn$p.anc.wei[r.idxs] # Assign weights
        weights <- weights / sum(weights) # Normalise weights
        #### 4 #### Extract particle locations location
        simn.draws <- simn$p.anc.mat[,,idxs,drop=FALSE] # Define subset of draws
        draws <- matrix(0,simn$p.num*length(idxs),simn$dimen,dimnames=list(NULL,sprintf("dim.%i",1:simn$dimen))) # Define draws matrix
        for(i in 1:simn$dimen){ # Assign draws dimension at a time
            simn.draws.dim <- simn.draws[i,,,drop=FALSE] # Extract dimension
            draws[,i] <- simn.draws.dim[r.idxs]} # Assign dimension to draw matrix
        #### 5 #### If symmetric, symmetrise
        if(is.null(symm)==FALSE){draws <- rbind(draws,-draws); weights <- c(weights/2,weights/2)}}
    if(is.null(even.weights)==FALSE){if(even.weights==TRUE){p.idxs <- sample(length(weights),length(weights),replace=TRUE,prob=weights); weights <- rep(1/length(weights),length(weights)); draws <- draws[p.idxs,,drop=FALSE]}}
    #### 6 #### Output ergodic average
    list(idxs=idxs,times=times,weights=weights,draws=draws)}

#############################################
#############################################
#### 6 - Extend existing scale algorithm
#############################################
#############################################
scale_extend <- function(prev.simn,t.inc,T.extend){
    ### Extend current simulation
    new.simn <- scale_exact(p.num=prev.simn$p.num,t.inc=t.inc,T.fin=T.extend+prev.simn$T.fin,ss.phi=ss.phi,ss.phiC=prev.simn$ss.phiC,dimen=prev.simn$dimen,scale.transform=scale.transform,un.scale.transform=un.scale.transform,T.start=prev.simn$T.fin,ss.size=prev.simn$ss.size,ess.thresh=prev.simn$ess.thresh,resamp.method=prev.simn$resamp.method,neg.wei.mech=prev.simn$neg.wei.mech,prev.simn=prev.simn,progress.check=prev.simn$progress.check,phi.record=prev.simn$phi.record,resamp.freq=prev.simn$resamp.freq,theta=prev.simn$theta,p.path.renew=prev.simn$p.path.renew)
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
    list(p.num=new.simn$p.num,p.idx=new.simn$p.idx,log.p.wei=new.simn$log.p.wei,p.mat=new.simn$p.mat,p.layer=new.simn$p.layer,theta=new.simn$theta,p.dapts=new.simn$p.dapts,p.anc.idx=p.anc.idx,p.anc.wei=p.anc.wei,p.anc.mat=p.anc.mat,resamp.freq=new.simn$resamp.freq,resamp.method=new.simn$resamp.method,neg.wei.mech=new.simn$neg.wei.mech,p.anc.resamp=p.anc.resamp,p.anc.ess=p.anc.ess,p.anc.neg=p.anc.neg,p.phi.hist=p.phi.hist,anc.times=anc.times,T.start=anc.times[1],t.inc=new.simn$t.inc,T.fin=new.simn$T.fin,dimen=new.simn$dimen,p.mu=new.simn$p.mu,ss.size=new.simn$ss.size,ess.thresh=new.simn$ess.thresh,ss.phi=new.simn$ss.phi,ss.phiC=new.simn$ss.phiC,p.cyc.arr=new.simn$p.cyc.arr,p.pass.arr=new.simn$p.pass.arr,scale.transform=prev.simn$scale.transform,un.scale.transform=prev.simn$un.scale.transform,progress.check=new.simn$progress.check,phi.record=new.simn$phi.record,p.path.renew=new.simn$p.path.renew)} # Output list



