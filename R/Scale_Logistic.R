############################################
############################################
######       Scale Algorithm         #######
######        Logistic Data          #######
######  Last Updated: 18/04/2018 MP  #######
############################################
############################################

### SCALE_LOGISTIC
scale_logistic <- function(fnm="logistic_default.RData",p.num=2^10,t.inc=0.1,T.extend=0.1,run.length=100000,ss.size=10,ss.on=TRUE,seed.default=1,default.data=TRUE){
    set.seed(seed.default)
    if(default.data==TRUE){big.logistic.example()}else{user.logistic.example()}
    time.elapse <- system.time(simn <<- scale_exact(p.num=p.num,t.inc=t.inc,T.fin=T.extend,ss.phi=ss.phi,ss.phiC=ss.phiC,dimen,scale.transform,un.scale.transform,T.start=0,x.init=NULL,ss.size=ss.size,ess.thresh=0.5,resamp.method=resid.resamp,neg.wei.mech=scale_zero.wei,prev.simn=NULL,progress.check=FALSE,phi.record=FALSE,resamp.freq=p.num-1,theta=NULL,p.path.renew=p.path.renew))[3]
    save(file=fnm,simn=simn,time.elapse=time.elapse); print(time.elapse)
    while(simn$T.fin < run.length){time.inc <- system.time(simn <<- scale_extend(simn,t.inc,T.extend))[3]
    time.elapse <- time.elapse + time.inc
        save(file=fnm,simn,time.elapse); print(time.elapse)}}

#p.num - function
#t.inc - function
#T.extend - function
#ss.phi - *
#ss.phiC - *
#dimen - *
#scale.transform - *
#un.scale.transform - *
#ss.size - function
#resid.resamp / resamp.method - SCALE FILE
#scale_zero.wei - SCALE FILE
#p.path.renew - SCALE FILE
#run.length - function

